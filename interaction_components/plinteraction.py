#!/usr/bin/env python
# coding: utf-8
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import Descriptors

from collections import namedtuple
from operator import itemgetter
from interaction_components.utils import *

from interaction_components.utils import centroid, tilde_expansion, tmpfile, classify_by_name, get_atom_coords
from interaction_components.utils import cluster_doubles, is_lig, normalize_vector, vector, ring_is_planar
from interaction_components.utils import extract_pdbid, read_pdb, create_folder_if_not_exists, canonicalize
from interaction_components.utils import read, nucleotide_linkage, sort_members_by_importance, is_acceptor, is_donor
from interaction_components.utils import whichchain, whichatomname, whichrestype, whichresnumber, euclidean3d, int32_to_negative
from interaction_components.detection import halogen, pication, water_bridges, metal_complexation
from interaction_components.detection import hydrophobic_interactions, pistacking, hbonds, saltbridge
from interaction_components import config


def get_features(mol):
    donors, acceptors, hydrophobics = [], [], []
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atomics = [a.GetAtomicNum() for a in atom.GetNeighbors()]
        if symbol in ["O", "N", "S"]:
            if atomics.count(1) >= 1:
                donors.append(atom.GetIdx())
            elif symbol in ["O", "S"] and atom.GetExplicitValence() <= 2:
                acceptors.append(atom.GetIdx())
            elif symbol == "N" and atomics.count(1) == 0 and atom.GetExplicitValence() <= 3:
                acceptors.append(atom.GetIdx())
        elif (atom.GetAtomicNum() == 6 and set(atomics).issubset({1, 6})):
            hydrophobics.append(atom.GetIdx())

    data = namedtuple("features", "donors acceptors hydrophobics")
    donors = list(set(donors))
    acceptors = list(set(acceptors))
    hydrophobics = list(set(hydrophobics))
    return data(donors=donors, acceptors=acceptors, hydrophobics=hydrophobics)


class Mol:
    def __init__(self, mol):
        self.mol = mol
        self.all_atoms = mol.GetAtoms()
        self.mol_conf = mol.GetConformers()[0]
        self.rings = None
        self.hydroph_atoms = None
        self.charged = None
        self.hbond_don_atom_pairs = None
        self.hbond_acc_atoms = None

    def find_hba(self):
        raise Exception("have to find hbond acceptors!")

    def find_hbd(self):
        raise Exception("have to find hbond donors!")

    def find_rings(self):
        """Find rings and return only aromatic.
        Rings have to be sufficiently planar OR be detected by OpenBabel as aromatic."""
        data = namedtuple(
            'aromatic_ring',
            'atoms orig_atoms atoms_orig_idx normal obj center type')
        rings = []
        aromatic_amino = ['TYR', 'TRP', 'HIS', 'PHE']

        ring_info = self.mol.GetRingInfo()

        rings_atom_idx = ring_info.AtomRings()
        for ring in rings_atom_idx:
            r_atoms = [self.mol.GetAtomWithIdx(idx) for idx in ring]
            r_atoms = sorted(r_atoms, key=lambda x: x.GetIdx())

            if 4 < len(ring) <= 6:
                res = list(set([whichrestype(a) for a in r_atoms]))
                if res[0] == "UNL":
                    ligand_orig_idx = ring
                    sort_order = np.argsort(np.array(ligand_orig_idx))
                    r_atoms = [r_atoms[i] for i in sort_order]
                if is_aromatic(r_atoms) or res[0] in aromatic_amino or ring_is_planar(
                        self.mol_conf, ring, r_atoms):
                    ring_type = '%s-membered' % len(ring)
                    ring_atms = get_coords(
                        self.mol_conf, [
                            r_atoms[0].GetIdx(), r_atoms[2].GetIdx(), r_atoms[4].GetIdx()])
                    ringv1 = vector(ring_atms[0], ring_atms[1])
                    ringv2 = vector(ring_atms[2], ring_atms[0])

                    atoms_orig_idx = [r_atom.GetIdx() for r_atom in r_atoms]
                    orig_atoms = r_atoms
                    rings.append(
                        data(
                            atoms=r_atoms,
                            orig_atoms=orig_atoms,
                            atoms_orig_idx=atoms_orig_idx,
                            normal=normalize_vector(
                                np.cross(
                                    ringv1,
                                    ringv2)),
                            obj=ring,
                            center=centroid(
                                get_coords(
                                    self.mol_conf,
                                    ring)),
                            type=ring_type))

        return rings

    def find_hal(self):
        """Look for halogen bond acceptors (Y-{O|P|N|S}, with Y=C,P,S)"""
        data = namedtuple(
            'hal_acceptor',
            'o o_orig_idx y y_orig_idx o_coords y_coords')
        a_set = []
        # All oxygens, nitrogen, sulfurs with neighboring carbon, phosphor,
        # nitrogen or sulfur
        for a in [
            at for at in self.mol.GetAtoms() if at.GetAtomicNum() in [
                8,
                7,
                16]]:
            n_atoms = [na for na in a.GetNeighbors() if na.GetAtomicNum() in [
                6, 7, 15, 16]]
            if len(n_atoms) == 1:  # Proximal atom
                o_orig_idx = a.GetIdx()
                y_orig_idx = n_atoms[0].GetIdx()
                o_coords = get_atom_coords(a)
                y_coords = get_atom_coords(n_atoms[0])
                a_set.append(data(o=a, o_orig_idx=o_orig_idx, y=n_atoms[0],
                                  y_orig_idx=y_orig_idx, o_coords=o_coords,
                                  y_coords=y_coords))
        return a_set

    def get_hydrophobic_atoms(self):
        return self.hydroph_atoms

    def get_hba(self):
        return self.hbond_acc_atoms

    def get_hbd(self):
        return [
            don_pair for don_pair in self.hbond_don_atom_pairs if don_pair.type == 'regular']

    def get_weak_hbd(self):
        return [
            don_pair for don_pair in self.hbond_don_atom_pairs if don_pair.type == 'weak']

    def get_pos_charged(self):
        return [charge for charge in self.charged if charge.type == 'positive']

    def get_neg_charged(self):
        return [charge for charge in self.charged if charge.type == 'negative']


class Ligand(Mol):
    def __init__(self, mol, mol_water):
        super(Ligand, self).__init__(mol)
        self.smiles = Chem.MolToSmiles(mol)
        self.heavy_atoms = mol.GetNumHeavyAtoms()  # Heavy atoms count
        self.features = get_features(mol)
        self.rings = self.find_rings()
        self.hydroph_atoms = self.hydrophobic_atoms()
        self.hbond_acc_atoms = self.find_hba()
        self.num_rings = len(self.rings)

        self.num_rot_bonds = Descriptors.NumRotatableBonds(mol)

        self.hbond_don_atom_pairs = self.find_hbd()

        self.charged = self.find_charged()
        self.centroid = np.round(
            self.mol_conf.GetPositions().mean(
                axis=0), 5).tolist()
        self.max_dist_to_center = max(
            (euclidean3d(
                self.centroid,
                get_coord(
                    self.mol_conf,
                    a.GetIdx())) for a in mol.GetAtoms()))

        self.water = []
        if mol_water is not None:
            self.mol_water = mol_water
            self.mol_water_conf = mol_water.GetConformers()[0]
            data = namedtuple('water', 'oxy oxy_orig_idx oxy_coords')
            for hoh in self.mol_water.GetAtoms():
                oxy = None
                if hoh.GetAtomicNum() == 8:
                    oxy = hoh
                if oxy is not None:
                    if euclidean3d(self.centroid, get_atom_coords(
                            oxy)) < self.max_dist_to_center + config.BS_DIST:
                        oxy_orig_idx = oxy.GetIdx()
                        oxy_coords = get_atom_coords(oxy)
                        self.water.append(
                            data(
                                oxy=oxy,
                                oxy_orig_idx=oxy_orig_idx,
                                oxy_coords=oxy_coords))

        self.halogenbond_don = self.find_hal()

        self.metal_binding = self.find_metal_binding()
        self.num_hba, self.num_hbd = len(
            self.hbond_acc_atoms), len(
            self.hbond_don_atom_pairs)
        self.num_hal = len(self.halogenbond_don)

    def hydrophobic_atoms(self):
        """Select all carbon atoms which have only carbons and/or hydrogens as direct neighbors."""
        atom_set = []
        data = namedtuple('hydrophobic', 'atom orig_atom orig_idx coords')
        atom_idx_set = self.features.hydrophobics
        atm = [self.mol.GetAtomWithIdx(idx) for idx in atom_idx_set]

        for atom in atm:
            orig_idx = atom.GetIdx()
            orig_atom = orig_idx
            atom_set.append(
                data(
                    atom=atom,
                    orig_atom=orig_atom,
                    orig_idx=orig_idx,
                    coords=get_atom_coords(atom)))
        return atom_set

    def find_hba(self):
        """Find all possible hydrogen bond acceptors"""
        data = namedtuple(
            'hbondacceptor',
            'a a_orig_atom a_orig_idx type coords')
        a_set = []

        atom_idx_set = self.features.acceptors
        atm = [self.mol.GetAtomWithIdx(idx) for idx in atom_idx_set]

        for atom_idx, atom in zip(atom_idx_set, atm):
            if atom.GetAtomicNum() not in [
                    9, 17, 35, 53]:  # Exclude halogen atoms
                a_orig_idx = atom_idx
                a_orig_atom = atom
                coords = get_atom_coords(atom)
                a_set.append(
                    data(
                        a=atom,
                        a_orig_atom=a_orig_atom,
                        a_orig_idx=a_orig_idx,
                        type='regular',
                        coords=coords))
        a_set = sorted(a_set, key=lambda x: x.a_orig_idx)
        return a_set

    def find_hbd(self):
        """Find all possible strong and weak hydrogen bonds donors (all hydrophobic C-H pairings)"""
        donor_pairs = []
        data = namedtuple(
            'hbonddonor',
            'd d_orig_atom d_orig_idx h type d_coords h_coords')

        donor_idxs = self.features.donors
        donor_atoms = [self.mol.GetAtomWithIdx(idx) for idx in donor_idxs]

        for donor_idx, donor_atom in zip(donor_idxs, donor_atoms):
            in_ring = False
            if not in_ring:
                for adj_atom in [
                        a for a in donor_atom.GetNeighbors() if a.GetAtomicNum() == 1]:
                    d_orig_idx = donor_idx
                    d_orig_atom = donor_atom
                    d_coords = get_atom_coords(donor_atom)
                    h_coords = get_atom_coords(adj_atom)
                    donor_pairs.append(
                        data(
                            d=donor_atom,
                            d_orig_atom=d_orig_atom,
                            d_orig_idx=d_orig_idx,
                            h=adj_atom,
                            type='regular',
                            d_coords=d_coords,
                            h_coords=h_coords))

        for carbon in self.hydroph_atoms:
            for adj_atom in [
                    a for a in carbon.atom.GetNeighbors() if a.GetAtomicNum() == 1]:
                d_orig_idx = carbon.atom.GetIdx()
                d_orig_atom = carbon.atom
                d_coords = get_atom_coords(carbon.atom)
                h_coords = get_atom_coords(adj_atom)
                donor_pairs.append(
                    data(
                        d=carbon,
                        d_orig_atom=d_orig_atom,
                        d_orig_idx=d_orig_idx,
                        h=adj_atom,
                        type='weak',
                        d_coords=d_coords,
                        h_coords=h_coords))
        donor_pairs = sorted(
            donor_pairs, key=lambda x: (
                x.d_orig_idx, x.h.GetIdx()))
        return donor_pairs

    def find_hal(self):
        """Look for halogen bond donors (X-C, with X=F, Cl, Br, I)"""
        data = namedtuple(
            'hal_donor',
            'x orig_x x_orig_idx c c_orig_idx x_coords c_coords')
        a_set = []
        for a in self.all_atoms:
            if self.is_functional_group(a, 'halocarbon'):
                n_atoms = [
                    na for na in a.GetNeighbors() if na.GetAtomicNum() == 6]
                x_orig_idx = a.GetIdx()
                orig_x = x_orig_idx
                c_orig_idx = [na.GetIdx() for na in n_atoms]
                x_coords = get_atom_coords(a)
                c_coords = get_atom_coords(n_atoms[0])
                a_set.append(data(x=a, orig_x=orig_x, x_orig_idx=x_orig_idx,
                                  c=n_atoms[0], c_orig_idx=c_orig_idx,
                                  x_coords=x_coords, c_coords=c_coords))
        if len(a_set) != 0:
            #print(f'ligand contains {len(a_set)} halogen atom(s)')
            pass
        return a_set

    def find_charged(self):
        """Identify all positively charged groups in a ligand. This search is not exhaustive, as the cases can be quite
        diverse. The typical cases seem to be protonated amines, quaternary ammoinium and sulfonium
        as mentioned in 'Cation-pi interactions in ligand recognition and catalysis' (Zacharias et al., 2002)).
        Identify negatively charged groups in the ligand.
        """
        data = namedtuple(
            'lcharge',
            'atoms orig_atoms atoms_orig_idx type center fgroup')
        a_set = []
        for a in self.mol.GetAtoms():
            a_orig_idx = a.GetIdx()
            a_orig = a.GetIdx()
            if self.is_functional_group(a, 'quartamine'):
                a_set.append(
                    data(
                        atoms=[
                            a, ], orig_atoms=[
                            a_orig, ], atoms_orig_idx=[
                            a_orig_idx, ], type='positive', center=list(
                            get_coord(
                                self.mol_conf, a.GetIdx())), fgroup='quartamine'))
            elif self.is_functional_group(a, 'tertamine'):
                a_set.append(
                    data(
                        atoms=[
                            a, ], orig_atoms=[
                            a_orig, ], atoms_orig_idx=[
                            a_orig_idx, ], type='positive', center=list(
                            get_coord(
                                self.mol_conf, a.GetIdx())), fgroup='tertamine'))
            if self.is_functional_group(a, 'sulfonium'):
                a_set.append(
                    data(
                        atoms=[
                            a, ], orig_atoms=[
                            a_orig, ], atoms_orig_idx=[
                            a_orig_idx, ], type='positive', center=list(
                            get_coord(
                                self.mol_conf, a.GetIdx())), fgroup='sulfonium'))
            if self.is_functional_group(a, 'phosphate'):
                a_contributing = [a, ]
                a_contributing_orig_idx = [a_orig_idx, ]
                [a_contributing.append(neighbor)
                 for neighbor in a.GetNeighbors()]
                [a_contributing_orig_idx.append(neighbor.GetIdx())
                 for neighbor in a_contributing]
                orig_contributing = [idx for idx in a_contributing_orig_idx]
                a_set.append(
                    data(
                        atoms=a_contributing,
                        orig_atoms=orig_contributing,
                        atoms_orig_idx=a_contributing_orig_idx,
                        type='negative',
                        center=list(
                            get_coord(
                                self.mol_conf,
                                a.GetIdx())),
                        fgroup='phosphate'))
            if self.is_functional_group(a, 'sulfonicacid'):
                a_contributing = [a, ]
                a_contributing_orig_idx = [a_orig_idx, ]
                [a_contributing.append(neighbor) for neighbor in a.GetNeighbors(
                ) if neighbor.GetAtomicNum() == 8]
                [a_contributing_orig_idx.append(
                    neighbor.GetIdx()) for neighbor in a_contributing]
                orig_contributing = a_contributing_orig_idx
                a_set.append(
                    data(
                        atoms=a_contributing,
                        orig_atoms=orig_contributing,
                        atoms_orig_idx=a_contributing_orig_idx,
                        type='negative',
                        center=list(
                            get_coord(
                                self.mol_conf,
                                a.GetIdx())),
                        fgroup='sulfonicacid'))
            elif self.is_functional_group(a, 'sulfate'):
                a_contributing = [a, ]
                a_contributing_orig_idx = [a_orig_idx, ]
                [a_contributing_orig_idx.append(
                    neighbor.GetIdx()) for neighbor in a_contributing]
                [a_contributing.append(neighbor)
                 for neighbor in a.GetNeighbors()]
                orig_contributing = a_contributing_orig_idx
                a_set.append(
                    data(
                        atoms=a_contributing,
                        orig_atoms=orig_contributing,
                        atoms_orig_idx=a_contributing_orig_idx,
                        type='negative',
                        center=get_coord(
                            self.mol_conf,
                            a.GetIdx()),
                        fgroup='sulfate'))
            if self.is_functional_group(a, 'carboxylate'):
                a_contributing = [
                    neighbor for neighbor in a.GetNeighbors() if neighbor.GetAtomicNum() == 8]
                a_contributing_orig_idx = [
                    neighbor.GetIdx() for neighbor in a_contributing]
                orig_contributing = a_contributing_orig_idx
                a_set.append(data(atoms=a_contributing,
                                  orig_atoms=orig_contributing,
                                  atoms_orig_idx=a_contributing_orig_idx,
                                  type='negative',
                                  center=centroid([get_coord(self.mol_conf,
                                                             a.GetIdx()) for a in a_contributing]),
                                  fgroup='carboxylate'))
            elif self.is_functional_group(a, 'guanidine'):
                a_contributing = [
                    neighbor for neighbor in a.GetNeighbors() if neighbor.GetAtomicNum() == 7]
                a_contributing_orig_idx = [
                    neighbor.GetIdx() for neighbor in a_contributing]
                orig_contributing = a_contributing_orig_idx
                a_set.append(
                    data(
                        atoms=a_contributing,
                        orig_atoms=orig_contributing,
                        atoms_orig_idx=a_contributing_orig_idx,
                        type='positive',
                        center=get_coord(
                            self.mol_conf,
                            a.GetIdx()),
                        fgroup='guanidine'))
        return a_set

    def is_functional_group(self, atom, group):
        """Given a pybel atom, look up if it belongs to a function group"""
        n_atoms = [a_neighbor.GetAtomicNum()
                   for a_neighbor in atom.GetNeighbors()]

        if group in [
            'quartamine',
                'tertamine'] and atom.GetAtomicNum() == 7:  # Nitrogen
            # It's a nitrogen, so could be a protonated amine or quaternary
            # ammonium
            if '1' not in n_atoms and len(n_atoms) == 4:
                # It's a quat. ammonium (N with 4 residues != H)
                return True if group == 'quartamine' else False
            elif atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and len(n_atoms) >= 3:
                # It's sp3-hybridized, so could pick up an hydrogen
                return True if group == 'tertamine' else False
            else:
                return False

        if group in ['sulfonium', 'sulfonicacid',
                     'sulfate'] and atom.GetAtomicNum() == 16:  # Sulfur
            if '1' not in n_atoms and len(
                    n_atoms) == 3:  # It's a sulfonium (S with 3 residues != H)
                return True if group == 'sulfonium' else False
            elif n_atoms.count(8) == 3:  # It's a sulfonate or sulfonic acid
                return True if group == 'sulfonicacid' else False
            elif n_atoms.count(8) == 4:  # It's a sulfate
                return True if group == 'sulfate' else False

        if group == 'phosphate' and atom.GetAtomicNum() == 15:  # Phosphor
            if set(n_atoms) == {8}:  # It's a phosphate
                return True

        if group in [
            'carboxylate',
                'guanidine'] and atom.GetAtomicNum() == 6:  # It's a carbon atom
            if n_atoms.count(8) == 2 and n_atoms.count(
                    6) == 1:  # It's a carboxylate group
                return True if group == 'carboxylate' else False
            elif n_atoms.count(7) == 3 and len(n_atoms) == 3:  # It's a guanidine group
                nitro_partners = []
                for nitro in atom.GetNeighbors():
                    nitro_partners.append(
                        len([b_neighbor for b_neighbor in nitro.GetNeighbors()]))
                if min(
                        nitro_partners) == 1:  # One nitrogen is only connected to the carbon, can pick up a H
                    return True if group == 'guanidine' else False

        if group == 'halocarbon' and atom.GetAtomicNum() in [
                9, 17, 35, 53]:  # Halogen atoms
            n_atoms = [
                na for na in atom.GetNeighbors() if na.GetAtomicNum() == 6]
            if len(n_atoms) == 1:  # Halocarbon
                return True
        else:
            return False

    def find_metal_binding(self):
        """Looks for atoms that could possibly be involved in binding a metal ion.
        This can be any water oxygen, as well as oxygen from carboxylate, phophoryl, phenolate, alcohol;
        nitrogen from imidazole; sulfur from thiolate.
        """
        a_set = []
        data = namedtuple(
            'metal_binding',
            'atom orig_atom atom_orig_idx type fgroup restype resnr reschain location coords')
        for oxygen in self.water:
            a_set.append(
                data(
                    atom=oxygen.oxy,
                    atom_orig_idx=oxygen.oxy_orig_idx,
                    type='O',
                    fgroup='water',
                    restype="HOH",
                    resnr=oxygen.oxy.GetPDBResidueInfo().GetResidueNumber(),
                    reschain="W",
                    coords=oxygen.oxy_coords,
                    location='water',
                    orig_atom=oxygen.oxy_orig_idx))

        for a in self.mol.GetAtoms():
            a_orig_idx = a.GetIdx()
            n_atoms = a.GetNeighbors()
            # All atomic numbers of neighboring atoms
            n_atoms_atomicnum = [n.GetAtomicNum() for n in a.GetNeighbors()]
            if a.GetAtomicNum() == 8:  # Oxygen
                if n_atoms_atomicnum.count('1') == 1 and len(
                        n_atoms_atomicnum) == 2:  # Oxygen in alcohol (R-[O]-H)
                    a_set.append(
                        data(
                            atom=a,
                            atom_orig_idx=a_orig_idx,
                            type='O',
                            coords=get_atom_coords(a),
                            fgroup='alcohol',
                            restype="l",
                            resnr=1,
                            reschain="L",
                            location='ligand',
                            orig_atom=a_orig_idx))
                if True in [
                        n.GetIsAromatic() for n in n_atoms] and not a.GetIsAromatic():  # Phenolate oxygen
                    a_set.append(
                        data(
                            atom=a,
                            atom_orig_idx=a_orig_idx,
                            type='O',
                            coords=get_atom_coords(a),
                            fgroup='phenolate',
                            restype="l",
                            resnr=1,
                            reschain="L",
                            location='ligand',
                            orig_atom=a_orig_idx))
            if a.GetAtomicNum() == 6:  # It's a carbon atom
                if n_atoms_atomicnum.count(8) == 2 and n_atoms_atomicnum.count(
                        6) == 1:  # It's a carboxylate group
                    for neighbor in [
                            n for n in n_atoms if n.GetAtomicNum() == 8]:
                        neighbor_orig_idx = neighbor.GetIdx()
                        a_set.append(
                            data(
                                atom=neighbor,
                                atom_orig_idx=neighbor_orig_idx,
                                type='O',
                                fgroup='carboxylate',
                                restype="l",
                                resnr=1,
                                reschain="L",
                                location='ligand',
                                orig_atom=a_orig_idx,
                                coords=get_atom_coords(neighbor)))
            if a.GetAtomicNum() == 15:  # It's a phosphor atom
                if n_atoms_atomicnum.count(8) >= 3:  # It's a phosphoryl
                    for neighbor in [
                            n for n in n_atoms if n.GetAtomicNum() == 8]:
                        neighbor_orig_idx = neighbor.GetIdx()
                        a_set.append(
                            data(
                                atom=neighbor,
                                atom_orig_idx=neighbor_orig_idx,
                                type='O',
                                fgroup='phosphoryl',
                                restype="l",
                                resnr=1,
                                reschain="L",
                                location='ligand',
                                orig_atom=a_orig_idx,
                                coords=get_atom_coords(neighbor)))
                if n_atoms_atomicnum.count(
                        8) == 2:  # It's another phosphor-containing group #@todo (correct name?)
                    for neighbor in [
                            n for n in n_atoms if n.GetAtomicNum() == 8]:
                        neighbor_orig_idx = neighbor.GetIdx()
                        a_set.append(
                            data(
                                atom=neighbor,
                                atom_orig_idx=neighbor_orig_idx,
                                type='O',
                                fgroup='phosphor.other',
                                restype="l",
                                resnr=1,
                                reschain="L",
                                location='ligand',
                                orig_atom=a_orig_idx,
                                coords=get_atom_coords(neighbor)))
            if a.GetAtomicNum() == 7:  # It's a nitrogen atom
                if n_atoms_atomicnum.count(
                        6) == 2:  # It's imidazole/pyrrole or similar
                    a_set.append(
                        data(
                            atom=a,
                            atom_orig_idx=a_orig_idx,
                            type='N',
                            coords=get_atom_coords(a),
                            fgroup='imidazole/pyrrole',
                            restype="l",
                            resnr=1,
                            reschain="L",
                            location='ligand',
                            orig_atom=a_orig_idx))
            if a.GetAtomicNum() == 16:  # It's a sulfur atom
                if True in [n.GetIsAromatic()
                            for n in n_atoms] and not a.GetIsAromatic():  # Thiolate
                    a_set.append(
                        data(
                            atom=a,
                            atom_orig_idx=a_orig_idx,
                            type='S',
                            coords=get_atom_coords(a),
                            fgroup='thiolate',
                            restype="l",
                            resnr=1,
                            reschain="L",
                            location='ligand',
                            orig_atom=a_orig_idx))
                if set(n_atoms_atomicnum) == {
                        26}:  # Sulfur in Iron sulfur cluster
                    a_set.append(
                        data(
                            atom=a,
                            atom_orig_idx=a_orig_idx,
                            type='S',
                            coords=get_atom_coords(a),
                            fgroup='iron-sulfur.cluster',
                            restype="l",
                            resnr=1,
                            reschain="L",
                            location='ligand',
                            orig_atom=a_orig_idx))
        return a_set


class Protein(Mol):
    def __init__(self, mol):
        super(Protein, self).__init__(mol)
        self.rings = self.find_rings()
        self.hydroph_atoms = self.hydrophobic_atoms()
        self.hbond_acc_atoms = self.find_hba()
        self.hbond_don_atom_pairs = self.find_hbd()
        self.residues = residue_order(mol)

        self.charged = self.find_charged()
        self.halogenbond_acc = self.find_hal()
        self.metal_binding = self.find_metal_binding()

        self.atom_prop_dict = config.atom_prop_dict

        self.metals = []
        data = namedtuple('metal', 'm orig_m m_orig_idx m_coords')
        for a in [a for a in self.all_atoms if a.GetSymbol().upper()
                  in config.METAL_IONS]:
            m_orig_idx = a.GetIdx()
            orig_m = m_orig_idx
            self.metals.append(
                data(
                    m=a,
                    m_orig_idx=m_orig_idx,
                    orig_m=orig_m,
                    m_coords=get_atom_coords(a)))

    def hydrophobic_atoms(self):
        """Select all carbon atoms which have only carbons and/or hydrogens as direct neighbors."""
        atom_set = []
        data = namedtuple('hydrophobic', 'atom orig_atom orig_idx coords')
        atm = [a for a in self.all_atoms if (a.GetAtomicNum() == 6 and set(
            [natom.GetAtomicNum() for natom in a.GetNeighbors()]).issubset({1, 6}))]
        for atom in atm:
            orig_idx = atom.GetIdx()
            orig_atom = orig_idx
            atom_set.append(
                data(
                    atom=atom,
                    orig_atom=orig_atom,
                    orig_idx=orig_idx,
                    coords=get_atom_coords(atom)))
        return atom_set

    def find_hba(self):
        data = namedtuple(
            'hbondacceptor',
            'a a_orig_atom a_orig_idx type coords')
        a_set = []
        for atom in self.all_atoms:
            if is_acceptor(atom):
                a_orig_idx = atom.GetIdx()
                a_orig_atom = atom
                a_set.append(
                    data(
                        a=atom,
                        a_orig_atom=a_orig_atom,
                        a_orig_idx=a_orig_idx,
                        type='regular',
                        coords=get_atom_coords(atom)))
        a_set = sorted(a_set, key=lambda x: x.a_orig_idx)
        return a_set

    def find_hbd(self):
        donor_pairs = []
        data = namedtuple(
            'hbonddonor',
            'd d_orig_atom d_orig_idx h type d_coords h_coords')
        for donor in [a for a in self.all_atoms if is_donor(a)]:
            in_ring = False
            if not in_ring:
                for adj_atom in [
                        a for a in donor.GetNeighbors() if a.GetAtomicNum() == 1]:
                    d_orig_idx = donor.GetIdx()
                    d_orig_atom = donor
                    d_coords = get_atom_coords(donor)
                    h_coords = get_atom_coords(adj_atom)
                    donor_pairs.append(
                        data(
                            d=donor,
                            d_orig_atom=d_orig_atom,
                            d_orig_idx=d_orig_idx,
                            h=adj_atom,
                            type='regular',
                            d_coords=d_coords,
                            h_coords=h_coords))

        for carbon in self.hydroph_atoms:
            for adj_atom in [
                    a for a in carbon.atom.GetNeighbors() if a.GetAtomicNum() == 1]:
                d_orig_idx = carbon.atom.GetIdx()
                d_orig_atom = carbon.atom
                d_coords = get_atom_coords(carbon.atom)
                h_coords = get_atom_coords(adj_atom)
                donor_pairs.append(
                    data(
                        d=carbon,
                        d_orig_atom=d_orig_atom,
                        d_coords=d_coords,
                        h_coords=h_coords,
                        d_orig_idx=d_orig_idx,
                        h=adj_atom,
                        type='weak'))
        donor_pairs = sorted(
            donor_pairs, key=lambda x: (
                x.d_orig_idx, x.h.GetIdx()))
        return donor_pairs

    def find_hal(self):
        """Look for halogen bond acceptors (Y-{O|P|N|S}, with Y=C,P,S)"""
        data = namedtuple(
            'hal_acceptor',
            'o o_orig_idx y y_orig_idx o_coords y_coords')
        a_set = []
        # All oxygens, nitrogen, sulfurs with neighboring carbon, phosphor,
        # nitrogen or sulfur
        for a in [at for at in self.all_atoms if at.GetAtomicNum() in [
                8, 7, 16]]:
            n_atoms = [na for na in a.GetNeighbors() if na.GetAtomicNum() in [
                6, 7, 15, 16]]
            if len(n_atoms) == 1:  # Proximal atom
                o_orig_idx = a.GetIdx()
                y_orig_idx = n_atoms[0].GetIdx()
                o_coords = get_atom_coords(a)
                y_coords = get_atom_coords(n_atoms[0])
                a_set.append(data(o=a, o_orig_idx=o_orig_idx, y=n_atoms[0],
                                  y_orig_idx=y_orig_idx, o_coords=o_coords,
                                  y_coords=y_coords))
        return a_set

    def find_charged(self):
        """Looks for positive charges in arginine, histidine or lysine, for negative in aspartic and glutamic acid."""
        data = namedtuple(
            'pcharge',
            'atoms atoms_orig_idx type center restype resnr reschain')
        a_set = []
        # Iterate through all residue, exclude those in chains defined as
        # peptides

        for res in self.residues:
            a_contributing = []
            a_contributing_orig_idx = []
            if res.residue_name in ("ARG", "HIS", "LYS"):
                for a in res.residue_atoms:
                    if a.GetSymbol() == "N" and a.GetPDBResidueInfo().GetName().strip(" ") != "N":
                        a_contributing.append(a)
                        a_contributing_orig_idx.append(a.GetIdx())
                if not len(a_contributing) == 0:
                    a_set.append(
                        data(
                            atoms=a_contributing,
                            atoms_orig_idx=a_contributing_orig_idx,
                            type='positive',
                            center=centroid(
                                [
                                    get_atom_coords(ac) for ac in a_contributing]),
                            restype=res.residue_name,
                            resnr=res.residue_number,
                            reschain=res.residue_chain))
            if res.residue_name in ("GLU", "ASP"):
                for a in res.residue_atoms:
                    if a.GetSymbol() == "O" and a.GetPDBResidueInfo().GetName().strip(" ") != "O":
                        a_contributing.append(a)
                        a_contributing_orig_idx.append(a.GetIdx())
                if not len(a_contributing) == 0:
                    a_set.append(
                        data(
                            atoms=a_contributing,
                            atoms_orig_idx=a_contributing_orig_idx,
                            type='negative',
                            center=centroid(
                                [
                                    get_atom_coords(ac) for ac in a_contributing]),
                            restype=res.residue_name,
                            resnr=res.residue_number,
                            reschain=res.residue_chain))
        return a_set

    def find_metal_binding(self):
        """Looks for atoms that could possibly be involved in chelating a metal ion.
        This can be any main chain oxygen atom or oxygen, nitrogen and sulfur from specific amino acids"""
        data = namedtuple(
            'metal_binding',
            'atom atom_orig_idx type restype resnr reschain location coords')
        a_set = []
        for res in self.residues:
            restype, resnr = res.residue_name, res.residue_number
            reschain = 'P'
            if restype in ("ASP", "GLU", "SER", "THR", "TYR"):
                for a in res.residue_atoms:
                    if a.GetSymbol() == "O" and a.GetPDBResidueInfo().GetName().strip(" ") != "O":
                        atom_orig_idx = a.GetIdx()
                        a_set.append(
                            data(
                                atom=a,
                                atom_orig_idx=atom_orig_idx,
                                type='O',
                                restype=restype,
                                resnr=resnr,
                                reschain=reschain,
                                coords=get_atom_coords(a),
                                location='protein.sidechain'))
            if restype == 'HIS':  # Look for nitrogen here
                for a in res.residue_atoms:
                    if a.GetSymbol() == "N" and a.GetPDBResidueInfo().GetName().strip(" ") != "N":
                        atom_orig_idx = a.GetIdx()
                        a_set.append(
                            data(
                                atom=a,
                                atom_orig_idx=atom_orig_idx,
                                type='N',
                                restype=restype,
                                resnr=resnr,
                                reschain=reschain,
                                coords=get_atom_coords(a),
                                location='protein.sidechain'))
            if restype == 'CYS':  # Look for sulfur here
                for a in res.residue_atoms:
                    if a.GetSymbol() == "S":
                        atom_orig_idx = a.GetIdx()
                        a_set.append(
                            data(
                                atom=a,
                                atom_orig_idx=atom_orig_idx,
                                type='S',
                                restype=restype,
                                resnr=resnr,
                                reschain=reschain,
                                coords=get_atom_coords(a),
                                location='protein.sidechain'))

            for a in res.residue_atoms:  # All main chain oxygens
                if a.GetSymbol() == "O" and a.GetPDBResidueInfo().GetName().strip(" ") == "O":
                    atom_orig_idx = a.GetIdx()
                    a_set.append(
                        data(
                            atom=a,
                            atom_orig_idx=atom_orig_idx,
                            type='O',
                            restype=res.residue_name,
                            resnr=res.residue_number,
                            reschain=reschain,
                            coords=get_atom_coords(a),
                            location='protein.mainchain'))
        return a_set


class PLInteraction:
    """Class to store a ligand, a protein and their interactions."""

    def __init__(self, lig_obj, bs_obj, pdbid):
        """Detect all interactions when initializing"""
        self.ligand = lig_obj
        self.pdbid = pdbid
        self.protein = bs_obj

        # #@todo Refactor code to combine different directionality

        self.saltbridge_lneg = saltbridge(
            self.protein.get_pos_charged(),
            self.ligand.get_neg_charged(),
            True)
        self.saltbridge_pneg = saltbridge(
            self.ligand.get_pos_charged(),
            self.protein.get_neg_charged(),
            False)

        self.all_hbonds_ldon = hbonds(self.protein.get_hba(),
                                      self.ligand.get_hbd(), False, 'strong')
        self.all_hbonds_pdon = hbonds(self.ligand.get_hba(),
                                      self.protein.get_hbd(), True, 'strong')

        self.hbonds_ldon = self.refine_hbonds_ldon(
            self.all_hbonds_ldon, self.saltbridge_lneg, self.saltbridge_pneg)
        self.hbonds_pdon = self.refine_hbonds_pdon(
            self.all_hbonds_pdon, self.saltbridge_lneg, self.saltbridge_pneg)

        self.pistacking = pistacking(self.protein.rings, self.ligand.rings)

        self.all_pi_cation_laro = pication(
            self.ligand.rings, self.protein.get_pos_charged(), True)
        self.pication_paro = pication(
            self.protein.rings,
            self.ligand.get_pos_charged(),
            False)

        self.pication_laro = self.refine_pi_cation_laro(
            self.all_pi_cation_laro, self.pistacking)

        self.all_hydrophobic_contacts = hydrophobic_interactions(
            self.protein.get_hydrophobic_atoms(), self.ligand.get_hydrophobic_atoms())
        self.hydrophobic_contacts = self.refine_hydrophobic(
            self.all_hydrophobic_contacts, self.pistacking)
        self.halogen_bonds = halogen(
            self.protein.halogenbond_acc,
            self.ligand.halogenbond_don)
        self.all_water_bridges = water_bridges(
            self.protein.get_hba(),
            self.ligand.get_hba(),
            self.protein.get_hbd(),
            self.ligand.get_hbd(),
            self.ligand.water)

        self.water_bridges = self.refine_water_bridges(
            self.all_water_bridges, self.hbonds_ldon, self.hbonds_pdon)

        self.metal_complexes = metal_complexation(
            self.protein.metals,
            self.ligand.metal_binding,
            self.protein.metal_binding)

        self.all_itypes = self.saltbridge_lneg + \
            self.saltbridge_pneg + self.hbonds_pdon
        self.all_itypes = self.all_itypes + self.hbonds_ldon + \
            self.pistacking + self.pication_laro + self.pication_paro
        self.all_itypes = self.all_itypes + self.hydrophobic_contacts + \
            self.halogen_bonds + self.water_bridges
        self.all_itypes = self.all_itypes + self.metal_complexes

        self.no_interactions = all(len(i) == 0 for i in self.all_itypes)
        self.unpaired_hba, self.unpaired_hbd, self.unpaired_hal = self.find_unpaired_ligand()
        self.unpaired_hba_orig_idx = [atom.GetIdx()
                                      for atom in self.unpaired_hba]
        self.unpaired_hbd_orig_idx = [atom.GetIdx()
                                      for atom in self.unpaired_hbd]
        self.unpaired_hal_orig_idx = [atom.GetIdx()
                                      for atom in self.unpaired_hal]
        self.num_unpaired_hba, self.num_unpaired_hbd = len(
            self.unpaired_hba), len(self.unpaired_hbd)
        self.num_unpaired_hal = len(self.unpaired_hal)

        # Exclude empty chains (coming from ligand as a target, from metal
        # complexes)
        self.interacting_chains = sorted(list(set([i.reschain for i in self.all_itypes
                                                   if i.reschain not in [' ', None]])))

        # Get all interacting residues, excluding ligand and water molecules
        self.interacting_res = list(set([''.join([str(i.resnr), str(
            i.reschain)]) for i in self.all_itypes if i.restype not in ['LIG', 'HOH']]))
        if len(self.interacting_res) != 0:
            interactions_list = []
            num_saltbridges = len(self.saltbridge_lneg + self.saltbridge_pneg)
            num_hbonds = len(self.hbonds_ldon + self.hbonds_pdon)
            num_pication = len(self.pication_laro + self.pication_paro)
            num_pistack = len(self.pistacking)
            num_halogen = len(self.halogen_bonds)
            num_waterbridges = len(self.water_bridges)
            if num_saltbridges != 0:
                interactions_list.append('%i salt bridge(s)' % num_saltbridges)
            if num_hbonds != 0:
                interactions_list.append('%i hydrogen bond(s)' % num_hbonds)
            if num_pication != 0:
                interactions_list.append(
                    '%i pi-cation interaction(s)' %
                    num_pication)
            if num_pistack != 0:
                interactions_list.append('%i pi-stacking(s)' % num_pistack)
            if num_halogen != 0:
                interactions_list.append('%i halogen bond(s)' % num_halogen)
            if num_waterbridges != 0:
                interactions_list.append(
                    '%i water bridge(s)' %
                    num_waterbridges)
            if not len(interactions_list) == 0:
                # raise RuntimeWarning(f'complex uses {interactions_list}')
                #print(f'complex uses {interactions_list}')
                pass
        else:
            # raise RuntimeWarning('no interactions for this ligand')
            print('no interactions for this ligand')

    def find_unpaired_ligand(self):
        """Identify unpaired functional in groups in ligands, involving H-Bond donors, acceptors, halogen bond donors.
        """
        unpaired_hba, unpaired_hbd, unpaired_hal = [], [], []
        # Unpaired hydrogen bond acceptors/donors in ligand (not used for
        # hydrogen bonds/water, salt bridges/mcomplex)
        involved_atoms = [hbond.a.GetIdx() for hbond in self.hbonds_pdon] + \
            [hbond.d.GetIdx() for hbond in self.hbonds_ldon]
        [[involved_atoms.append(atom.GetIdx()) for atom in sb.negative.atoms]
         for sb in self.saltbridge_lneg]
        [[involved_atoms.append(atom.GetIdx()) for atom in sb.positive.atoms]
         for sb in self.saltbridge_pneg]
        [involved_atoms.append(wb.a.GetIdx())
         for wb in self.water_bridges if wb.protisdon]
        [involved_atoms.append(wb.d.GetIdx())
         for wb in self.water_bridges if not wb.protisdon]
        [involved_atoms.append(mcomplex.target.atom.GetIdx(
        )) for mcomplex in self.metal_complexes if mcomplex.location == 'ligand']
        for atom in [hba.a for hba in self.ligand.get_hba()]:
            if atom.GetIdx() not in involved_atoms:
                unpaired_hba.append(atom)
        for atom in [hbd.d for hbd in self.ligand.get_hbd()]:
            if atom.GetIdx() not in involved_atoms:
                unpaired_hbd.append(atom)

        # unpaired halogen bond donors in ligand (not used for the previous +
        # halogen bonds)
        [involved_atoms.append(atom.don.x.GetIdx())
         for atom in self.halogen_bonds]
        for atom in [haldon.x for haldon in self.ligand.halogenbond_don]:
            if atom.GetIdx() not in involved_atoms:
                unpaired_hal.append(atom)
        return unpaired_hba, unpaired_hbd, unpaired_hal

    def refine_hydrophobic(self, all_h, pistacks):
        """Apply several rules to reduce the number of hydrophobic interactions."""
        sel = {}
        # 1. Rings interacting via stacking can't have additional hydrophobic
        # contacts between each other.
        for pistack, h in itertools.product(pistacks, all_h):
            h1, h2 = h.bsatom.GetIdx(), h.ligatom.GetIdx()
            brs, lrs = [p1.GetIdx() for p1 in pistack.proteinring.atoms], [
                p2.GetIdx() for p2 in pistack.ligandring.atoms]
            if h1 in brs and h2 in lrs:
                sel[(h1, h2)] = "EXCLUDE"
        hydroph = [
            h for h in all_h if not (
                h.bsatom.GetIdx(),
                h.ligatom.GetIdx()) in sel]
        sel2 = {}
        #  2. If a ligand atom interacts with several binding site atoms in the same residue,
        #  keep only the one with the closest distance
        for h in hydroph:
            if not (h.ligatom.GetIdx(), h.resnr) in sel2:
                sel2[(h.ligatom.GetIdx(), h.resnr)] = h
            else:
                if sel2[(h.ligatom.GetIdx(), h.resnr)].distance > h.distance:
                    sel2[(h.ligatom.GetIdx(), h.resnr)] = h
        hydroph = [h for h in sel2.values()]
        hydroph_final = []
        bsclust = {}
        # 3. If a protein atom interacts with several neighboring ligand atoms,
        # just keep the one with the closest dist
        for h in hydroph:
            if h.bsatom.GetIdx() not in bsclust:
                bsclust[h.bsatom.GetIdx()] = [h, ]
            else:
                bsclust[h.bsatom.GetIdx()].append(h)

        idx_to_h = {}
        for bs in [a for a in bsclust if len(bsclust[a]) == 1]:
            hydroph_final.append(bsclust[bs][0])

        # A list of tuples with the idx of an atom and one of its neighbours is
        # created
        for bs in [a for a in bsclust if not len(bsclust[a]) == 1]:
            tuples = []
            all_idx = [i.ligatom.GetIdx() for i in bsclust[bs]]
            for b in bsclust[bs]:
                idx = b.ligatom.GetIdx()
                neigh = [na for na in b.ligatom.GetNeighbors()]
                for n in neigh:
                    n_idx = n.GetIdx()
                    if n_idx in all_idx:
                        if n_idx < idx:
                            tuples.append((n_idx, idx))
                        else:
                            tuples.append((idx, n_idx))
                        idx_to_h[idx] = b

            tuples = list(set(tuples))
            tuples = sorted(tuples, key=itemgetter(1))
            # Cluster connected atoms (i.e. find hydrophobic patches)
            clusters = cluster_doubles(tuples)

            for cluster in clusters:
                min_dist = float('inf')
                min_h = None
                for atm_idx in cluster:
                    h = idx_to_h[atm_idx]
                    if h.distance < min_dist:
                        min_dist = h.distance
                        min_h = h
                hydroph_final.append(min_h)
        before, reduced = len(all_h), len(hydroph_final)
        if not before == 0 and not before == reduced:
            # raise RuntimeWarning(f'reduced number of hydrophobic contacts from {before} to {reduced}')
            #print(f'reduced number of hydrophobic contacts from {before} to {reduced}')
            pass
        return hydroph_final

    def refine_hbonds_ldon(self, all_hbonds, salt_lneg, salt_pneg):
        """Refine selection of hydrogen bonds. Do not allow groups which already form salt bridges to form H-Bonds."""
        i_set = {}
        for hbond in all_hbonds:
            i_set[hbond] = False
            for salt in salt_pneg:
                protidx, ligidx = [
                    at.GetIdx() for at in salt.negative.atoms], [
                    at.GetIdx() for at in salt.positive.atoms]
                if hbond.d.GetIdx() in ligidx and hbond.a.GetIdx() in protidx:
                    i_set[hbond] = True
            for salt in salt_lneg:
                protidx, ligidx = [
                    at.GetIdx() for at in salt.positive.atoms], [
                    at.GetIdx() for at in salt.negative.atoms]
                if hbond.d.GetIdx() in ligidx and hbond.a.GetIdx() in protidx:
                    i_set[hbond] = True

        # Allow only one hydrogen bond per donor, select interaction with
        # larger donor angle
        second_set = {}
        hbls = [k for k in i_set.keys() if not i_set[k]]
        for hbl in hbls:
            if hbl.d.GetIdx() not in second_set:
                second_set[hbl.d.GetIdx()] = (hbl.angle, hbl)
            else:
                if second_set[hbl.d.GetIdx()][0] < hbl.angle:
                    second_set[hbl.d.GetIdx()] = (hbl.angle, hbl)
        return [hb[1] for hb in second_set.values()]

    def refine_hbonds_pdon(self, all_hbonds, salt_lneg, salt_pneg):
        """Refine selection of hydrogen bonds. Do not allow groups which already form salt bridges to form H-Bonds with
        atoms of the same group.
        """
        i_set = {}
        for hbond in all_hbonds:
            i_set[hbond] = False
            for salt in salt_lneg:
                protidx, ligidx = [
                    at.GetIdx() for at in salt.positive.atoms], [
                    at.GetIdx() for at in salt.negative.atoms]
                if hbond.a.GetIdx() in ligidx and hbond.d.GetIdx() in protidx:
                    i_set[hbond] = True
            for salt in salt_pneg:
                protidx, ligidx = [
                    at.GetIdx() for at in salt.negative.atoms], [
                    at.GetIdx() for at in salt.positive.atoms]
                if hbond.a.GetIdx() in ligidx and hbond.d.GetIdx() in protidx:
                    i_set[hbond] = True

        # Allow only one hydrogen bond per donor, select interaction with
        # larger donor angle
        second_set = {}
        hbps = [k for k in i_set.keys() if not i_set[k]]
        for hbp in hbps:
            if hbp.d.GetIdx() not in second_set:
                second_set[hbp.d.GetIdx()] = (hbp.angle, hbp)
            else:
                if second_set[hbp.d.GetIdx()][0] < hbp.angle:
                    second_set[hbp.d.GetIdx()] = (hbp.angle, hbp)
        return [hb[1] for hb in second_set.values()]

    def refine_pi_cation_laro(self, all_picat, stacks):
        """Just important for constellations with histidine involved. If the histidine ring is positioned in stacking
        position to an aromatic ring in the ligand, there is in most cases stacking and pi-cation interaction reported
        as histidine also carries a positive charge in the ring. For such cases, only report stacking.
        """
        i_set = []
        for picat in all_picat:
            exclude = False
            for stack in stacks:
                if whichrestype(
                        stack.proteinring.atoms[0]) == 'HIS' and picat.ring.obj == stack.ligandring.obj:
                    exclude = True
            if not exclude:
                i_set.append(picat)
        return i_set

    def refine_water_bridges(self, wbridges, hbonds_ldon, hbonds_pdon):
        """A donor atom already forming a hydrogen bond is not allowed to form a water bridge. Each water molecule
        can only be donor for two water bridges, selecting the constellation with the omega angle closest to 110 deg."""
        donor_atoms_hbonds = [hb.d.GetIdx()
                              for hb in hbonds_ldon + hbonds_pdon]
        wb_dict = {}
        wb_dict2 = {}
        omega = 110.0

        # Just one hydrogen bond per donor atom
        for wbridge in [
                wb for wb in wbridges if wb.d.GetIdx() not in donor_atoms_hbonds]:
            if (wbridge.water.GetIdx(), wbridge.a.GetIdx()) not in wb_dict:
                wb_dict[(wbridge.water.GetIdx(), wbridge.a.GetIdx())] = wbridge
            else:
                if abs(omega - wb_dict[(wbridge.water.GetIdx(),
                                        wbridge.a.GetIdx())].w_angle) < abs(omega - wbridge.w_angle):
                    wb_dict[(wbridge.water.GetIdx(),
                             wbridge.a.GetIdx())] = wbridge
        for wb_tuple in wb_dict:
            water, acceptor = wb_tuple
            if water not in wb_dict2:
                wb_dict2[water] = [
                    (abs(
                        omega -
                        wb_dict[wb_tuple].w_angle),
                        wb_dict[wb_tuple]),
                ]
            elif len(wb_dict2[water]) == 1:
                wb_dict2[water].append(
                    (abs(omega - wb_dict[wb_tuple].w_angle), wb_dict[wb_tuple]))
                wb_dict2[water] = sorted(wb_dict2[water], key=lambda x: x[0])
            else:
                if wb_dict2[water][1][0] < abs(
                        omega - wb_dict[wb_tuple].w_angle):
                    wb_dict2[water] = [ wb_dict2[water][0],
                        (wb_dict[wb_tuple].w_angle, wb_dict[wb_tuple])]

        filtered_wb = []
        for fwbridges in wb_dict2.values():
            [filtered_wb.append(fwb[1]) for fwb in fwbridges]
        return filtered_wb


def get_interactions(mol_protein, mol_ligand, mol_waters=None, pdbid=None):
    data = namedtuple("interaction", "lig prot interactions")
    lig = Ligand(mol_ligand, mol_waters)
    prot = Protein(mol_protein)

    interactions = PLInteraction(lig, prot, pdbid)
    return data(lig=lig, prot=prot, interactions=interactions)
