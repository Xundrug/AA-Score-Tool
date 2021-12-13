try:
    from openbabel import pybel
except:
    import pybel
from scipy.spatial.distance import  cdist, euclidean
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def get_atom_coords(atom):
    mol = atom.GetOwningMol()
    conf = mol.GetConformers()[0]
    pos = conf.GetAtomPosition(atom.GetIdx())
    return (pos.x, pos.y, pos.z)

def get_atoms_coords(atoms):
    coords = np.zeros((len(atoms), 3))
    for idx, atom in enumerate(atoms):
        coord = get_atom_coords(atom)
        coords[idx, :] = coord
    return coords

def is_sidechain(atom):
    res = atom.GetPDBResidueInfo()
    atom_name = res.GetName().strip(" ")
    if atom_name in ("C", "CA", "N", "O", "H"):
        return False
    else:
        return True

#Eletronic
class Electrostatic(object):
    def __init__(self, residue, mol_lig, mol_prot):
        self.mol_lig = mol_lig
        self.pmol_lig = pybel.readstring("pdb", Chem.MolToPDBBlock(mol_lig))
        self.pmol_prot = pybel.readstring("pdb", Chem.MolToPDBBlock(mol_prot))
        self.residue = residue
        self.residue_atoms = [a for a in residue.residue_atoms if a.GetAtomicNum() != 1]
        self.side_atoms, self.main_atoms = self._classify_atoms(self.residue_atoms)
        self.lig_atoms = [a for a in mol_lig.GetAtoms() if a.GetAtomicNum() != 1]
        self.side_charges = self.get_partial_charge(self.side_atoms, self.pmol_prot)
        self.main_charges = self.get_partial_charge(self.main_atoms, self.pmol_prot)
        self.lig_charges = self.get_partial_charge(self.lig_atoms, self.pmol_lig)
        
        self.side_ele_same, self.side_ele_opposite = self.calc_eletronic(self.side_atoms, self.side_charges)
        self.main_ele_same, self.main_ele_opposite = self.calc_eletronic(self.main_atoms, self.main_charges)
        
    def _classify_atoms(self, atoms):
        side_atoms, main_atoms = [], []
        for a in atoms:
            if is_sidechain(a):
                side_atoms.append(a)
            else:
                main_atoms.append(a)
        return side_atoms, main_atoms

    def get_partial_charge(self, atoms, pmol):
        charges = []
        aidxs = [at.GetIdx() for at in atoms]
        patoms = [pmol.atoms[aidx] for aidx in aidxs]
        for atom in patoms:
            charge = atom.partialcharge
            charges.append(charge)
        return charges

    def get_ele_matrix(self, atoms, charges):
        ele_matrix = np.zeros((len(atoms), len(self.lig_atoms)), dtype=np.float)
        for idxp, p_charge in enumerate(charges):
           for idxl, l_charge in enumerate(self.lig_charges):
               ele = p_charge * l_charge 
               ele_matrix[idxp, idxl] = ele
        return ele_matrix

    def get_dist_matrix(self, atoms):
        residue_coords = get_atoms_coords(atoms)
        lig_coords = get_atoms_coords(self.lig_atoms)
        dist_matrix = cdist(residue_coords, lig_coords, 'euclidean')
        return dist_matrix

    def calc_eletronic(self, atoms, charges):
        ele_matrix = self.get_ele_matrix(atoms, charges)
        dist_matrix = self.get_dist_matrix(atoms)
        ele = ele_matrix / dist_matrix
        
        ele_same, ele_opposite = 0, 0
        for i in range(ele.shape[0]):
            for j in range(ele.shape[1]):
                e = ele[i][j]
                if e <=0:
                    ele_opposite += e
                else:
                    ele_same += e
        return ele_same, ele_opposite



