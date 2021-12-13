import numpy as np
from scipy.spatial.distance import cdist
from schrodinger.structure import StructureReader, StructureWriter, PDBWriter
import argparse

def get_pocket_atomidx(protein_coords, ligand_coords):
    distacne_matrics = cdist(protein_coords, ligand_coords)
    pocket_atoms_id = list(set(np.where(distacne_matrics <= 5.5)[0]))
    pocket_atoms_idx = [idx+1 for idx in pocket_atoms_id]
    return pocket_atoms_idx

def get_residue_atomidx(pocket_atoms_idx, prot):
    atoms_idx = []
    for idx in pocket_atoms_idx:
        atom = prot.atom[idx]
        residue = atom.getResidue()
        atoms_idx.extend(residue.getAtomIndices())
    atoms_idx = list(set(atoms_idx))
    atoms_idx = sorted(atoms_idx,  reverse=False)
    return atoms_idx

def get_pocket(protein_file, ligand_file, pocket_file):
    prot = next(StructureReader(protein_file))
    lig = next(StructureReader(ligand_file))
    protein_coords = prot.getXYZ()
    ligand_coords = lig.getXYZ()
    pocket_atoms_idx = get_pocket_atomidx(protein_coords, ligand_coords)
    atoms_idx = get_residue_atomidx(pocket_atoms_idx, prot)
    pocket = prot.extract(atoms_idx, copy_props=False)
    with PDBWriter(pocket_file) as writer:
        writer.append(pocket)
    return

if __name__=="__main__":
    import sys
    protein_file=sys.argv[1]
    ligand_file=sys.argv[2]
    pocket_file=sys.argv[3]
    get_pocket(protein_file, ligand_file, pocket_file)
