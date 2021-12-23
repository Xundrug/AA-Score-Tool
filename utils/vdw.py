import numpy as np
from scipy.spatial.distance import  cdist, euclidean
from rdkit.Chem import AllChem

atom_radius = {"N": 1.8, "O": 1.7, "S": 2.0, "P": 2.1, "F": 1.5, "Cl": 1.8,
               "Br": 2.0, "I": 2.2, "C": 1.9, "H": 0.0, "Zn": 0.5, "B": 1.8}


def get_mol_coords(mol):
    conf = mol.GetConformers()[0]
    coords = conf.GetPositions()
    return coords

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

def accelerate_vdw(d0_matrix, dist_matrix):
    vdw = np.sum(np.power((d0_matrix / dist_matrix), 8) - 2 * np.power((d0_matrix / dist_matrix), 4))
    return vdw

def get_d0_matrix(atoms, mol_lig):
    atom_radius_keys =  atom_radius.keys()

    residue_symbols = [a.GetSymbol() for a in atoms]
    ligand_symbols = [a.GetSymbol() for a in mol_lig.GetAtoms()]

    d0_matrix = np.zeros((len(residue_symbols), len(ligand_symbols)), dtype=np.float)
    for idxp, elemp in enumerate(residue_symbols):
        for idxl, eleml in enumerate(ligand_symbols):
            d0 = atom_radius.get(elemp, 0.0) + atom_radius.get(eleml, 0.0)
            d0_matrix[idxp, idxl] = d0
    return d0_matrix

def is_sidechain(atom):
    res = atom.GetPDBResidueInfo()
    atom_name = res.GetName().strip(" ")
    if atom_name in ("C", "CA", "N", "O", "H"):
        return False
    else:
        return True

def calc_vdw_chain(atoms, lig_coords, mol_lig):
    atoms_coords = get_atoms_coords(atoms)
    dist_matrix = cdist(atoms_coords, lig_coords, 'euclidean')
    d0_matrix = get_d0_matrix(atoms, mol_lig)
    vdw = accelerate_vdw(d0_matrix, dist_matrix)
    return vdw

def calc_vdw(residue, mol_lig):
    residue_atoms = [a for a in residue.residue_atoms if a.GetAtomicNum() != 1]
    side_atoms, main_atoms = [], []
    for atom in residue_atoms:
        if is_sidechain(atom):
            side_atoms.append(atom)
        else:
            main_atoms.append(atom)

    #mol_lig = AllChem.RemoveHs(mol_lig)
    lig_coords = get_mol_coords(mol_lig)
    side_vdw = calc_vdw_chain(side_atoms, lig_coords, mol_lig)
    main_vdw = calc_vdw_chain(main_atoms, lig_coords, mol_lig)
    return main_vdw, side_vdw

