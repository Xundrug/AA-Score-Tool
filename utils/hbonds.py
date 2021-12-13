#!/usr/bin/env python
# coding: utf-8

from rdkit import Chem
from rdkit.Chem import AllChem
from collections import namedtuple

import numpy as np

def get_neighbor_bond_type(atom):
    names = {Chem.rdchem.BondType.SINGLE: 's',
             Chem.rdchem.BondType.DOUBLE: 'd',
             Chem.rdchem.BondType.TRIPLE: 't',
             Chem.rdchem.BondType.AROMATIC: 'a',
             Chem.rdchem.BondType.OTHER: 'o'}
    bonds = atom.GetBonds()
    bond_types = []
    for bond in bonds:
        bond_types.append(names[bond.GetBondType()])
    return bond_types

def is_amine(atom):
    n_atomics = [a.GetAtomicNum() for a in atom.GetNeighbors()]
    if atom.GetAtomicNum() == 7 and n_atomics.count(1) >=1:
        return True
    else:
        return False

def is_carbonyl(atom):
    n_atomics = [a.GetAtomicNum() for a in atom.GetNeighbors()]
    bond_types = get_neighbor_bond_type(atom)
    if (atom.GetAtomicNum() == 8 
        and bond_types.count('d') == 1 
        and n_atomics.count(1) == 0
        and n_atomics.count(6) == 1):
        return True
    else:
        return False

def is_hydroxy(atom):
    n_atomics = [a.GetAtomicNum() for a in atom.GetNeighbors()]
    bond_types = get_neighbor_bond_type(atom)
    if (atom.GetAtomicNum() == 8 
        and bond_types.count('s') == 2
        and n_atomics.count(1) == 1):
        return True
    else:
        return False

def is_carboxylate(atom):
    neighbors = atom.GetNeighbors()
    if atom.GetAtomicNum() == 8 and len(neighbors) == 1:
        n_atom = neighbors[0]
        if n_atom.GetAtomicNum() == 6:
            n_neighbors = n_atom.GetNeighbors()
            n_atomics = [a.GetAtomicNum() for a in n_neighbors]
            n_neighbors_num = [len(a.GetNeighbors()) for a in n_neighbors]
            if n_atomics.count(8) == 2 and n_neighbors_num.count(1) == 2:
                return True
    return False

def is_sulfurH0(atom):
    bond_types = get_neighbor_bond_type(atom)
    atomics = [a.GetAtomicNum() for a in atom.GetNeighbors()]
    if atom.GetAtomicNum() == 16 and bond_types.count("s") == 2 and atomics.count(1) == 0:
        return True
    return False

def is_sulfurH1(atom):
    bond_types = get_neighbor_bond_type(atom)
    atomics = [a.GetAtomicNum() for a in atom.GetNeighbors()]
    if atom.GetAtomicNum() == 16 and bond_types.count("s") == 2 and atomics.count(1) == 1:
        return True
    return False

def is_pyridine(atom):
    aromatics = [a.GetIsAromatic() for a in atom.GetNeighbors()]
    if atom.GetAtomicNum() == 7 and atom.GetIsAromatic() and aromatics.count(True) == 2:
        return True
    else:
        return False

def is_ether(atom):
    bond_types = get_neighbor_bond_type(atom)
    n_atomic = [a.GetAtomicNum() for a in atom.GetNeighbors()]
    if atom.GetAtomicNum() == 8 and bond_types.count('s') == 2 and n_atomic.count(6) == 2:
        return True
    else:
        return False

def is_phosphoric_o(atom):
    n_atomic = [a.GetAtomicNum() for a in atom.GetNeighbors()]
    if atom.GetAtomicNum() == 8 and len(n_atomic) == 1 and n_atomic.count(15) == 1:
        return True
    else:
        return False

def is_fluorine(atom):
    if atom.GetAtomicNum() == 9:
        return True
    else:
        return False

def is_water(atom):
    n_atomic = [a.GetAtomicNum() for a in atom.GetNeighbors()]
    if atom.GetAtomicNum() == 8 and n_atomic.count(1) == 2 and len(n_atomic) == 2:
        return True
    else:
        return False

def get_group(atom):
    if is_amine(atom):
        group = "amine"
    elif is_carbonyl(atom):
        group = "carbonyl"
    elif is_carboxylate(atom):
        group = "carboxylate"
    elif is_ether(atom):
        group = "ether"
    elif is_hydroxy(atom):
        group = "hydroxy"
    elif is_pyridine(atom):
        group = "pyridine"
    elif is_sulfurH0(atom):
        group = "sulfurH0"
    elif is_sulfurH1(atom):
        group = "sulfurH1"
    elif is_phosphoric_o(atom):
        group = "phosphoric"
    elif is_fluorine(atom):
        group = "fluorine"
    elif is_water(atom):
        group = "water"
    else:
        raise RuntimeError("group of atom can't be classify, {} {}".format(atom.GetSymbol(), atom.GetIdx()))
    return group

radii_angle_dict ={
    "carbonyl_amine": [2.9, 135], "amine_carbonyl": [2.9, 135], 
    "carbonyl_hydroxy": [2.8, 135], "hydroxy_carbonyl": [2.8, 135],
    "carboxylate_amine": [2.9, 135], "amine_carboxylate": [2.9, 135],
    "carboxylate_hydroxy": [2.8, 135], "hydroxy_carboxylate": [2.8, 135],
    "amine_pyridine": [2.7, 120], "pyridine_amine": [2.7, 120], 
    "hydroxy_pyridine": [2.9, 120], "pyridine_hydroxy": [2.9, 120], 
    "amine_ether": [2.9, 109.5], "ether_amine": [2.9, 109.5], 
    "hydroxy_ether": [2.8, 109.5], "ether_hydroxy": [2.8, 109.5],
    "sulfurH0_amine": [3.2, 109.5], "amine_sulfurH0": [3.2, 109.5],
    "sulfurH0_hydroxy": [3.2, 109.5], "hydroxy_sulfurH0": [3.2, 109.5],
    "amine_hydroxy": [2.9, 135], "hydroxy_amine": [2.9, 135],
    "hydroxy_hydroxy": [2.8, 120], 
    "phosphoric_amine": [2.9, 135], "amine_phosphoric": [2.9, 135],
    "phosphoric_hydroxy": [2.8, 135], "hydroxy_phosphoric": [2.8, 135],
    "sulfurH1_carbonyl": [2.9, 135], "carbonyl_sulfurH1": [2.9, 135],
    "sulfurH1_carboxylate": [2.8, 135], "carboxylate_sulfurH1": [2.8, 135],
    "sulfurH1_pyridine": [2.9, 120], "pyridine_sulfurH1": [2.9, 120],
    "sulfurH1_ether": [2.8, 109.5], "ether_sulfurH1": [2.8, 109.5],
    "amine_fluorine": [2.9, 180], "fluorine_amine": [2.9, 180],
    "hydroxy_fluorine": [2.9, 180], "fluorine_hydroxy": [2.9, 180],
    "water_carbonyl": [2.8, 135], "carbonyl_water": [2.8, 135],
    "carboxylate_water": [2.8, 135], "water_carboxylate": [2.8, 135],
    "water_pyridine": [2.9, 120], "pyridine_water": [2.9, 120],
    "water_ether": [2.8, 109.5], "ether_water": [2.8, 109.5],
    "sulfurH0_water": [3.2, 109.5], "water_sulfurH0": [3.2, 109.5],
    "amine_water": [2.9, 135], "water_amine": [2.9, 135],
    "phosphoric_water": [2.8, 135], "water_phosphoric": [2.8, 135],
    "water_fluorine": [2.9, 180], "fluorine_water": [2.9, 180],
}

def get_radii_angle(atom1, atom2):
    data = namedtuple('radii_angle', 'radii0 angle0')
    group1 = get_group(atom1)
    group2 = get_group(atom2)
    label = group1 + "_" + group2
    
    param = radii_angle_dict[label]
    return data(radii0=param[0], angle0=param[1])

def calc_hbond_strength(hb):
    dist = hb.distance_ad
    return  -(1 / (1 + np.power(dist / 2.6, 6))) / 0.58


