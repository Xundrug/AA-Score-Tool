from interaction_components.plinteraction import get_interactions
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem

from utils.hbonds import calc_hbond_strength
from utils.hydrophobic import calc_hydrophobic
from utils.vdw import calc_vdw
from utils.electrostatic import Electrostatic
import os
import sys
import argparse

residue_names = [
    "HIS",
    "ASP",
    "ARG",
    "PHE",
    "ALA",
    "CYS",
    "GLY",
    "GLN",
    "GLU",
    "LYS",
    "LEU",
    "MET",
    "ASN",
    "SER",
    "TYR",
    "THR",
    "ILE",
    "TRP",
    "PRO",
    "VAL"]


def is_sidechain(atom):
    res = atom.GetPDBResidueInfo()
    atom_name = res.GetName().strip(" ")
    if atom_name in ("C", "CA", "N", "O", "H"):
        return False
    else:
        return True


def create_dict():
    interaction_dict = {}
    for name in residue_names:
        interaction_dict.update({name + "_side": 0})
        interaction_dict.update({name + "_main": 0})
    return interaction_dict


def calc_hb(hbonds, hb_dict):
    for hb in hbonds:
        restype = hb.restype
        sidechain = hb.sidechain
        energy = calc_hbond_strength(hb)
        if restype == "HIN":
            restype = "HIS"
        if restype == "ACE":
            continue
        if sidechain:
            key = restype + "_side"
        else:
            key = restype + "_main"
        hb_dict[key] += energy
    return


def calc_hbonds_descriptor(interactions):
    hb_dict = create_dict()
    calc_hb(interactions.all_hbonds_ldon, hb_dict)
    calc_hb(interactions.all_hbonds_pdon, hb_dict)
    return hb_dict


def calc_hydrophybic_descriptor(interactions):
    hc_dict = create_dict()
    for hc in interactions.hydrophobic_contacts:
        restype = hc.restype
        sidechain = hc.sidechain
        energy = calc_hydrophobic(hc)
        if restype[:2] == "HI" and restype not in residue_names:
            restype = "HIS"
        if restype == "ACE":
            continue
        if sidechain:
            key = restype + "_side"
        else:
            key = restype + "_main"
        hc_dict[key] += energy
    return hc_dict


def calc_vdw_descriptor(result, mol_lig):
    prot = result.prot
    residues = prot.residues
    vdw_dict = create_dict()
    for res in residues:
        main_vdw, side_vdw = calc_vdw(res, mol_lig)
        restype = res.residue_name
        if restype[:2] == "HI" and restype not in residue_names:
            restype = "HIS"
        if restype not in residue_names:
            continue
        vdw_dict[restype + "_side"] += main_vdw
        vdw_dict[restype + "_main"] += side_vdw
    return vdw_dict


def calc_ele_descriptor(result, mol_lig, mol_prot):
    prot = result.prot
    residues = prot.residues
    ele_same_dict = create_dict()
    ele_opposite_dict = create_dict()
    for res in residues:
        ele = Electrostatic(res, mol_lig, mol_prot)
        restype = res.residue_name
        if restype[:2] == "HI" and restype not in residue_names:
            restype = "HIS"
        if restype not in residue_names:
            continue
        ele_same_dict[restype + "_side"] += ele.side_ele_same
        ele_same_dict[restype + "_main"] += ele.main_ele_same

        ele_opposite_dict[restype + "_side"] += ele.side_ele_opposite
        ele_opposite_dict[restype + "_main"] += ele.main_ele_opposite
    return ele_same_dict, ele_opposite_dict


def calc_desolvation_descriptor(result, mol_prot, mol_lig):
    dehyd = Dehydration(mol_prot, mol_lig)
    prot = result.prot

    dehyd_energy = 0
    origin = "protein"
    for at in mol_prot.GetAtoms():
        dehyd_energy += dehyd.calc_atom_dehyd(at, origin)
    origin = "ligand"
    for at in mol_lig.GetAtoms():
        dehyd_energy += dehyd.calc_atom_dehyd(at, origin)
    return dehyd_energy


def calc_metal_complexes(metal):
    dist = metal.distance
    if dist < 2.0:
        return -1.0
    elif 2.0 <= dist < 3.0:
        return -3.0 + dist
    else:
        return 0.0


def calc_metal_descriptor(interactions):
    ml_energy = 0
    for ml in interactions.metal_complexes:
        if ml.target.location != "ligand":
            continue
        energy = calc_metal_complexes(ml)
        ml_energy += energy
    return ml_energy


def calc_pistacking_descriptor(interactions):
    T_pistacking_energy, P_pistacking_energy = 0, 0
    for pis in interactions.pistacking:
        if pis.type == "T":
            T_pistacking_energy += -1
        else:
            P_pistacking_energy += -1
    return T_pistacking_energy, P_pistacking_energy


def calc_pication_laro(interactions):
    pic_dict = create_dict()
    for pic in interactions.pication_laro:
        restype = pic.restype
        sidechain = is_sidechain(pic.charge.atoms[0])
        energy = -1
        if restype[:2] == "HI" and restype not in residue_names:
            restype = "HIS"
        if restype == "ACE":
            continue
        if sidechain:
            key = restype + "_side"
        else:
            key = restype + "_main"
        pic_dict[key] += energy
    return pic_dict


def calc_pication_descriptor(interactions):
    paro_pication_energy, laro_pication_energy = 0, 0
    for pic in interactions.pication_paro:
        paro_pication_energy += -1
    pic_dict = calc_pication_laro(interactions)
    return paro_pication_energy, pic_dict


class Model:
    def __init__(self, arr):
        self.arr = arr

    def predict(self, data):
        data = np.array(data)
        return np.sum(self.arr * data)-0.999


def load_model():
    param = np.load("models/model-final.npy")
    clf = Model(param)
    return clf


def merge_descriptors(
        hb_dict,
        hc_dict,
        vdw_dict,
        ele_same_dict,
        ele_opposite_dict,
        pic_dict,
        metal_ligand,
        tpp_energy,
        ppp_energy,
        ppc_energy,
        rotat):
    line = []
    descriptors = [hb_dict, vdw_dict, ele_same_dict, ele_opposite_dict]
    for des in descriptors:
        for v in des.values():
            line.append(v)
    line.append(sum(hc_dict.values()))
    line.append(sum(pic_dict.values()))
    line.append(metal_ligand)
    line.append(tpp_energy)
    line.append(ppp_energy)
    line.append(ppc_energy)
    line.append(rotat)
    descriptors = np.array(line)
    return descriptors


def calc_score(mol_lig, mol_prot, clf):
    result = get_interactions(mol_prot, mol_lig)
    interactions = result.interactions

    hb_dict = calc_hbonds_descriptor(interactions)
    hc_dict = calc_hydrophybic_descriptor(interactions)
    vdw_dict = calc_vdw_descriptor(result, mol_lig)
    ele_same_dict, ele_opposite_dict = calc_ele_descriptor(
        result, mol_lig, mol_prot)
    metal_ligand = calc_metal_descriptor(interactions)
    tpp_energy, ppp_energy = calc_pistacking_descriptor(interactions)
    ppc_energy, pic_dict = calc_pication_descriptor(interactions)
    rotat = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol_lig)

    descriptors = merge_descriptors(
        hb_dict,
        hc_dict,
        vdw_dict,
        ele_same_dict,
        ele_opposite_dict,
        pic_dict,
        metal_ligand,
        tpp_energy,
        ppp_energy,
        ppc_energy,
        rotat)
    score = clf.predict(descriptors)
    return score

def get_format(ligand_file):
    file_format = os.path.basename( ligand_file ).split(".")[1]
    return file_format

def calc_batch(mol_prot, mol_ligs, output_file, clf):
    for mol_lig in mol_ligs:
        name = mol_lig.GetProp("_Name")
        score = calc_score(mol_lig, mol_prot, clf)

        if output_file:
            with open(output_file, "a") as f:
                f.write(name + "\t" + str(score) + "\n")
        else:
            print( name, score )
    return

def calc_single(mol_prot, mol_lig, output_file, clf):
    name = mol_lig.GetProp("_Name")
    score = calc_score(mol_lig, mol_prot, clf)

    if output_file:
        with open(output_file, "a") as f:
            f.write(name + "\t" + str(score) + "\n")
    else:
        print( name, score )
    return

def calc_single(mol_prot, mol_lig, output_file, clf):
    if not mol_lig:
        raise RuntimeError("RDKit parse the file error")
    name = mol_lig.GetProp("_Name")
    score = calc_score(mol_lig, mol_prot, clf)

    if output_file:
        with open(output_file, "a") as f:
            f.write(name + "\t" + str(score) + "\n")
    else:
        print( name, score )
    return

def predict_dG(mol_prot, mol_lig, output_file=None):
    clf = load_model()
    
    name = mol_lig.GetProp("_Name")
    score = calc_score(mol_lig, mol_prot, clf)

    if output_file:
        with open(output_file, "a") as f:
            f.write(name + "\t" + str(score) + "\n")
    else:
        return name, score

def func():
    parser = argparse.ArgumentParser(description='parse AA Score prediction parameters')
    parser.add_argument('--Rec', type=str, help='the file of binding pocket, only support PDB format')
    parser.add_argument('--Lig', type=str, help='the file of ligands, support mol2, mol, sdf, PDB')
    parser.add_argument('--Out', type=str, help='the output file for recording scores', default=None)
    args = parser.parse_args()
    protein_file = args.Rec
    ligand_file = args.Lig
    output_file = args.Out
    
    clf = load_model()
    mol_prot = Chem.MolFromPDBFile(protein_file, removeHs=False)
    lig_format = get_format(ligand_file)
    if lig_format not in ["sdf", "mol2", "mol", "pdb"]:
        raise RuntimeError("ligand format {} is not supported".format(lig_format))

    if lig_format == "sdf":
        mol_ligs = Chem.SDMolSupplier(ligand_file, removeHs=False)
        calc_batch(mol_prot, mol_ligs, output_file, clf)
    elif lig_format == "mol2":
        mol_lig = Chem.MolFromMol2File(ligand_file, removeHs=False)
        calc_single(mol_prot, mol_lig, output_file, clf)
    elif lig_format == "mol":
        mol_lig = Chem.MolFromMolFile(ligand_file, removeHs=False)
        calc_single(mol_prot, mol_lig, output_file, clf)
    elif lig_format == "pdb":
        mol_lig = Chem.MolFromPDBFile(ligand_file, removeHs=False)
        calc_single(mol_prot, mol_lig, output_file, clf)
    return

if __name__ == "__main__":
    func()
