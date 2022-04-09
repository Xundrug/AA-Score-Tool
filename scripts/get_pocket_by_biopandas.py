from random import randint
import pandas as pd
import os,glob
import numpy as np
import copy
from biopandas.pdb import PandasPdb
from rdkit import Chem
pd.options.mode.chained_assignment = None
import warnings
warnings.simplefilter(action='ignore',category=FutureWarning)
warnings.simplefilter(action='ignore',category=UserWarning)

def deal_pro_and_lig(pro_file, lig_file, pdb_id):
    ppdb = PandasPdb()
    ppdb.read_pdb(pro_file)
    pro = ppdb.df['ATOM']

    pro_cut, pros_near_lig = select_cut_residue(pro, lig, cut=5.5)
    
    ppdb.df['ATOM'] = pro_cut
    
    newmolname = str(randint(1,1000000)).zfill(10)
    name = 'temp/pocket_{}_{}.pdb'.format(pdb_id, newmolname)
    
    ppdb.to_pdb(path=name, records=['ATOM'])
    pmol = Chem.MolFromPDBFile(name, removeHs=False) 
    return pmol, name

def conc(a,b,c):
    return [a,b,c]

def add_xyz(ligu):
    ligu['xyz'] = ligu.apply(lambda row: conc(row['x_coord'],row['y_coord'],row['z_coord']),axis=1)
    return ligu

def cal_dist(a,b):
    a1 = np.array(a)
    b1 = np.array(b)
    dist = np.linalg.norm(a1-b1)
    dist =round(dist,2)
    return dist

def get_min_dist(am, ligu):
    ligu['pro_xyz']=[am]*ligu.shape[0]
    ligu['dist']= ligu.apply(lambda row: cal_dist(row['xyz'], row['pro_xyz']), axis =1)
    md = min(ligu['dist'])
    return md

class GetPocket:
    def __init__(self, ligand_file, protein_file, pdb_id):
        """
        ligand_file: format mol
        protein_file: format pdb
        """
        self.ligand_file = ligand_file
        self.protein_file = protein_file
        self.pdb_id =  pdb_id
        self.ligand_mol = Chem.MolFromMolFile(ligand_file, removeHs=False)
        self.pocket_mol,  self.temp_file = self.process_pro_and_lig()
        self.pocket_path = os.path.basename(protein_file).split(".")[0] + "_pocket.pdb"
        Chem.MolToPDBFile(self.pocket_mol, self.pocket_path)
        
    def process_pro_and_lig(self):
        ppdb = PandasPdb()
        ppdb.read_pdb(self.protein_file)
        protein_biop = ppdb.df['ATOM']
        
        pro_cut, pros_near_lig = self.select_cut_residue(protein_biop, self.ligand_mol, cut=5.5)

        ppdb.df['ATOM'] = pro_cut
        newmolname = str(randint(1,1000000)).zfill(10)
        name = 'pocket_{}_{}.pdb'.format(pdb_id, newmolname)

        ppdb.to_pdb(path=name, records=['ATOM'])
        pmol = Chem.MolFromPDBFile(name, removeHs=False) 
        return pmol, name
    
    def select_cut_residue(self, protein_biop, ligand_mol, cut=5.0):
        """
        pro: biopandas DataFrame
        lig: rdkit mol
        """
        pro = self.cal_pro_min_dist(protein_biop, ligand_mol)

        pro['chain_rid'] = pro.apply(lambda row: 
                                     str(row['chain_id'])+str(row['residue_number']), axis=1)
        pros = pro[pro['min_dist'] < cut]
        pros_near_lig = copy.deepcopy(pros)
        use_res = list(set(list(pros['chain_rid'])))
        pro= pro[pro['chain_rid'].isin(use_res)]
        pro = pro.drop(['chain_rid'],axis=1)
        return pro, pros_near_lig
    
    def get_ligu(self, ligand_mol):
        mol_ligand_conf = ligand_mol.GetConformers()[0]
        pos = mol_ligand_conf.GetPositions()
        df = pd.DataFrame(pos)
        df.columns = ["x_coord", "y_coord","z_coord"]
        return df
    
    def cal_pro_min_dist(self, protein_biop, ligand_mol):
        protein_biop =add_xyz(protein_biop)
        
        ligu = self.get_ligu(ligand_mol)
        ligu = add_xyz(ligu)
        protein_biop['min_dist']= protein_biop.apply(lambda row: get_min_dist(row['xyz'], ligu), axis=1)
        return protein_biop


if __name__=="__main__":
    import sys
    protein_file=sys.argv[1]
    ligand_file=sys.argv[2]
    pdb_id=sys.argv[3]
    get_pocket = GetPocket(ligand_file, protein_file, pdb_id)




