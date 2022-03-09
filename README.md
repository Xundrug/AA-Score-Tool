# AA-Score
The protein-ligand scoring function plays an important role in computer-aided drug discovery, which is heavily used in virtual screening and lead optimization. In this study, we developed a new empirical protein-ligand scoring function, which is a linear combination of empirical energy components, including hydrogen bond, van der Waals, electrostatic, hydrophobic, π-stacking, π-cation, and metal-ligand interaction. Different from previous empirical scoring functions, AA-Score uses several amino acid-specific empirical interaction components. We tested AA-Score on several test sets. The resulting performance shows AA-Score performs well on scoring, docking, and ranking compared with other widely used traditional scoring functions. Our results suggest that AA-Score gains substantial improvements from using detailed protein-ligand interaction components.

## Requirements

* Python 3.6
* Openbabel >= 3.0
* RDKit >= 2019 (http://www.rdkit.org/docs/Install.html)
* numpy 1.18.1
* scipy
* pandas 0.25.3
* py3dmol
* schrodinger (optional)
* biopandas

You also can create the python environment by conda configure file:
```
conda env create -f environment.yaml 
```

## Usage
This is a example pdb for model prediction, including protein and ligand file:
```
data/5otc
```
This tool only support command line. First, you should prepare the protein-ligand complex file, including add hydrogens, remove water and assign bond order, etc. Then, you should use the scripts to select binding pocket from protein by referring to the atom position of the ligand. 
```
get_pocket_by_biopandas.py or get_pocket_by_schrodinger.py 
```
The purpose of `AA_Score.py` is to predict binding affinity by AA-Score:
```
AA_Score.py
usage: AA_Score.py [-h] [--Rec REC] [--Lig LIG] [--Out OUT]

parse AA Score prediction parameters

optional arguments:
  -h, --help  show this help message and exit
  --Rec REC   the file of binding pocket, only support PDB format
  --Lig LIG   the file of ligands, support mol2, mol, sdf, PDB
  --Out OUT   the output file for recording scores
```
If you want to print the result in the screen, please ignore the parameter `--out`:
```
python AA_Score.py --Rec data/2reg/pocket.pdb --Lig data/2reg/Lig.sdf
2REG - minimized_ligand1 -5.790370939906665
```

The purpose of `AA_fp.ipynb` is to calculate interaction fingerprint by AA-Score.

The purpose of the code is to calculate the binding affinity by AA-Score using jupyter:
```
AA_Score.ipynb
```
The purpose of the code is to visualize different interactions by AA-Score, including hydrogen bond, pi-pi stacking, pi-cation, salt bridge:
```
show_interactions.ipynb
```

