## Scripts
These scripts are used to extract the binding pocket and assign the bond order for ligands. If you have a liscense of schrodinger, you can use the script `get_pocket_by_schrodinger.py`; otherwise, you should use the script `get_pocket_by_biopandas.py`.

## Usage
The purpose of `get_pocket_by_biopandas.py` is to get the binding pocket:
```
python get_pocket_by_biopandas.py
```
The purpose of `get_pocket_by_schrodinger.py` is to get the binding pocket by schrodinger python api (requires a license of schrodinger):
```
run get_pocket_by_schrodinger.py ../data/2reg/Rec.pdb ../data/2reg/Lig.sdf pocket.pdb
```
The purpose of `repair_ligand_by_schrodinger.py` is to assign bond order for ligand by schrodinger python api (requires a license of schrodinger):
```
usage: repair_ligand_by_schrodinger.py [-h] in_file out_file

positional arguments:
  in_file     Input structure file (Maestro, PDB, or SD, mol2, mol)
  out_file    Input structure file (PDB)

optional arguments:
  -h, --help  show this help message and exit
```

