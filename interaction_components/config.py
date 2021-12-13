__version__ = '2.1.3'
__maintainer__ = 'PharmAI GmbH (2020) - www.pharm.ai - hello@pharm.ai'

import logging

DEFAULT_LOG_LEVEL = logging.INFO
VERBOSE = False  # Set verbose mode
QUIET = False  # Set verbose mode
SILENT = False  # Set verbose mode
MAXTHREADS = 1  # Maximum number of main threads for binding site visualization
XML = False
TXT = False
PICS = False
PYMOL = False
STDOUT = False
RAWSTRING = False  # use raw strings for input / output
OUTPATH = './'
BASEPATH = './'
BREAKCOMPOSITE = False  # Break up composite ligands with covalent bonds
ALTLOC = False  # Consider alternate locations
PLUGIN_MODE = False  # Special mode for PLIP in Plugins (e.g. PyMOL)
NOFIX = False  # Turn off fixing of errors in PDB files
NOFIXFILE = False  # Turn off writing to files for fixed PDB structures
PEPTIDES = []  # Definition which chains should be considered as peptide ligands
INTRA = None
KEEPMOD = False
DNARECEPTOR = False
OUTPUTFILENAME = "report"  # Naming for the TXT and XML report files
NOPDBCANMAP = False  # Skip calculation of mapping canonical atom order: PDB atom order
# Do not add hydrogen bonds (in case already present in the structure)
NOHYDRO = False

# Configuration file for Protein-Ligand Interaction Profiler (PLIP)
# Set thresholds for detection of interactions

# Thresholds for detection (global variables)
BS_DIST = 7.5  # Determines maximum distance to include binding site residues
# Determines allowed deviation from planarity in aromatic rings
AROMATIC_PLANARITY = 5.0
MIN_DIST = 0.5  # Minimum distance for all distance thresholds
# Some distance thresholds were extended (max. 1.0A) if too restrictive
# too account for low-quality structures
HYDROPH_DIST_MAX = 4.0  # Distance cutoff for detection of hydrophobic contacts
# Max. distance between hydrogen bond donor and acceptor (Hubbard &
# Haider, 2001) + 0.6 A
HBOND_DIST_MAX = 4.1
#HBOND_DON_ANGLE_MIN = 100
# Min. angle at the hydrogen bond donor (Hubbard & Haider, 2001) + 10
HBOND_DON_ANGLE_MIN = 120
# Max. distance for parallel or offset pistacking (McGaughey, 1998)
PISTACK_DIST_MAX = 5.5
# Max. Deviation from parallel or perpendicular orientation (in degrees)
PISTACK_ANG_DEV = 30
# Maximum offset of the two rings (corresponds to the radius of benzene +
# 0.5 A)
PISTACK_OFFSET_MAX = 2.0
# Max. distance between charged atom and aromatic ring center (Gallivan
# and Dougherty, 1999)
PICATION_DIST_MAX = 6.0
# Max. distance between centers of charge for salt bridges (Barlow and
# Thornton, 1983) + 1.5
SALTBRIDGE_DIST_MAX = 5.5
# Max. distance between oxy. and halogen (Halogen bonds in biological
# molecules., Auffinger)+0.5
HALOGEN_DIST_MAX = 4.0
# Optimal acceptor angle (Halogen bonds in biological molecules., Auffinger)
HALOGEN_ACC_ANGLE = 120
# Optimal donor angle (Halogen bonds in biological molecules., Auffinger)
HALOGEN_DON_ANGLE = 165
HALOGEN_ANGLE_DEV = 30  # Max. deviation from optimal angle
# Min. distance between water oxygen and polar atom (Jiang et al., 2005) -0.1
WATER_BRIDGE_MINDIST = 2.5
# Max. distance between water oxygen and polar atom (Jiang et al., 2005) +0.5
WATER_BRIDGE_MAXDIST = 4.1
# Min. angle between acceptor, water oxygen and donor hydrogen (Jiang et
# al., 2005) - 9
WATER_BRIDGE_OMEGA_MIN = 71
# Max. angle between acceptor, water oxygen and donor hydrogen (Jiang et
# al., 2005)
WATER_BRIDGE_OMEGA_MAX = 140
# WATER_BRIDGE_THETA_MIN = 100  # Min. angle between water oxygen, donor
# hydrogen and donor atom (Jiang et al., 2005)
WATER_BRIDGE_THETA_MIN = 120
# Max. distance between metal ion and interacting atom (Harding, 2001)
METAL_DIST_MAX = 3.0

# Other thresholds
MAX_COMPOSITE_LENGTH = 200  # Filter out ligands with more than 200 fragments

#########
# Names #
#########

# Names of RNA and DNA residues to be considered (detection by name)
RNA = ['U', 'A', 'C', 'G']
DNA = ['DT', 'DA', 'DC', 'DG']

#############
# Whitelist #
#############

# Metal cations which can be complexed

METAL_IONS = [
    'CA',
    'CO',
    'MG',
    'MN',
    'FE',
    'CU',
    'ZN',
    'FE2',
    'FE3',
    'FE4',
    'LI',
    'NA',
    'K',
    'RB',
    'SR',
    'CS',
    'BA',
    'CR',
    'NI',
    'FE1',
    'NI',
    'RU',
    'RU1',
    'RH',
    'RH1',
    'PD',
    'AG',
    'CD',
    'LA',
    'W',
    'W1',
    'OS',
    'IR',
    'PT',
    'PT1',
    'AU',
    'HG',
    'CE',
    'PR',
    'SM',
    'EU',
    'GD',
    'TB',
    'YB',
    'LU',
    'AL',
    'GA',
    'IN',
    'SB',
    'TL',
    'PB']

##############
# Blacklists #
##############

# Other Ions/Atoms (not yet supported)
anions = ['CL', 'IOD', 'BR']
other = ['MO', 'RE', 'HO']
UNSUPPORTED = anions + other

# BioLiP list of suspicious ligands from http://zhanglab.ccmb.med.umich.edu/BioLiP/ligand_list (2014-07-10)
# Add ligands here to get warnings for possible artifacts.
biolip_list = ['ACE', 'HEX', 'TMA', 'SOH', 'P25', 'CCN', 'PR', 'PTN', 'NO3', 'TCN', 'BU1', 'BCN', 'CB3', 'HCS', 'NBN',
               'SO2', 'MO6', 'MOH', 'CAC', 'MLT', 'KR', '6PH', 'MOS', 'UNL', 'MO3', 'SR', 'CD3', 'PB', 'ACM', 'LUT',
               'PMS', 'OF3', 'SCN', 'DHB', 'E4N', '13P', '3PG', 'CYC', 'NC', 'BEN', 'NAO', 'PHQ', 'EPE', 'BME', 'TB',
               'ETE', 'EU', 'OES', 'EAP', 'ETX', 'BEZ', '5AD', 'OC2', 'OLA', 'GD3', 'CIT', 'DVT', 'OC6', 'MW1', 'OC3',
               'SRT', 'LCO', 'BNZ', 'PPV', 'STE', 'PEG', 'RU', 'PGE', 'MPO', 'B3P', 'OGA', 'IPA', 'LU', 'EDO', 'MAC',
               '9PE', 'IPH', 'MBN', 'C1O', '1PE', 'YF3', 'PEF', 'GD', '8PE', 'DKA', 'RB', 'YB', 'GGD', 'SE4', 'LHG',
               'SMO', 'DGD', 'CMO', 'MLI', 'MW2', 'DTT', 'DOD', '7PH', 'PBM', 'AU', 'FOR', 'PSC', 'TG1', 'KAI', '1PG',
               'DGA', 'IR', 'PE4', 'VO4', 'ACN', 'AG', 'MO4', 'OCL', '6UL', 'CHT', 'RHD', 'CPS', 'IR3', 'OC4', 'MTE',
               'HGC', 'CR', 'PC1', 'HC4', 'TEA', 'BOG', 'PEO', 'PE5', '144', 'IUM', 'LMG', 'SQU', 'MMC', 'GOL', 'NVP',
               'AU3', '3PH', 'PT4', 'PGO', 'ICT', 'OCM', 'BCR', 'PG4', 'L4P', 'OPC', 'OXM', 'SQD', 'PQ9', 'BAM', 'PI',
               'PL9', 'P6G', 'IRI', '15P', 'MAE', 'MBO', 'FMT', 'L1P', 'DUD', 'PGV', 'CD1', 'P33', 'DTU', 'XAT', 'CD',
               'THE', 'U1', 'NA', 'MW3', 'BHG', 'Y1', 'OCT', 'BET', 'MPD', 'HTO', 'IBM', 'D01', 'HAI', 'HED', 'CAD',
               'CUZ', 'TLA', 'SO4', 'OC5', 'ETF', 'MRD', 'PT', 'PHB', 'URE', 'MLA', 'TGL', 'PLM', 'NET', 'LAC', 'AUC',
               'UNX', 'GA', 'DMS', 'MO2', 'LA', 'NI', 'TE', 'THJ', 'NHE', 'HAE', 'MO1', 'DAO', '3PE', 'LMU', 'DHJ',
               'FLC', 'SAL', 'GAI', 'ORO', 'HEZ', 'TAM', 'TRA', 'NEX', 'CXS', 'LCP', 'HOH', 'OCN', 'PER', 'ACY', 'MH2',
               'ARS', '12P', 'L3P', 'PUT', 'IN', 'CS', 'NAW', 'SB', 'GUN', 'SX', 'CON', 'C2O', 'EMC', 'BO4', 'BNG',
               'MN5', '__O', 'K', 'CYN', 'H2S', 'MH3', 'YT3', 'P22', 'KO4', '1AG', 'CE', 'IPL', 'PG6', 'MO5', 'F09',
               'HO', 'AL', 'TRS', 'EOH', 'GCP', 'MSE', 'AKR', 'NCO', 'PO4', 'L2P', 'LDA', 'SIN', 'DMI', 'SM', 'DTD',
               'SGM', 'DIO', 'PPI', 'DDQ', 'DPO', 'HCA', 'CO5', 'PD', 'OS', 'OH', 'NA6', 'NAG', 'W', 'ENC', 'NA5',
               'LI1', 'P4C', 'GLV', 'DMF', 'ACT', 'BTB', '6PL', 'BGL', 'OF1', 'N8E', 'LMT', 'THM', 'EU3', 'PGR', 'NA2',
               'FOL', '543', '_CP', 'PEK', 'NSP', 'PEE', 'OCO', 'CHD', 'CO2', 'TBU', 'UMQ', 'MES', 'NH4', 'CD5', 'HTG',
               'DEP', 'OC1', 'KDO', '2PE', 'PE3', 'IOD', 'NDG', 'CL', 'HG', 'F', 'XE', 'TL', 'BA', 'LI', 'BR', 'TAU',
               'TCA', 'SPD', 'SPM', 'SAR', 'SUC', 'PAM', 'SPH', 'BE7', 'P4G', 'OLC', 'OLB', 'LFA', 'D10', 'D12', 'DD9',
               'HP6', 'R16', 'PX4', 'TRD', 'UND', 'FTT', 'MYR', 'RG1', 'IMD', 'DMN', 'KEN', 'C14', 'UPL', 'CMJ', 'ULI',
               'MYS', 'TWT', 'M2M', 'P15', 'PG0', 'PEU', 'AE3', 'TOE', 'ME2', 'PE8', '6JZ', '7PE', 'P3G', '7PG', 'PG5',
               '16P', 'XPE', 'PGF', 'AE4', '7E8', '7E9', 'MVC', 'TAR', 'DMR', 'LMR', 'NER', '02U', 'NGZ', 'LXB', 'A2G',
               'BM3', 'NAA', 'NGA', 'LXZ', 'PX6', 'PA8', 'LPP', 'PX2', 'MYY', 'PX8', 'PD7', 'XP4', 'XPA', 'PEV', '6PE',
               'PEX', 'PEH', 'PTY', 'YB2', 'PGT', 'CN3', 'AGA', 'DGG', 'CD4', 'CN6', 'CDL', 'PG8', 'MGE', 'DTV', 'L44',
               'L2C', '4AG', 'B3H', '1EM', 'DDR', 'I42', 'CNS', 'PC7', 'HGP', 'PC8', 'HGX', 'LIO', 'PLD', 'PC2', 'PCF',
               'MC3', 'P1O', 'PLC', 'PC6', 'HSH', 'BXC', 'HSG', 'DPG', '2DP', 'POV', 'PCW', 'GVT', 'CE9', 'CXE', 'C10',
               'CE1', 'SPJ', 'SPZ', 'SPK', 'SPW', 'HT3', 'HTH', '2OP', '3NI', 'BO3', 'DET', 'D1D', 'SWE', 'SOG']

atom_prop_dict = {'TRP_CD2': (), 'TRP_CD1': (), 'LYS_N': ('Donor',), 'LYS_O': ('Acceptor',), 'LYS_C': ('NegIonizable',),
                  'LEU_C': (), 'ASP_OXT': ('NegIonizable', 'Acceptor'), 'LEU_N': ('Donor',), 'LEU_O': ('Acceptor',), 'PHE_CG': (),
                  'SER_CB': (), 'SER_CA': (), 'HIS_CD2': ('PosIonizable',), 'GLU_OE1': ('NegIonizable', 'Acceptor'),
                  'LYS_NZ': ('Donor', 'PosIonizable'), 'ASN_OD1': ('Acceptor',), 'TYR_CD1': (), 'TRP_C': (), 'TRP_N': ('Donor',),
                  'TRP_O': ('Acceptor',), 'LEU_OXT': ('NegIonizable', 'Acceptor'), 'MET_CG': (), 'MET_CE': (), 'MET_CB': (), 'MET_CA': (),
                  'THR_OG1': ('Donor', 'Acceptor'), 'ALA_OXT': ('NegIonizable', 'Acceptor'), 'VAL_CG2': (), 'VAL_CB': (), 'VAL_CA': (), 'VAL_C': (),
                  'VAL_O': ('Acceptor',), 'VAL_N': ('Donor',), 'ASN_CB': (), 'ARG_CZ': ('PosIonizable',), 'ASN_CA': (), 'ASN_CG': (),
                  'HIS_ND1': ('Donor', 'PosIonizable', 'Acceptor'), 'ARG_CA': (), 'MET_SD': (), 'ARG_CD': (), 'ARG_CG': (), 'TYR_CE2': (),
                  'TYR_CE1': (), 'PHE_CZ': (), 'TYR_O': ('Acceptor',), 'TYR_N': ('Donor',), 'TYR_C': (), 'MET_OXT': ('NegIonizable', 'Acceptor'),
                  'SER_C': (), 'SER_N': ('Donor',), 'SER_O': ('Acceptor',), 'TRP_CG': (), 'TRP_CB': (), 'TRP_CA': (), 'PHE_CE2': (),
                  'PHE_CE1': (), 'ARG_CB': (), 'ASP_N': ('Donor',), 'ASP_O': ('Acceptor',), 'ASP_C': (), 'ASP_OD1': ('NegIonizable', 'Acceptor'),
                  'ASP_OD2': ('NegIonizable', 'Donor', 'Acceptor'), 'SER_OG': ('Donor', 'Acceptor'), 'TRP_NE1': ('Donor',), 'PHE_CB': (),
                  'PHE_CA': (), 'VAL_CG1': (), 'GLY_CA': (), 'TRP_CH2': (), 'TYR_OH': ('Donor', 'Acceptor'), 'ARG_OXT': ('NegIonizable', 'Acceptor'),
                  'THR_CB': (), 'THR_CA': (), 'ARG_NE': ('Donor', 'PosIonizable'), 'ILE_OXT': ('NegIonizable', 'Acceptor'),
                  'GLY_OXT': ('NegIonizable', 'Acceptor'), 'TRP_CZ3': (), 'TRP_CZ2': (), 'GLN_NE2': ('Donor',), 'CYS_OXT': ('NegIonizable', 'Acceptor'), 'LEU_CB': (), 'GLU_OXT': ('NegIonizable', 'Acceptor'), 'LEU_CG': (), 'PRO_O': ('Acceptor',), 'PRO_N': ('Donor',), 'PRO_C': (),
                  'ILE_CD1': (), 'ASN_OXT': ('NegIonizable', 'Acceptor'), 'ARG_N': ('Donor',), 'ARG_O': ('Acceptor',), 'ARG_C': (), 'ILE_CA': (),
                  'ILE_CB': (), 'ARG_NH1': ('Donor', 'PosIonizable'), 'ARG_NH2': ('Donor', 'PosIonizable'), 'ALA_O': ('Acceptor',), 'GLN_CA': (),
                  'TRP_CE2': (), 'TRP_CE3': (), 'ALA_C': (), 'GLN_CG': (), 'ALA_N': ('Donor',), 'GLN_CD': (), 'GLN_CB': (), 'GLN_OE1': ('Acceptor',),
                  'TYR_OXT': ('NegIonizable', 'Acceptor'), 'LEU_CD2': (), 'PHE_OXT': ('NegIonizable', 'Acceptor'), 'GLN_C': (), 'GLN_O': ('Acceptor',),
                  'GLN_N': ('Donor',), 'VAL_OXT': ('NegIonizable', 'Acceptor'), 'PRO_OXT': ('NegIonizable', 'Acceptor'), 'GLU_N': ('Donor',),
                  'GLU_O': ('Acceptor',), 'LYS_CA': (), 'LYS_CB': (), 'HIS_CE1': ('PosIonizable',), 'LYS_CD': (), 'LYS_CE': (), 'LYS_CG': (),
                  'HIS_OXT': ('NegIonizable', 'Acceptor'), 'GLU_OE2': ('NegIonizable', 'Donor', 'Acceptor'), 'PHE_O': ('Acceptor',),
                  'PHE_N': ('Donor',), 'PHE_C': (), 'CYS_C': (), 'GLU_CB': (), 'GLU_CA': (), 'GLU_CG': (), 'GLU_CD': ('NegIonizable',),
                  'THR_OXT': ('NegIonizable', 'Acceptor'), 'ASN_C': (), 'HIS_CG': ('PosIonizable',), 'HIS_CA': (), 'HIS_CB': (),
                  'LYS_OXT': ('NegIonizable', 'Acceptor'), 'LEU_CA': (), 'ASN_N': ('Donor',), 'ASN_O': ('Acceptor',), 'CYS_O': ('Acceptor',),
                  'CYS_N': ('Donor',), 'TYR_CD2': (), 'ALA_CA': (), 'ALA_CB': (), 'GLU_C': (), 'GLN_OXT': ('NegIonizable', 'Acceptor'),
                  'ASN_ND2': ('Donor',), 'THR_CG2': (), 'PHE_CD1': (), 'PHE_CD2': (), 'ASP_CG': ('NegIonizable',), 'ASP_CA': (), 'ASP_CB': (),
                  'MET_N': ('Donor',), 'MET_O': ('Acceptor',), 'CYS_SG': ('Donor',), 'SER_OXT': ('NegIonizable', 'Acceptor'),
                  'THR_N': ('Donor',), 'THR_O': ('Acceptor',), 'ILE_N': ('Donor',), 'ILE_O': ('Acceptor',), 'THR_C': (), 'ILE_C': (),
                  'MET_C': (), 'ILE_CG1': (), 'ILE_CG2': (), 'GLY_N': ('Donor',), 'GLY_O': ('Acceptor',), 'HIS_N': ('Donor',), 'HIS_O': ('Acceptor',),
                  'HIS_C': (), 'GLY_C': (), 'LEU_CD1': (), 'PRO_CB': (), 'PRO_CA': (), 'PRO_CG': (), 'PRO_CD': (), 'TRP_OXT': ('NegIonizable', 'Acceptor'),
                  'HIS_NE2': ('Donor', 'PosIonizable'), 'CYS_CA': (), 'CYS_CB': (), 'TYR_CZ': (), 'TYR_CG': (), 'TYR_CB': (), 'TYR_CA': ()}
