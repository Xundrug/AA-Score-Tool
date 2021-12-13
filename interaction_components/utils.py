from rdkit import Chem
from rdkit.Chem import AllChem
import itertools
from interaction_components import config
import numpy as np
from collections import namedtuple
from scipy.spatial.distance import euclidean

atom_prop_dict = config.atom_prop_dict
biolip_list = config.biolip_list


def normalize_vector(v):
    """Take a vector and return the normalized vector
    :param v: a vector v
    :returns : normalized vector v
    """
    norm = np.linalg.norm(v)
    return v / norm if not norm == 0 else v


def whichrestype(atom):
    res_info = atom.GetPDBResidueInfo()
    return res_info.GetResidueName().strip(" ") if res_info is not None else None


def whichatomname(atom):
    res_info = atom.GetPDBResidueInfo()
    return res_info.GetName().strip(" ") if res_info is not None else None


def whichresnumber(atom):
    res_info = atom.GetPDBResidueInfo()
    return res_info.GetResidueNumber() if res_info is not None else None


def whichchain(atom):
    res_info = atom.GetPDBResidueInfo()
    return res_info.GetChainId() if res_info is not None else None


def euclidean3d(v1, v2):
    """Faster implementation of euclidean distance for the 3D case."""
    if not (len(v1) == 3 and len(v2) == 3):
        return None
    return euclidean(v1, v2)
    # return np.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] -
    # v2[2]) ** 2)


def is_aromatic(r_atoms):
    for atom in r_atoms:
        if atom.GetIsAromatic():
            return True
        else:
            return False


def vector(p1, p2):
    """Vector from p1 to p2.
    :param p1: coordinates of point p1
    :param p2: coordinates of point p2
    :returns : numpy array with vector coordinates
    """
    return None if len(p1) != len(p2) else np.array(
        [p2[i] - p1[i] for i in range(len(p1))])


def vecangle(v1, v2, deg=True):
    """Calculate the angle between two vectors
    :param v1: coordinates of vector v1
    :param v2: coordinates of vector v2
    :param deg: whether to return degrees or radians
    :returns : angle in degree or rad
    """
    if np.array_equal(v1, v2):
        return 0.0
    dm = np.dot(v1, v2)
    cm = np.linalg.norm(v1) * np.linalg.norm(v2)
    angle = np.arccos(dm / cm)  # Round here to prevent floating point errors
    return np.degrees([angle, ])[0] if deg else angle


def centroid(coo):
    """Calculates the centroid from a 3D point cloud and returns the coordinates
    :param coo: Array of coordinate arrays
    :returns : centroid coordinates as list
    """
    return list(map(np.mean,
                    (([c[0] for c in coo]),
                     ([c[1] for c in coo]),
                        ([c[2] for c in coo]))))


def get_atom_coords(atom):
    mol = atom.GetOwningMol()
    conf = mol.GetConformers()[0]
    pos = conf.GetAtomPosition(atom.GetIdx())
    return (pos.x, pos.y, pos.z)


def get_coords(mol_conf, idxs):
    coords = []
    for idx in idxs:
        pos = mol_conf.GetAtomPosition(idx)
        coords.append((pos.x, pos.y, pos.z))
    return coords


def get_coord(mol_conf, idx):
    pos = mol_conf.GetAtomPosition(idx)
    return (pos.x, pos.y, pos.z)


def ring_is_planar(mol_conf, ring, r_atoms):
    """Given a set of ring atoms, check if the ring is sufficiently planar
    to be considered aromatic"""
    normals = []
    for a in r_atoms:
        a_pos = mol_conf.GetAtomPosition(a.GetIdx())
        a_coord = (a_pos.x, a_pos.y, a_pos.z)

        adj = a.GetNeighbors()
        # Check for neighboring atoms in the ring
        n_neighs_idx = [neigh.GetIdx()
                        for neigh in adj if neigh.GetIdx() in ring]
        n_coords = get_coords(mol_conf, n_neighs_idx)
        vec1, vec2 = vector(a_coord, n_coords[0]), vector(a_coord, n_coords[1])
        normals.append(np.cross(vec1, vec2))
    # Given all normals of ring atoms and their neighbors, the angle between
    # any has to be 5.0 deg or less
    for n1, n2 in itertools.product(normals, repeat=2):
        arom_angle = vecangle(n1, n2)
        if all([arom_angle > config.AROMATIC_PLANARITY,
                arom_angle < 180.0 - config.AROMATIC_PLANARITY]):
            return False
    return True


def residue_order(mol):
    residue_dict = {}
    for atom in mol.GetAtoms():
        res = atom.GetPDBResidueInfo()
        res_num = res.GetResidueNumber()
        res_name = res.GetResidueName()
        res_chain = res.GetChainId()
        key = str(res_num) + "_" + res_name + "_" + res_chain
        residue_dict.setdefault(key, [])
        residue_dict[key].append(atom)

    data = namedtuple(
        'residue',
        'residue_number residue_name residue_chain residue_atoms')
    residues = []
    for key, value in residue_dict.items():
        atoms = value
        residue_number = int(key.split("_")[0])
        residue_name = list(
            set([a.GetPDBResidueInfo().GetResidueName() for a in atoms]))
        residue_chain = key.split("_")[2]
        if len(residue_name) != 1:
            raise RuntimeError("get residue iterator error")
        residues.append(
            data(
                residue_number=residue_number,
                residue_name=residue_name[0],
                residue_chain=residue_chain,
                residue_atoms=atoms))
    return residues


def is_sidechain(atom):
    res = atom.GetPDBResidueInfo()
    atom_name = res.GetName().strip(" ")
    if atom_name in ("C", "CA", "N", "O", "H"):
        return False
    else:
        return True


def projection(pnormal1, ppoint, tpoint):
    """Calculates the centroid from a 3D point cloud and returns the coordinates
    :param pnormal1: normal of plane
    :param ppoint: coordinates of point in the plane
    :param tpoint: coordinates of point to be projected
    :returns : coordinates of point orthogonally projected on the plane
    """
    # Choose the plane normal pointing to the point to be projected
    pnormal2 = [coo * (-1) for coo in pnormal1]
    d1 = euclidean3d(tpoint, np.array(pnormal1) + np.array(ppoint))
    d2 = euclidean3d(tpoint, np.array(pnormal2) + np.array(ppoint))
    pnormal = pnormal1 if d1 < d2 else pnormal2
    # Calculate the projection of tpoint to the plane
    sn = -np.dot(pnormal, vector(ppoint, tpoint))
    sd = np.dot(pnormal, pnormal)
    sb = sn / sd
    return [c1 + c2 for c1, c2 in zip(tpoint, [sb * pn for pn in pnormal])]


def tilde_expansion(folder_path):
    """Tilde expansion, i.e. converts '~' in paths into <value of $HOME>."""
    return os.path.expanduser(
        folder_path) if '~' in folder_path else folder_path


def tmpfile(prefix, direc):
    """Returns the path to a newly created temporary file."""
    return tempfile.mktemp(prefix=prefix, suffix='.pdb', dir=direc)


def classify_by_name(names):
    """Classify a (composite) ligand by the HETID(s)"""
    if len(names) > 3:  # Polymer
        if len(set(config.RNA).intersection(set(names))) != 0:
            ligtype = 'RNA'
        elif len(set(config.DNA).intersection(set(names))) != 0:
            ligtype = 'DNA'
        else:
            ligtype = "POLYMER"
    else:
        ligtype = 'SMALLMOLECULE'

    for name in names:
        if name in config.METAL_IONS:
            if len(names) == 1:
                ligtype = 'ION'
            else:
                if "ION" not in ligtype:
                    ligtype += '+ION'
    return ligtype


def cluster_doubles(double_list):
    """Given a list of doubles, they are clustered if they share one element
    :param double_list: list of doubles
    :returns : list of clusters (tuples)
    """
    location = {}  # hashtable of which cluster each element is in
    clusters = []
    # Go through each double
    for t in double_list:
        a, b = t[0], t[1]
        # If they both are already in different clusters, merge the clusters
        if a in location and b in location:
            if location[a] != location[b]:
                if location[a] < location[b]:
                    clusters[location[a]] = clusters[location[a]].union(
                        clusters[location[b]])  # Merge clusters
                    clusters = clusters[:location[b]] + \
                        clusters[location[b] + 1:]
                else:
                    clusters[location[b]] = clusters[location[b]].union(
                        clusters[location[a]])  # Merge clusters
                    clusters = clusters[:location[a]] + \
                        clusters[location[a] + 1:]
                # Rebuild index of locations for each element as they have
                # changed now
                location = {}
                for i, cluster in enumerate(clusters):
                    for c in cluster:
                        location[c] = i
        else:
            # If a is already in a cluster, add b to that cluster
            if a in location:
                clusters[location[a]].add(b)
                location[b] = location[a]
            # If b is already in a cluster, add a to that cluster
            if b in location:
                clusters[location[b]].add(a)
                location[a] = location[b]
            # If neither a nor b is in any cluster, create a new one with a and
            # b
            if not (b in location and a in location):
                clusters.append(set(t))
                location[a] = len(clusters) - 1
                location[b] = len(clusters) - 1
    return map(tuple, clusters)


def is_lig(hetid):
    """Checks if a PDB compound can be excluded as a small molecule ligand"""
    h = hetid.upper()
    return not (h == 'HOH' or h in config.UNSUPPORTED)


def extract_pdbid(string):
    """Use regular expressions to get a PDB ID from a string"""
    p = re.compile("[0-9][0-9a-z]{3}")
    m = p.search(string.lower())
    try:
        return m.group()
    except AttributeError:
        return "UnknownProtein"


def read_pdb(pdbfname, as_string=False):
    """Reads a given PDB file and returns a Pybel Molecule."""
    pybel.ob.obErrorLog.StopLogging()  # Suppress all OpenBabel warnings
    if os.name != 'nt':  # Resource module not available for Windows
        maxsize = resource.getrlimit(resource.RLIMIT_STACK)[-1]
        resource.setrlimit(
            resource.RLIMIT_STACK, (min(
                2 ** 28, maxsize), maxsize))
    sys.setrecursionlimit(10 ** 5)  # increase Python recursion limit
    return readmol(pdbfname, as_string=as_string)


def create_folder_if_not_exists(folder_path):
    """Creates a folder if it does not exists."""
    folder_path = tilde_expansion(folder_path)
    folder_path = "".join(
        [folder_path, '/']) if not folder_path[-1] == '/' else folder_path
    direc = os.path.dirname(folder_path)
    if not folder_exists(direc):
        os.makedirs(direc)


def canonicalize(lig, preserve_bond_order=False):
    """Get the canonical atom order for the ligand."""
    atomorder = None
    # Get canonical atom order

    lig = pybel.ob.OBMol(lig.OBMol)
    if not preserve_bond_order:
        for bond in pybel.ob.OBMolBondIter(lig):
            if bond.GetBondOrder() != 1:
                bond.SetBondOrder(1)
    lig.DeleteData(pybel.ob.StereoData)
    lig = pybel.Molecule(lig)
    testcan = lig.write(format='can')
    try:
        pybel.readstring('can', testcan)
        reference = pybel.readstring('can', testcan)
    except IOError:
        testcan, reference = '', ''
    if testcan != '':
        reference.removeh()
        # isomorphs now holds all isomorphisms within the molecule
        isomorphs = get_isomorphisms(reference, lig)
        if not len(isomorphs) == 0:
            smi_dict = {}
            smi_to_can = isomorphs[0]
            for x in smi_to_can:
                smi_dict[int(x[1]) + 1] = int(x[0]) + 1
            atomorder = [smi_dict[x + 1] for x in range(len(lig.atoms))]
        else:
            atomorder = None
    return atomorder


def read(fil):
    """Returns a file handler and detects gzipped files."""
    if os.path.splitext(fil)[-1] == '.gz':
        return gzip.open(fil, 'rb')
    elif os.path.splitext(fil)[-1] == '.zip':
        zf = zipfile.ZipFile(fil, 'r')
        return zf.open(zf.infolist()[0].filename)
    else:
        return open(fil, 'r')


def nucleotide_linkage(residues):
    """Support for DNA/RNA ligands by finding missing covalent linkages to stitch DNA/RNA together."""

    nuc_covalent = []
    #######################################
    # Basic support for RNA/DNA as ligand #
    #######################################
    nucleotides = ['A', 'C', 'T', 'G', 'U', 'DA', 'DC', 'DT', 'DG', 'DU']
    dna_rna = {}  # Dictionary of DNA/RNA residues by chain
    covlinkage = namedtuple(
        "covlinkage",
        "id1 chain1 pos1 conf1 id2 chain2 pos2 conf2")
    # Create missing covlinkage entries for DNA/RNA
    for ligand in residues:
        resname, chain, pos = ligand
        if resname in nucleotides:
            if chain not in dna_rna:
                dna_rna[chain] = [(resname, pos), ]
            else:
                dna_rna[chain].append((resname, pos))
    for chain in dna_rna:
        nuc_list = dna_rna[chain]
        for i, nucleotide in enumerate(nuc_list):
            if not i == len(nuc_list) - 1:
                name, pos = nucleotide
                nextnucleotide = nuc_list[i + 1]
                nextname, nextpos = nextnucleotide
                newlink = covlinkage(
                    id1=name,
                    chain1=chain,
                    pos1=pos,
                    conf1='',
                    id2=nextname,
                    chain2=chain,
                    pos2=nextpos,
                    conf2='')
                nuc_covalent.append(newlink)

    return nuc_covalent


def sort_members_by_importance(members):
    """Sort the members of a composite ligand according to two criteria:
    1. Split up in main and ion group. Ion groups are located behind the main group.
    2. Within each group, sort by chain and position."""
    main = [x for x in members if x[0] not in config.METAL_IONS]
    ion = [x for x in members if x[0] in config.METAL_IONS]
    sorted_main = sorted(main, key=lambda x: (x[1], x[2]))
    sorted_main = sorted(main, key=lambda x: (x[1], x[2]))
    sorted_ion = sorted(ion, key=lambda x: (x[1], x[2]))
    return sorted_main + sorted_ion


def int32_to_negative(int32):
    """Checks if a suspicious number (e.g. ligand position) is in fact a negative number represented as a
    32 bit integer and returns the actual number.
    """
    dct = {}
    # Special case in some structures (note, this is just a workaround)
    if int32 == 4294967295:
        return -1
    for i in range(-1000, -1):
        dct[np.uint32(i)] = i
    if int32 in dct:
        return dct[int32]
    else:
        return int32


def is_donor(atom):
    if atom.GetSymbol() in ["H", "C"]:
        return False
    restype = whichrestype(atom)
    atomname = whichatomname(atom)
    # if restype in biolip_list:
    #    return False
    res_name = restype + "_" + atomname
    if res_name not in atom_prop_dict.keys():
        return False
    atom_prop = atom_prop_dict[res_name]
    if "Donor" in atom_prop:
        return True
    else:
        return False


def is_acceptor(atom):
    if atom.GetSymbol() in["H", "C"]:
        return False
    restype = whichrestype(atom)
    if restype == "HIN":
        restype = "HIS"
    atomname = whichatomname(atom)
    res_name = restype + "_" + atomname
    if res_name not in atom_prop_dict.keys():
        #print("warning atom: idx {} res_name {}".format(atom.GetIdx(), res_name))
        return False
    atom_prop = atom_prop_dict[res_name]
    if "Acceptor" in atom_prop:
        return True
    else:
        return False
