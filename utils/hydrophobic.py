from collections import namedtuple

radius = {"N": 1.8, "O": 1.7, "S": 2.0, "P": 2.1, "F": 1.5, "Cl": 1.8,
          "Br": 2.0, "I": 2.2, "C": 1.9, "H": 0.0, "Zn": 0.5, "B": 1.8,
          "Si": 1.8, "As": 1.8, "Se": 1.8}

def get_hyd_strength(dist, patom, latom):
    p_radii = radius.get(patom.GetSymbol(), 0.5)
    l_radii = radius.get(latom.GetSymbol(), 0.5)

    d1 = p_radii + l_radii
    if dist <= d1 + 0.5:
        return -1.0
    elif d1 + 0.5 < dist < d1 + 2.0:
        return -0.66666 * (d1 + 2.0 - dist)
    else:
        return 0

def calc_hydrophobic(hydro):
    dist = hydro.distance
    patom = hydro.bsatom
    latom = hydro.ligatom
    
    hydrophobic_energy = get_hyd_strength(dist, patom, latom)
    return hydrophobic_energy 

