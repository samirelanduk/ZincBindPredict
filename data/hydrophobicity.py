from atomium import Model, Atom

partial_charges = {
 "ALA": {
  "C": 0.526,
  "CA": 0.215,
  "CB": 0.031,
  "HN": 0.248,
  "N": -0.52,
  "O": -0.5
 },
 "ARG": {
  "C": 0.526,
  "CA": 0.237,
  "CB": 0.049,
  "CD": 0.111,
  "CG": 0.058,
  "CZ": 0.813,
  "HN": 0.248,
  "HN11": 0.3615,
  "HN12": 0.3615,
  "HN21": 0.3615,
  "HN22": 0.3615,
  "HNE": 0.294,
  "N": -0.52,
  "NE": -0.493,
  "NH1": -0.6345,
  "NH2": -0.6345,
  "O": -0.5
 },
 "ASN": {
  "C": 0.526,
  "CA": 0.217,
  "CB": 0.003,
  "CG": 0.675,
  "HN": 0.248,
  "HND1": 0.344,
  "HND2": 0.344,
  "N": -0.52,
  "ND2": -0.867,
  "O": -0.5,
  "OD1": -0.47
 },
 "ASP": {
  "C": 0.526,
  "CA": 0.246,
  "CB": -0.208,
  "CG": 0.62,
  "HN": 0.248,
  "N": -0.52,
  "O": -0.5,
  "OD1": -0.706,
  "OD2": -0.706
 },
 "CYS": {
  "C": 0.526,
  "CA": 0.146,
  "CB": 0.1,
  "HN": 0.248,
  "HSG": 0.135,
  "LP1": -0.481,
  "LP2": -0.481,
  "N": -0.52,
  "O": -0.5,
  "SG": 0.827
 },
 "CYX": {
  "C": 0.526,
  "CA": 0.088,
  "CB": 0.143,
  "HN": 0.248,
  "LP1": -0.4045,
  "LP2": -0.4045,
  "N": -0.52,
  "O": -0.5,
  "SG": 0.824
 },
 "GLN": {
  "C": 0.526,
  "CA": 0.21,
  "CB": 0.053,
  "CD": 0.675,
  "CG": -0.043,
  "HN": 0.248,
  "HNE1": 0.344,
  "HNE2": 0.344,
  "N": -0.52,
  "NE2": -0.867,
  "O": -0.5,
  "OE1": -0.47
 },
 "GLU": {
  "C": 0.526,
  "CA": 0.246,
  "CB": 0.0,
  "CD": 0.62,
  "CG": -0.208,
  "HN": 0.248,
  "N": -0.52,
  "O": -0.5,
  "OE1": -0.706,
  "OE2": -0.706
 },
 "GLY": {
  "C": 0.526,
  "CA": 0.246,
  "HN": 0.248,
  "N": -0.52,
  "O": -0.5
 },
 "HIS": {
  "C": 0.526,
  "CA": 0.219,
  "CB": 0.06,
  "CD2": 0.145,
  "CE1": 0.384,
  "CG": 0.089,
  "HN": 0.248,
  "HND": 0.32,
  "N": -0.52,
  "ND1": -0.444,
  "NE2": -0.527,
  "O": -0.5
 },
 "HID": {
  "C": 0.526,
  "CA": 0.219,
  "CB": 0.06,
  "CD2": 0.145,
  "CE1": 0.384,
  "CG": 0.089,
  "HN": 0.248,
  "HND": 0.32,
  "N": -0.52,
  "ND1": -0.444,
  "NE2": -0.527,
  "O": -0.5
 },
 "HIE": {
  "C": 0.526,
  "CA": 0.219,
  "CB": 0.06,
  "CD2": 0.122,
  "CE1": 0.384,
  "CG": 0.112,
  "HN": 0.248,
  "HNE": 0.32,
  "N": -0.52,
  "ND1": -0.527,
  "NE2": -0.444,
  "O": -0.5
 },
 "HIP": {
  "C": 0.526,
  "CA": 0.195,
  "CB": 0.211,
  "CD2": 0.353,
  "CE1": 0.719,
  "CG": 0.103,
  "HN": 0.248,
  "HND": 0.478,
  "HNE": 0.486,
  "N": -0.52,
  "ND1": -0.613,
  "NE2": -0.686,
  "O": -0.5
 },
 "ILE": {
  "C": 0.526,
  "CA": 0.199,
  "CB": 0.03,
  "CD1": -0.001,
  "CG1": 0.017,
  "CG2": 0.001,
  "HN": 0.248,
  "N": -0.52,
  "O": -0.5
 },
 "LEU": {
  "C": 0.526,
  "CA": 0.204,
  "CB": 0.016,
  "CD1": -0.014,
  "CD2": -0.014,
  "CG": 0.054,
  "HN": 0.248,
  "N": -0.52,
  "O": -0.5
 },
 "LYS": {
  "C": 0.526,
  "CA": 0.227,
  "CB": 0.039,
  "CD": 0.048,
  "CE": 0.218,
  "CG": 0.053,
  "HN": 0.248,
  "HNZ1": 0.311,
  "HNZ2": 0.311,
  "HNZ3": 0.311,
  "N": -0.52,
  "NZ": -0.272,
  "O": -0.5
 },
 "MET": {
  "C": 0.526,
  "CA": 0.137,
  "CB": 0.037,
  "CE": 0.007,
  "CG": 0.09,
  "HN": 0.248,
  "LP1": -0.381,
  "LP2": -0.381,
  "N": -0.52,
  "O": -0.5,
  "SD": 0.737
 },
 "PHE": {
  "C": 0.526,
  "CA": 0.214,
  "CB": 0.038,
  "CD1": -0.011,
  "CD2": -0.011,
  "CE1": 0.004,
  "CE2": 0.004,
  "CG": 0.011,
  "CZ": -0.003,
  "HN": 0.248,
  "N": -0.52,
  "O": -0.5
 },
 "PRO": {
  "C": 0.526,
  "CA": 0.112,
  "CB": -0.001,
  "CD": 0.084,
  "CG": 0.036,
  "N": -0.257,
  "O": -0.5
 },
 "SER": {
  "C": 0.526,
  "CA": 0.292,
  "CB": 0.194,
  "HN": 0.248,
  "HOG": 0.31,
  "N": -0.52,
  "O": -0.5,
  "OG": -0.55
 },
 "THR": {
  "C": 0.526,
  "CA": 0.268,
  "CB": 0.211,
  "CG2": 0.007,
  "HN": 0.248,
  "HOG": 0.31,
  "N": -0.52,
  "O": -0.5,
  "OG1": -0.55
 },
 "TRP": {
  "C": 0.526,
  "CA": 0.248,
  "CB": 0.02,
  "CD1": 0.117,
  "CD2": -0.275,
  "CE2": 0.0,
  "CE3": 0.145,
  "CG": 0.046,
  "CH2": 0.034,
  "CZ2": 0.029,
  "CZ3": -0.082,
  "HN": 0.248,
  "HNE": 0.294,
  "N": -0.52,
  "NE1": -0.33,
  "O": -0.5
 },
 "TYR": {
  "C": 0.526,
  "CA": 0.245,
  "CB": 0.022,
  "CD1": -0.035,
  "CD2": -0.035,
  "CE1": 0.1,
  "CE2": 0.1,
  "CG": -0.001,
  "CZ": -0.121,
  "HN": 0.248,
  "HOH": 0.339,
  "N": -0.52,
  "O": -0.5,
  "OH": -0.368
 },
 "VAL": {
  "C": 0.526,
  "CA": 0.201,
  "CB": 0.033,
  "CG1": 0.006,
  "CG2": 0.006,
  "HN": 0.248,
  "N": -0.52,
  "O": -0.5
 }
}


def solvation(model, x, y, z, radius, pc=False, het=True, metal=True):
    """Determines the average solvation within a given sphere of an atomium
    model. By default, all atoms within the radius will be considered, but you
    can opt to exlcude heteroatoms (atoms not part of a chain residue) if you so
    desire.
    :param Model model: The atomium model to examine.
    :param x: The x-coordinate of the centre of the sphere.
    :param y: The y-coordinate of the centre of the sphere.
    :param z: The z-coordinate of the centre of the sphere.
    :param radius: The radius of the sphere.
    :param bool pc: If ``True``, atomic partial charges will be used instead of\
    atomic solvation parameters (squared).
    :param bool het: If ``False``, only atoms that have a residue will be\
    considered.
    :param bool metal: If ``False``, only non-metal atoms will be considered.
    :raises TypeError: if the model is not an atomium model object.
    :raises TypeError: if the coordinates are not numeric.
    :raises TypeError: if the radius is not numeric.
    :raises ValueError: if the radius is negative.
    :rtype: ``float``"""

    if not isinstance(model, Model):
        raise TypeError("{} is not a Model".format(model))
    if any(not isinstance(c, (int, float)) for c in (x, y, z)):
        raise TypeError("({}, {}, {}) not valid coordinate".format(x, y, z))
    if not isinstance(radius, (int, float)):
        raise TypeError("{} is not a valid radius".format(radius))
    if radius < 0:
        raise ValueError("{} is not a valid radius".format(radius))

    sphere = model.atoms_in_sphere((x, y, z), radius, is_metal=metal)
    solvations = ([atom_partial_charge(atom) ** 2 for atom in sphere]
     if pc else [atom_solvation(atom) for atom in sphere])
    return sum(solvations) / len(sphere) if len(solvations) else 0


def atom_solvation(atom):
    """Returns the atomic solvation parameter of an atomium atom. The atomic
    solvation parameters are taken from Yamashita et al (1990).
    :param Atom atom: an atomium atom object.
    :rtype: ``float``"""

    specials = {
     "O": {"GLU": ["OE1", "OE2"], "ASP": ["OD1", "OD2"]},
     "N": {"HIS": ["ND1", "NE2"], "ARG": ["NH1", "NH2"]}
    }
    if atom.element == "C": return 18
    if atom.element == "S": return -5
    if atom.element in specials:
        if atom.charge != 0:
            return -37 if atom.element == "O" else -38
        if atom.het and atom.het.name in specials[atom.element]:
            if atom.name in specials[atom.element][atom.het.name]:
                return -23 if atom.element == "O" else -23.5
        return -9
    return 0


def atom_partial_charge(atom):
    """Returns the atomic partial charge of an atomium atom.
    :param Atom atom: an atomium atom object.
    :rtype: ``float``"""

    if atom.charge != 0: return atom.charge
    if atom.residue is not None and atom.residue.name in partial_charges:
        if atom.name in partial_charges[atom.residue.name]:
            return partial_charges[atom.residue.name][atom.name]
    return 0


def hydrophobic_contrast(model, x, y, z, radius, pc=False, het=True, metal=True):
    """Determines the hydrophobic contrast within a sphere - a measure of
    how heterogenous the hydrophobicity is within the sphere.
    A homogeneous sphere will evaluate to zero, a sphere with a region of high
    hydrophilic atoms enclosed within a region of high hydrophobic regions will
    have a high positive value, and the converse will have a high negative
    value.
    :param Model model: The atomium model to examine.
    :param x: The x-coordinate of the centre of the sphere.
    :param y: The y-coordinate of the centre of the sphere.
    :param z: The z-coordinate of the centre of the sphere.
    :param radius: The radius of the sphere.
    :param bool het: If ``False``, only atoms that have a residue will be\
    considered.
    :param bool metal: If ``False``, only non-metal atoms will be considered.
    :raises TypeError: if the model is not an atomium model object.
    :raises TypeError: if the coordinates are not numeric.
    :raises TypeError: if the radius is not numeric.
    :raises ValueError: if the radius is negative.
    :rtype: ``float``"""

    

    if not isinstance(model, Model):
        raise TypeError("{} is not a Model".format(model))
    if any(not isinstance(c, (int, float)) for c in (x, y, z)):
        raise TypeError("({}, {}, {}) not valid coordinate".format(x, y, z))

    if not isinstance(radius, (int, float)):
        raise TypeError("{} is not a valid radius".format(radius))
    
    if radius < 0:
        raise ValueError("{} is not a valid radius".format(radius))

    sphere = model.atoms_in_sphere((x, y, z), radius, is_metal=metal)
    if len(sphere) == 0: return 0
    average_solvation = solvation(model, x, y, z, radius, pc=pc, het=het, metal=metal)
    sum_, r2 = 0, 0
    for atom in sphere:
        distance = atom.distance_to((x, y, z))
        solv = ((atom_partial_charge(atom)) ** 2) if pc else atom_solvation(atom)
        sum_ += solv * (distance ** 2)
        r2 += (distance ** 2)
    r2 /= len(sphere)
    return sum_ - (len(sphere) * average_solvation * r2)