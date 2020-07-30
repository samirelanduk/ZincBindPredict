"""Functions used in the three zincbindpredict packages - data, predict and
server."""

import random
from itertools import combinations
import math
import numpy as np

def sequence_site_to_vector(sequence):
    """Takes a sequence site and turns it into a feature vector. The site will
    be a string where the binding residues are upper case and everything else is
    in lower case."""

    site = {}
    residues = [i for i, char in enumerate(sequence) if char.isupper()]
    for i, residues in enumerate(zip(residues[:-1], residues[1:]), start=1):
        site[f"gap{i}"] = residues[1] - residues[0] - 1
    site["hydrophobicity"] = average_hydrophobicity(sequence)
    return site


def split_family(family):
    """Takes a family such as 'C3H1' and splits it into subfamilies such as 'C3'
    and 'H1'."""

    subfamilies, subfamily = [], ""
    for char in family:
        if char.isalpha() and subfamily:
            subfamilies.append([subfamily[0], int(subfamily[1:])])
            subfamily = ""
        subfamily += char
    subfamilies.append([subfamily[0], int(subfamily[1:])])
    return subfamilies


def average_hydrophobicity(sequence):
    """Takes a sequence, looks at the residues on either side of the upper case
    binding residues, and works out the average of their hydrophobicities."""

    scale = {
        "A": 0.17, "R": 0.81, "N": 0.42, "D": 1.23, "C": -0.24, "E": 2.02,
        "Q": 0.58, "G": 0.01, "H": 0.96, "I": -0.31, "L": -0.56, "K": 0.99,
        "M": -0.23, "F": -1.13, "P": 0.45, "S": 0.13, "T": 0.14, "W": -1.85,
        "Y": -0.94, "V": 0.07
    }
    scores = [score for sublist in [
        [
            scale.get(sequence[index - 1].upper()) if index != 0 else None,
            scale.get(sequence[index + 1].upper()) if\
             index != len(sequence) - 1 else None
        ] for index, char in enumerate(sequence) if char.isupper()
    ] for score in sublist if score is not None]
    return round(sum(scores) / len(scores), 3)


def structure_family_site_to_vector(residues):
    """Converts a set of residues into a dict of values ready to be classified
    by the models."""

    sample = {}
    alphas, betas = [], []
    try:
        for res1, res2 in combinations(residues, 2):
            alphas.append(res1.atom(name="CA").distance_to(res2.atom(name="CA")))
            betas.append(res1.atom(name="CB").distance_to(res2.atom(name="CB")))
        sample["ca_mean"] = round(sum(alphas) / len(alphas), 3)
        sample["ca_std"] = round(np.std(alphas), 3)
        sample["ca_min"] = round(min(alphas), 3)
        sample["ca_max"] = round(max(alphas), 3)
        sample["cb_mean"] = round(sum(betas) / len(betas), 3)
        sample["cb_std"] = round(np.std(betas), 3)
        sample["cb_min"] = round(min(betas), 3)
        sample["cb_max"] = round(max(betas), 3)
        sample["helix"] = len([r for r in residues if r.helix])
        sample["strand"] = len([r for r in residues if r.strand])
        return sample
    except Exception as e: return None


def structure_half_family_site_to_vector(residues):
    """Converts a set of residues representing a half-site into a dict of values
    ready to be classified by the models."""

    sample = {}
    alphas, betas = [], []
    try:
        for res1, res2 in combinations(residues, 2):
            alphas.append(res1.atom(name="CA").distance_to(res2.atom(name="CA")))
            betas.append(res1.atom(name="CB").distance_to(res2.atom(name="CB")))
        sample["ca_mean"] = round(sum(alphas) / len(alphas), 3)
        sample["ca_std"] = round(np.std(alphas), 3)
        sample["ca_min"] = round(min(alphas), 3)
        sample["ca_max"] = round(max(alphas), 3)
        sample["cb_mean"] = round(sum(betas) / len(betas), 3)
        sample["cb_std"] = round(np.std(betas), 3)
        sample["cb_min"] = round(min(betas), 3)
        sample["cb_max"] = round(max(betas), 3)
        sample["helix"] = len([r for r in residues if r.helix])
        sample["strand"] = len([r for r in residues if r.strand])
        return sample
    except Exception as e: return None


def model_to_grid(model):
    """Takes a model and returns a grid of values around the electronegative
    atoms."""

    grid = set()
    for residue in model.residues():
        for atom in residue.atoms(element__regex="N|O|S"):
            location = tuple(round(n) for n in atom.location)
            for n in range(3):
                for modifier in [-3, 3]:
                    altered_location = list(location)
                    altered_location[n] += modifier
                    grid.add(tuple(altered_location))
    return tuple(grid)


def location_to_vector(point, model):
    """Converts a location in space within a model into a dict of values ready
    to be classified by the models."""
    
    vector = {}
    nearby_atoms = model.atoms_in_sphere(point, 8, element__ne="ZN")
    vector["8_atom_count"] = count = len(nearby_atoms)
    vector["8_ched_ratio"] = round(len([
        a for a in nearby_atoms if a.element in "SNO"
    ]) / count, 2) if count else 0
    vector["center_offset"] = round(math.sqrt(
        (((sum(a.location[0] for a in nearby_atoms) / count) - point[0]) ** 2) +
        (((sum(a.location[1] for a in nearby_atoms) / count) - point[1]) ** 2) +
        (((sum(a.location[2] for a in nearby_atoms) / count) - point[2]) ** 2)
    ), 2) if count else 0
    distant_atoms = model.atoms_in_sphere(point, 16, element__ne="ZN")
    vector["16_atom_count"] = count = len(distant_atoms)
    vector["16_ched_ratio"] = round(len([
        a for a in distant_atoms if a.element in "SNO"
    ]) / count, 2) if count else 0
    return vector


def half_location_to_vector(point, model, chain=None):
    """Converts a location in space within a model into a dict of values ready
    to be classified by the models - only the chain specified will be examined
    for surrounding atoms."""
    
    vector = {}
    nearby_atoms = model.atoms_in_sphere(
        point, 8, element__ne="ZN", chain__id=chain
    ) if chain else model.atoms_in_sphere(point, 8, element__ne="ZN")
    vector["8_atom_count"] = count = len(nearby_atoms)
    vector["8_ched_ratio"] = round(len([
        a for a in nearby_atoms if a.element in "SNO"
    ]) / count, 2) if count else 0
    vector["center_offset"] = round(math.sqrt(
        (((sum(a.location[0] for a in nearby_atoms) / count) - point[0]) ** 2) +
        (((sum(a.location[1] for a in nearby_atoms) / count) - point[1]) ** 2) +
        (((sum(a.location[2] for a in nearby_atoms) / count) - point[2]) ** 2)
    ), 2) if count else 0
    distant_atoms = model.atoms_in_sphere(
        point, 16, element__ne="ZN", chain__id=chain
    ) if chain else model.atoms_in_sphere(point, 16, element__ne="ZN")
    vector["16_atom_count"] = count = len(distant_atoms)
    vector["16_ched_ratio"] = round(len([
        a for a in distant_atoms if a.element in "SNO"
    ]) / count, 2) if count else 0
    return vector

