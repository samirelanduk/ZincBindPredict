"""Functions used in the three zincbindpredict packages - data, predict and
server."""

import random

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



   