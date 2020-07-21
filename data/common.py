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



   