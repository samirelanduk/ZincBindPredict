"""Functions used in the three zincbindpredict packages - data, predict and
server."""

import random
from itertools import combinations
import math
import numpy as np
import sys
sys.path.append("../zincbindpredict")
from data.hydrophobicity import *

def sequence_site_to_sample(sequence):
    """Takes a sequence site and turns it into a feature sample. The site will
    be a string where the binding residues are upper case and everything else is
    in lower case."""

    site = {}
    residues = [i for i, char in enumerate(sequence) if char.isupper()]
    for i, residues in enumerate(zip(residues[:-1], residues[1:]), start=1):
        site[f"gap{i}"] = residues[1] - residues[0] - 1
        site[f"hydrophobicity{i}"] = average_hydrophobicity(
            sequence, span=[residues[0], residues[1]]
        )
    site["hydrophobicity_window_1"] = average_hydrophobicity(sequence, window=1)
    site["hydrophobicity_window_3"] = average_hydrophobicity(sequence, window=3)
    site["hydrophobicity_window_5"] = average_hydrophobicity(sequence, window=5)
    site["charged_window_1"] = residue_count(sequence, "DERHK", window=1)
    site["charged_window_3"] = residue_count(sequence, "DERHK", window=3)
    site["charged_window_5"] = residue_count(sequence, "DERHK", window=5)
    for res in "ARNDCEQGHILKMFPSTWYV":
        site[f"{res}_window_1"] = residue_count(sequence, res, window=1)
        site[f"{res}_window_3"] = residue_count(sequence, res, window=3)
        site[f"{res}_window_5"] = residue_count(sequence, res, window=5)
    return site


def structure_site_to_sample(residues):
    """Converts a set of residues into a dict of values ready to be classified
    by the models."""

    sample = {}
    alphas, betas = [], []
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
    sample["hcf"] = hydrophobic_contrast_function(residues)
    return sample


def residue_count(sequence, residues, window=1, span=None):
    """Counts the number of residues around the binding residues in a sequence
    that match a list of residus to match."""
    
    counts = 0
    total_window = 0
    for index, char in enumerate(sequence):
        if char.isupper():
            sub_sequence = sequence[index - window: index + window + 1]
            for i, sub_char in enumerate(sub_sequence):
                    if i != (len(sub_sequence) - 1) / 2:
                        total_window += 1
                        if sub_char.upper() in residues.upper(): counts += 1
    return round(counts / total_window, 3) if total_window else 0