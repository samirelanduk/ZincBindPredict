from itertools import combinations
import numpy as np

def split_family(family):
    """Takes a family such as 'C3H1' and splits it into subfamilies such as 'C3'
    and 'H1'."""

    subfamilies, subfamily = [], ""
    for char in family:
        if char.isalpha() and subfamily:
            subfamilies.append(subfamily)
            subfamily = ""
        subfamily += char
    subfamilies.append(subfamily)
    return subfamilies


def model_to_residue_combos(model, family):
    """Takes an atomium model and returns all combinations of residues which
    match the family given."""

    subfamilies = split_family(family)
    residue_combos = []
    for subfamily in subfamilies:
        residues = model.residues(code=subfamily[0])
        residue_combos.append(tuple(combinations(residues, int(subfamily[1:]))))
    while len(residue_combos) > 1:
        residue_combos.insert(0, tuple(
         x + y for x in residue_combos[0] for y in residue_combos[1]
        ))
        for _ in range(2): residue_combos.pop(1)
    return residue_combos[0] if residue_combos else ()


def residues_to_sample(residues, site_id):
    """Converts a set of residues into a dict of values ready to be classified
    by the models."""

    sample = {}
    alphas, betas = [], []
    try:
        for res1, res2 in combinations(residues, 2):
            alphas.append(res1.atom(name="CA").distance_to(res2.atom(name="CA")))
            betas.append(res1.atom(name="CB").distance_to(res2.atom(name="CB")))
        sample["site"] = site_id
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
    except: return None