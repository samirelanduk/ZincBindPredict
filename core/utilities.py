from itertools import combinations, product
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


def model_to_residue_combos(model, family, limit=None):
    """Takes an atomium model and returns all combinations of residues which
    match the family given. You can limit the number of combinations it returns
    to prevent catastrophic consumption of residues if you want."""

    subfamilies = split_family(family)
    residue_combos = []
    for subfamily in subfamilies:
        residues = model.residues(code=subfamily[0])
        residue_combos.append(combinations(residues, int(subfamily[1:])))
    initial_combos = product(*residue_combos)
    count, combos = 1, []
    for combo in initial_combos:
        combos.append(tuple([item for sublist in combo for item in sublist]))
        if limit and count == limit: break
        count += 1
    return tuple(combos)


def sequence_to_residue_combos(sequence, family, limit=None):
    """"""

    subfamilies = split_family(family)
    residue_combos = []
    for subfamily in subfamilies:
        indices = [i for i, r in enumerate(sequence.upper()) if r == subfamily[0]]
        residue_combos.append(combinations(indices, int(subfamily[1:])))
    initial_combos = product(*residue_combos)
    count, combos = 1, []
    for combo in initial_combos:
        combos.append(tuple([item for sublist in combo for item in sublist]))
        if limit and count == limit: break
        count += 1
    return ["".join([char.upper() if i in combo else char for i, char in enumerate(sequence.lower())]) for combo in combos]
    

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


def sequence_to_sample(sequence, site_id):
    sample = {}
    sample["site"] = site_id
    caps = []
    for index, char in enumerate(sequence):
        if char.isupper(): caps.append(index)
    
    spacers = [second - first for first, second in zip(caps[:-1], caps[1:])]
    sample["min_spacer"], sample["max_spacer"] = min(spacers), max(spacers)
    for i, spacer in enumerate(spacers, start=1):
        sample[f"spacer_{i}"] = spacer
    return sample