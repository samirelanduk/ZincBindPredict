import math
import random
import os
from itertools import combinations, product
from functools import reduce
import numpy as np
import kirjava

def fetch_data(url, query, variables):
    """Gets data from the ZincBindDB API and formats it to remove all the edges
    and nodes and stuff."""

    def strip_dict(d):
        for key, value in d.items():
            if isinstance(value, dict) and "edges" in value:
                d[key] = [strip_dict(e["node"]) for e in value["edges"]]
        return d
    data = kirjava.execute(url, query, variables=variables)["data"]
    strip_dict(data)
    return list(data.values())[0]


def update_data_file(family, samples=None):
    """Saves samples to a CSV file if given, otherwise it just makes an empty
    file ready for samples to be saved to it."""

    path = f"data/{family}.csv"
    if samples:
        csv_lines = [",".join(samples[0].keys()) + "\n"]\
         if os.path.getsize(path) == 0 else []
        for sample in samples:
            csv_lines.append(",".join([str(v) for v in sample.values()]) + "\n")
        with open(path, "a") as f: f.writelines(csv_lines)
    else:
        with open(path, "w") as f: f.write("")


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


def get_residues_from_model(model, residues):
    """Takes a model, and some residues in JSON form, and then finds the
    matching residues in the model.
    
    It will throw an exception if it encounters any problems."""

    at_residues = []
    for res in residues:
        possibles = model.residues(id=res["atomiumId"])
        matching = [r for r in possibles if r.atom(name="CA").location == (
         res["atoms"][0]["x"], res["atoms"][0]["y"], res["atoms"][0]["z"]
        )]
        at_residues.append(matching[0])
    return at_residues


def count_combinations(model, family):
    """Takes a model and a family string, and returns the number of matching
    residue combinations in that model."""

    subfamilies = split_family(family)
    counts = [[len(model.residues(code=f[0])), int(f[1:])] for f in subfamilies]
    combo_counts = [(math.factorial(c[0]) / (
     math.factorial(c[1]) * math.factorial(c[0] - c[1]
    ))) if c[0] - c[1] >= 0 else 0 for c in counts]
    count = reduce(lambda x, y: x * y, combo_counts) if combo_counts else 0
    return int(count)


def model_to_residue_combos(model, family, limit=None, ignore=None):
    """Takes an atomium model and returns all combinations of residues which
    match the family given. You can limit the number of combinations it returns
    to prevent catastrophic consumption of residues if you want."""

    combo_count = count_combinations(model, family)
    indices = random.sample(range(combo_count), limit)\
     if limit and limit <= combo_count else range(combo_count)
    subfamilies = split_family(family)
    residue_combos = []
    for subfamily in subfamilies:
        residues = model.residues(code=subfamily[0])
        residue_combos.append(combinations(residues, int(subfamily[1:])))
    initial_combos = product(*residue_combos)
    count, combos = 0, []
    for combo in initial_combos:
        if count in indices:
            combo = tuple([item for sublist in combo for item in sublist])
            ids = set([res.id for res in combo])
            if not ignore or ids not in ignore:
                combos.append(combo)
        count += 1
    return  combos
    

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
    except Exception as e: return None


def sequence_to_residue_combos(sequence, family, limit=None):
    """Takes a sequence and returns all combinations of residues which
    match the family given. You can limit the number of combinations it returns
    to prevent catastrophic consumption of residues if you want."""

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
    return ["".join([char.upper() if i in combo else char for i, char
     in enumerate(sequence.lower())]) for combo in combos]


def sequence_to_sample(sequence, site_id):
    """Converts a sequence into a dict of values ready to be classified
    by the models."""

    sample = {}
    sample["site"] = site_id
    caps = []
    for index, char in enumerate(sequence):
        if char.isupper(): caps.append(index)
    spacers = [(second - first) - 1 for first, second in zip(caps[:-1], caps[1:])]
    sample["min_spacer"], sample["max_spacer"] = min(spacers), max(spacers)
    for i, spacer in enumerate(spacers, start=1):
        sample[f"spacer_{i}"] = spacer
    return sample