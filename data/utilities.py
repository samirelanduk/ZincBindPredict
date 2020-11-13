"""Utility functions used in the generation of datasets."""

import random
import os
import kirjava
from collections import Counter
import sys
sys.path.append("../zincbindpredict")
from data.common import split_family, structure_family_site_to_vector

def fetch_data(url, query, variables):
    """Gets data from the ZincBindDB API and formats it to remove all the edges
    and nodes and stuff."""

    def strip_dict(d):
        for key, value in d.items():
            if isinstance(value, dict) and "edges" in value:
                d[key] = [strip_dict(e["node"]) for e in value["edges"]]
            elif isinstance(value, dict):
                d[key] = strip_dict(d[key])
        return d
    data = kirjava.execute(url, query, variables=variables)["data"]
    strip_dict(data)
    return list(data.values())[0]


def save_csv(positives, negatives, name, path):
    """Takes a list of positive samples and a list of negative samples, and
    saves them to CSV."""

    if not positives and not negatives: return

    path = f"{path}{os.path.sep}{name}.csv"
    lines = []
    flag = "a+"
    if not os.path.exists(path):
        lines.append(",".join((positives + negatives)[0].keys()) + ",positive")
        #flag = "w"
    for index, samples in enumerate([positives, negatives]):
        for sample in samples:
            lines.append(",".join([str(v) for v in sample.values()] + [str(1 - index)]))
    with open(path, flag) as f:
        f.write("\n".join(lines) + "\n")


def random_sequence_family_input(sequence, family):
    """Takes a sequence string and gets a random binding site for a given
    family. If there is no matching site, returns None."""

    indices = []
    for subfamily in split_family(family):
        code, count = subfamily[0].lower(), subfamily[1]
        subfamily_indices = [index for index, char in
            enumerate(sequence.lower()) if char == code]
        if len(subfamily_indices) < count: return None
        indices += random.sample(subfamily_indices, count)
    return "".join([char.upper() if index in indices else char.lower()
        for index, char in enumerate(sequence)])


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


def create_negative_samples_for_model(model, family, positives):
    """Takes a model, and generates a number of residue combinations for a given
    family, excluding those flagges as positive."""

    negative_samples = []
    iter_count = 0
    while len(negative_samples) < len(positives) * 1:
        random_combination = random_structure_family_input(model, family)
        if not random_combination: break
        if set(random_combination) not in positives:
            vector = structure_family_site_to_vector(random_combination)
            if vector:
                negative_samples.append(vector)
        iter_count += 1
        if iter_count > 100 * len(positives): break
    return negative_samples


def random_structure_family_input(model, family):
    """Takes an atomium model and generates a random combination of residues
    matching a family, or None if that is impossible."""
    
    residues = []
    for subfamily in split_family(family):
        code, count = subfamily[0].upper(), subfamily[1]
        matching_residues = list(model.residues(code=code))
        if len(matching_residues) < count: return None
        residues += random.sample(matching_residues, count)
    return residues


def chars_to_family(chars):
    """Takes a list of characters and constructs a family from them. So, A1B2
    would be created from ['B', 'A', 'B'] for example."""

    counter = Counter(chars)
    return "".join(sorted([char + str(n) for char, n in counter.items()]))
    





















def update_data_file(family, kind, samples=None):
    """Saves samples to a CSV file if given, otherwise it just makes an empty
    file ready for samples to be saved to it."""

    path = f"data/csv/{kind}/{family}.csv"
    if samples:
        csv_lines = [",".join(samples[0].keys()) + "\n"]\
         if os.path.getsize(path) == 0 else []
        for sample in samples:
            csv_lines.append(",".join([str(v) for v in sample.values()]) + "\n")
        with open(path, "a") as f: f.writelines(csv_lines)
    else:
        with open(path, "w") as f: f.write("")





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
        
        '''stabiliser_contacts = set()
        hydrogen_bonds = set()
        for residue in residues:
            for atom in residue.atoms():
                nearby = atom.nearby_atoms(3.5)
                for nearby_atom in nearby:
                    if isinstance(nearby_atom.het, atomium.Residue)\
                     and nearby_atom.het is not residue:
                        stabiliser_contacts.add((atom, nearby_atom))
                        if atom.element in ["N", "O", "CL"] and nearby_atom.element in ["N", "O", "CL"]:
                            hydrogen_bonds.add((atom, nearby_atom))

        sample["contacts"] = len(stabiliser_contacts)
        sample["h_bonds"] = len(hydrogen_bonds)'''
        return sample
    except Exception as e: return None


def model_to_residue_combos(model, family, ignore=None):
    """Takes an atomium model and returns all combinations of residues which
    match the family given. You can limit the number of combinations it returns
    to prevent catastrophic consumption of residues if you want."""
    
    subfamilies = split_family(family)
    residues = []
    for subfamily in subfamilies:
        residues += list(model.residues(code=subfamily[0]))
    residues = {r: [r] for r in residues}

    for res in residues:
        for other_res in residues:
            if res is not other_res and res.atom(name="CA") and\
             other_res.atom(name="CA") and\
              res.atom(name="CA").distance_to(other_res.atom(name="CA")) <= 25:
                residues[res].append(other_res)
    combos = set()
    for res in residues:
        residue_combos = []
        for subfamily in subfamilies:
            sub_res = [r for r in residues[res] if r.code == subfamily[0]]
            residue_combos.append(combinations(sub_res, int(subfamily[1:])))
        initial_combos = product(*residue_combos)
        for combo in initial_combos:
            combo = frozenset([item for sublist in combo for item in sublist])
            ids = set([res.id for res in combo])
            if not ignore or ids not in ignore:
                combos.add(combo)
    return combos


def count_model_combinations(model, family):
    """Takes a model and a family string, and returns the number of matching
    residue combinations in that model."""

    subfamilies = split_family(family)
    counts = [[len(model.residues(code=f[0])), int(f[1:])] for f in subfamilies]
    combo_counts = [(math.factorial(c[0]) / (
     math.factorial(c[1]) * math.factorial(c[0] - c[1]
    ))) if c[0] - c[1] >= 0 else 0 for c in counts]
    count = reduce(lambda x, y: x * y, combo_counts) if combo_counts else 0
    return int(count)




def split_dataset(df):
    """Takes a pandas dataset and produces four variants on it. They all have
    the last column removed. The second one is limited to positive cases,
    the third to negative cases, and the fourth is unlabelled *and* has no IDs.
    """

    unlabelled = df.iloc[:, :-1]
    positives = df.loc[df["positive"] == 1].iloc[:, :-1]
    negatives = df.loc[df["positive"] == -1].iloc[:, :-1]
    core = df.iloc[:, 1:-1]
    return unlabelled, positives, negatives, core


def family_count(family):
    return sum([int(sub[1:]) for sub in split_family(family)])


def residue_sequence_count(sequence):
    return len([char for char in sequence if char.isupper()])


def sequence_to_sample(sequence, id):
    sample = {}
    residues = [[i, char] for i, char in enumerate(sequence) if char.isupper()]
    sample["id"] = id
    for i, residues in enumerate(zip(residues[:-1], residues[1:]), start=1):
        sample[f"gap{i}"] = residues[1][0] - residues[0][0]
    return sample


def sequence_to_residue_combos(sequence, family, limit=None, ignore=None):
    sequence = sequence.lower()
    combo_count = count_sequence_combinations(sequence, family)
    indices = random.sample(range(combo_count), limit)\
     if limit and limit <= combo_count else range(combo_count)
    subfamilies = split_family(family)
    residue_combos = []
    for subfamily in subfamilies:
        residues = [i for i, char in enumerate(sequence)
         if char == subfamily[0].lower()]
        residue_combos.append(combinations(residues, int(subfamily[1:])))
    initial_combos = product(*residue_combos)
    count, combos = 0, []
    for combo in initial_combos:
        if count in indices:
            combo = tuple([item for sublist in combo for item in sublist])
            if not ignore or combo not in ignore:
                combos.append(combo)
        count += 1
    return ["".join([
     c.upper() if i in combo else c for i, c in enumerate(sequence)
    ]) for combo in combos]


def count_sequence_combinations(sequence, family):
    """Takes a sequence and a family string, and returns the number of matching
    residue combinations in that sequence."""

    subfamilies = split_family(family)
    counts = [[len([
     c for c in sequence.upper() if c == f[0]
    ]), int(f[1:])] for f in subfamilies]
    combo_counts = [(math.factorial(c[0]) / (
     math.factorial(c[1]) * math.factorial(c[0] - c[1]
    ))) if c[0] - c[1] >= 0 else 0 for c in counts]
    count = reduce(lambda x, y: x * y, combo_counts) if combo_counts else 0
    return int(count)