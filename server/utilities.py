"""Utilility functions for enabling API access to the models."""

import sys
sys.path.append(".")
from itertools import combinations, product
import requests
import os
import json
from django.http import JsonResponse
from data.utilities import model_to_residue_combos, residues_to_sample
from data.utilities import sequence_to_residue_combos, sequence_to_sample

def parse_arguments():
    """Gets the JSON arguments passed in at the command line."""

    if len(sys.argv) < 2:
        print("Please provide JSON arguments")
        sys.exit(1)
    return json.loads(sys.argv[1])


def get_job_location(id):
    """Gets the location of a job file on disk."""

    return f"server{os.path.sep}jobs{os.path.sep}{id}.json"


def load_job(id):
    """Loads a job dictionary from JSON file."""

    with open(get_job_location(id)) as f:
        return json.load(f)


def save_job(job, status=None):
    """Saves a job dictionary to file as JSON."""

    if status: job["status"] = status
    with open(get_job_location(job["id"]), "w") as f:
        json.dump(job, f, indent=4)


def initialise_job(id, type, protein):
    """Creates an empty job dictionary."""

    return {
        "id": id,
        "status": "initialising",
        "type": type,
        "protein": protein,
        "sites": [], "rejected": []
    }


def get_sequence_families():
    """Get the families for which there are sequence models."""

    return ["H3", "C4", "C2H2"]







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


def model_to_family_inputs(model, family):
    """Takes a model and returns a list of family inputs."""

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
            
            combos.add(combo)
    return combos

    
def structure_family_site_to_vector(site):
    """Takes a structure family site object and turns it into a feature vector."""

    return []


def sequence_to_family_inputs(sequence, family):
    """Takes a sequence and returns a list of potential binding sites for a
    given family."""

    sequence, family = sequence.lower(), family.lower()
    subfamilies = split_family(family.lower())
    subfamily_combinations = []
    for subfamily in subfamilies:
        print(subfamily)
        subfamily_indices = [i for i, char in enumerate(sequence)
            if char == subfamily[0]]
        print(subfamily_indices)
        subfamily_combinations.append(
            list(combinations(subfamily_indices, subfamily[1]))
        )
    potential_sites = [[
        index for subsite in site for index in subsite
    ] for site in product(*subfamily_combinations)]
    return ["".join(
        char.upper() if i in site else char for i, char in enumerate(sequence)
    ) for site in potential_sites]


def sequence_site_to_vector(site):
    """Takes a sequence site object and turns it into a feature vector."""

    return []


def sequence_site_to_gql_object(site, family, probability):
    """Creates a sequence site object ready for the graphql server, from the
    sequence input object and its attributes."""
    
    residues = [{
        "name": sequence[index], "identifier": index
    } for index in site]
    return {family: family, probability: probability, residues: residues}