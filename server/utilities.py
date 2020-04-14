"""Utilility functions for enabling API access to the models."""

import sys
sys.path.append(".")
import time
from itertools import combinations, product
import requests
import atomium
import os
import json
from django.http import JsonResponse
from data.utilities import model_to_residue_combos, residues_to_sample
from data.utilities import sequence_to_residue_combos, sequence_to_sample

def initialize_job(protein):
    """Creates an empty job dictionary with an ID created from the current
    time."""

    return {
        "id": str(int(time.time() * 1000)),
        "status": "initializing",
        "protein": protein,
        "sites": [], "rejected": []
    }


def get_job_location(id):
    """Gets the location of a job file on disk."""

    return f"server{os.path.sep}jobs{os.path.sep}{id}.json"


def save_job(job, status=None):
    """Saves a job dictionary to file as JSON."""

    if status: job["status"] = status
    with open(get_job_location(job["id"]), "w") as f:
        json.dump(job, f, indent=4)


def save_structure_file(uploaded_file, job_id):
    """Takes an uploaded file and a job ID, and saves the uploaded locally. The
    new file name will be returned."""

    name = uploaded_file.name
    file_extension = ("." + name.split(".")[-1]) if "." in name else ""
    file_name = f"{job_id}{file_extension}"
    file_path = f"server{os.path.sep}jobs{os.path.sep}{file_name}"
    with open(file_path, "wb") as f:
        f.write(uploaded_file.read())
    return file_name




def parse_arguments():
    """Gets the JSON arguments passed in at the command line."""

    if len(sys.argv) < 2:
        print("Please provide JSON arguments")
        sys.exit(1)
    return json.loads(sys.argv[1])


def load_job(id):
    """Loads a job dictionary from JSON file."""

    with open(get_job_location(id)) as f:
        return json.load(f)


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




def get_sequence_families():
    """Get the families for which there are sequence models."""

    return ["H3", "C4", "C2H2"]


def sequence_to_family_inputs(sequence, family):
    """Takes a sequence and returns a list of potential binding sites for a
    given family."""

    sequence, family = sequence.lower(), family.lower()
    subfamilies = split_family(family)
    subfamily_combinations = []
    for subfamily in subfamilies:
        subfamily_indices = [i for i, char in enumerate(sequence)
            if char == subfamily[0]]
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




def get_model_for_job(filename):
    """Gets an atomium model from a given filename representing a local
    structure file."""

    return atomium.open(f"server{os.path.sep}jobs{os.path.sep}{filename}").model


def get_structure_families():
    """Get the families for which there are structure models."""

    return ["H3", "C4", "C2H2"]


def model_to_family_inputs(model, family):
    """Takes a model and returns a list of potential binding sites for a
    given family."""

    family = family.lower()
    subfamilies = split_family(family)
    subfamily_combinations = []
    for subfamily in subfamilies:
        residues = model.residues(code=subfamily[0].upper())
        subfamily_combinations.append(
            list(combinations(residues, subfamily[1]))
        )
    potential_sites = [[
        residue for subsite in site for residue in subsite
    ] for site in product(*subfamily_combinations)]
    
    return potential_sites

    
def structure_family_site_to_vector(site):
    """Takes a structure family site object and turns it into a feature vector."""

    return []


def get_structure_half_families():
    """Get the families for which there are structure half-site models."""

    return ["H2", "C2"]


def structure_family_half_site_to_vector(half_site):
    """Takes a structure family half site object and turns it into a feature vector."""

    return []