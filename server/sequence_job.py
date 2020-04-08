#! /usr/bin/env python3

"""This script finds zinc binding sites in a sequence."""

import sys
import json
import joblib
from itertools import combinations, product
from utilities import *

def sequence_to_family_inputs(sequence, family):
    """Takes a sequence and returns a list of family inputs."""

    sequence = sequence.lower()
    subfamilies = split_family(family)
    residue_combos = []
    for subfamily in subfamilies:
        residues = [i for i, char in enumerate(sequence)
            if char == subfamily[0].lower()]
        residue_combos.append(combinations(residues, int(subfamily[1:])))
    initial_combos = product(*residue_combos)
    return [c[0] for c in initial_combos]


def sequence_site_to_vector(site):
    """Takes a sequence site object and turns it into a feature vector."""
    return []

# Get arguments from JSON
if len(sys.argv) < 2:
    print("Please provide JSON arguments")
    sys.exit(1)
arguments = json.loads(sys.argv[1])
sequence = arguments["sequence"].upper()

# Open JSON
with open(get_job_location(arguments["job_id"])) as f: job = json.load(f)

families = ["H3", "C4", "C2H2"]

try:
    for family in families:
        # Update status
        job["status"] = f"Looking for {family} sites"
        save_job(job)

        # Find possible sites for this family
        possibles = list(sequence_to_family_inputs(sequence, family))

        # Convert possible sites to vectors
        vectors = [sequence_site_to_vector(possible) for possible in possibles]

        # Run vectors through models
        from random import random
        from time import sleep
        sleep(random() * 10)
        probabilities = [round(random(), 4) for _ in vectors]

        # Add sites to job object
        for site, probability in zip(possibles, probabilities):
            print(site, probability)
            l = job["sites"] if probability > 0.8 else job["rejected"]
            residues = [
                {"name": sequence[index], "identifier": index} for index in site
            ]
            l.append({"family": family, "probability": probability, "residues": residues})
        
        # Save job
        save_job(job)
    
    job["status"] = "complete"
    save_job(job)

except:
    job["status"] = "error"
    save_job(job)

"""# What is the sequence?
if len(sys.argv) < 3:
    print("Please provide sequence")
    sys.exit(1)
sequence = sys.argv[2].lower()

# Make sure there is a folder for this job
if not os.path.exists(f"server/jobs/{job_id}"):
    print("There's no job with that ID")
    sys.exit(1)
filename = get_job_structure_file(job_id)

# Create initial job JSON
job = {
 "job_id": job_id, "status": "initialising",
 "sequence": sequence, "families": {}
}
write_job(job)

# Determine what models are available
models = sorted([f for f in os.listdir("predict/models/sequence")
 if f.endswith("joblib")])
families = sorted(list(set([model.split("-")[0] for model in models])))
for family in families:
    job["families"][family] = {}
    for model in models:
        if model[:len(family)] == family:
            job["families"][family][
             model.split(".")[0].replace("Classifier", "")
            ] = {"status": "waiting"}
write_job(job)

# Go through families
for family in families:
    job["status"] = f"searching for {family} zinc binding sites"
    write_job(job)

    # Get potential sites
    combos = sequence_to_residue_combos(sequence, family)
    inputs = []
    for combo in combos:
        inputs.append(list(sequence_to_sample(combo, "X").values())[1:])
    
    # Get relevant models
    for model in models:
        if model[:len(family)] == family:

            # Find sites
            classifier = joblib.load(f"predict/models/sequence/{model}")
            if inputs:
                y = classifier.predict(inputs)
                prob = classifier.predict_proba(inputs)
                sites = [(combos[i], prob[i][1]) for i, o in enumerate(y) if o == 1]
            else:
                sites = []
            
            # Save sites to job
            model_name = model.split(".")[0].replace("Classifier", "")
            job["families"][family][model_name]["sites"] = [{
             "probability": round(prob, 3),
             "sequence": residues
            } for residues, prob in sites]
            job["families"][family][model_name]["status"] = "complete"
            job["families"][family][model_name]["metrics"] = {
             "validation_recall": round(classifier.validation_recall_, 3),
             "validation_precision": round(classifier.validation_precision_, 3),
             "test_recall": round(classifier.test_recall_, 3),
             "test_precision": round(classifier.test_precision_, 3)
            }
    
# Finish job
job["status"] = f"complete"
write_job(job)"""