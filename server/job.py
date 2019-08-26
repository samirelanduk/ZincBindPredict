#! /usr/bin/env python3

import sys
import atomium
import joblib
from utilities import *

# What is the job ID
if len(sys.argv) < 2:
    print("Please provide Job ID")
    sys.exit(1)
job_id = sys.argv[1]

# Make sure there is a folder for this job
if not os.path.exists(f"server/jobs/{job_id}"):
    print("There's no job with that ID")
    sys.exit(1)
filename = get_job_structure_file(job_id)

# Create initial job JSON
job = {
 "job_id": job_id, "status": "initialising",
 "structure": "unknown", "families": {}
}
write_job(job)

# Load structure
pdb = atomium.open(f"server/jobs/{job_id}/{filename}")
structure = pdb.model
job["structure"] = pdb.code or "unknown"
write_job(job)

# Determine what models are available
models = sorted([f for f in os.listdir("predict/models")
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
    combos = model_to_residue_combos(structure, family)
    inputs = []
    for combo in combos:
        inputs.append(list(residues_to_sample(combo, "X").values())[1:])
    
    # Get relevant models
    for model in models:
        if model[:len(family)] == family:

            # Find sites
            classifier = joblib.load(f"predict/models/{model}")
            if inputs:
                y = classifier.predict(inputs)
                sites = [combos[i] for i, o in enumerate(y) if o == 1]
            else:
                sites = []
            
            # Save sites to job
            model_name = model.split(".")[0].replace("Classifier", "")
            job["families"][family][model_name]["sites"] = [[
             {"id": res.id, "name": res.name} for res in site
            ] for site in sites]
            job["families"][family][model_name]["status"] = "complete"
            job["families"][family][model_name]["metrics"] = {
             "validation_recall": classifier.validation_recall_,
             "validation_precision": classifier.validation_precision_,
             "test_recall": classifier.test_recall_,
             "test_precision": classifier.test_precision_
            }
    
# Finish job
job["status"] = f"complete"
write_job(job)