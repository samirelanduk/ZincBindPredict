#! /usr/bin/env python3

import sys
import os
import joblib
import json
import atomium
from utilities import *

def update_status(job_id, status):
    with open(f"server/jobs/{job_id}/status.txt", "w") as f:
        f.write(status)


def update_results(job_id, results):
    with open(f"server/jobs/{job_id}/results.json", "w") as f:
        json.dump(results, f)


job_id = sys.argv[1]
filename = get_job_structure_file(job_id)
pdb = atomium.open(f"server/jobs/{job_id}/{filename}")
structure = pdb.model

update_status(job_id, "Identifying candidate residues")

# what models are there?
models = [f[:-7] for f in os.listdir("predict/models/structure") if f.endswith("joblib")]
families = sorted(list(set(m.split("-")[0] for m in models)))


results = {}
for model in models:
    family = model.split("-")[0]
    classifier = joblib.load(f"predict/models/structure/{model}.joblib")
    result = {
        "name": model, "sites": [],
        "validation_recall": round(classifier.validation_recall_, 3),
        "validation_precision": round(classifier.validation_precision_, 3),
        "test_recall": round(classifier.test_recall_, 3),
        "test_precision": round(classifier.test_precision_, 3)
    }
    combos = [list(c) for c in model_to_residue_combos(structure, family)]
    inputs = [list(residues_to_sample(combo, "X").values())[1:] for combo in combos]
    if inputs:
        y = classifier.predict(inputs)
        prob = classifier.predict_proba(inputs)
        sites = [(combos[i], prob[i][1]) for i, o in enumerate(y) if o == 1]
    else:
        sites = []
    for site in sites:
        result["sites"].append({
            "probability": round(site[1], 3), "family": family, "model": model,
            "residues": [{"id": res.id, "name": res.name} for res in site[0]],
        })
    results[model] = result
    update_results(job_id, results)

update_status(job_id, "complete")
    

'''# Make sure there is a folder for this job
if not os.path.exists(f"server/jobs/{job_id}"):
    print("There's no job with that ID")
    sys.exit(1)


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
models = sorted([f for f in os.listdir("predict/models/structure")
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
    combos = [list(c) for c in model_to_residue_combos(structure, family)]
    inputs = []
    for combo in combos:
        inputs.append(list(residues_to_sample(combo, "X").values())[1:])
    
    # Get relevant models
    for model in models:
        if model[:len(family)] == family:

            # Find sites
            classifier = joblib.load(f"predict/models/structure/{model}")
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
             "residues": [{"id": res.id, "name": res.name} for res in residues]
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
write_job(job)'''