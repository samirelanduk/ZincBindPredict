#! /usr/bin/env python3

import sys
import os
import joblib
import json
import atomium
import time
from utilities import *

def update_status(job_id, status):
    with open(f"server/jobs/{job_id}/status.txt", "w") as f:
        f.write(status)


def update_results(job_id, results):
    with open(f"server/jobs/{job_id}/results.json", "w") as f:
        json.dump(results, f)

start = time.time()
job_id = sys.argv[1]
assembly = sys.argv[2]
filename = get_job_structure_file(job_id)
pdb = atomium.open(f"server/jobs/{job_id}/{filename}")
try:
    structure = pdb.generate_assembly(int(assembly)) if int(assembly) else pdb.model
    structure.optimise_distances()
except:
    update_status(job_id, f"Job failed: invalid assembly")
    sys.exit()


# what models are there?
models = [f[:-7] for f in os.listdir("predict/models/structure") if f.endswith("joblib")]
families = sorted(list(set(m.split("-")[0] for m in models)))


results = {"time": time.time() - start}
for model in models:
    update_status(job_id, f"Running {model}")
    family = model.split("-")[0]
    classifier = joblib.load(f"predict/models/structure/{model}.joblib")
    result = {
        "name": model, "sites": [], "rejected": [],
        "recall": round(classifier._test_recall, 3),
        "precision": round(classifier._test_precision, 3),
    }
    combos = [list(c) for c in model_to_residue_combos(structure, family)]
    inputs = [list(residues_to_sample(combo, "X").values())[1:] for combo in combos]
    if inputs:
        y = classifier.predict(inputs)
        prob = classifier.predict_proba(inputs)
        sites = [(combos[i], prob[i][1], inputs[i]) for i, o in enumerate(y) if o == 1]
        rejected = [(combos[i], prob[i][1], inputs[i]) for i, o in enumerate(y) if o != 1]
    else:
        sites = []
        rejected = []
    for site in sites:
        result["sites"].append({
            "probability": round(site[1], 3), "family": family, "model": model,
            "residues": [{"id": res.id, "name": res.name} for res in site[0]],
            "vector": site[2]
        })
    for site in rejected:
        result["rejected"].append({
            "probability": round(site[1], 3), "family": family, "model": model,
            "residues": [{"id": res.id, "name": res.name} for res in site[0]],
            "vector": site[2]
        })
    results[model] = result
    results["time"] = time.time() - start
    update_results(job_id, results)
