#! /usr/bin/env python3

import sys
import os
import json
import joblib
import atomium
from itertools import combinations, product
from collections import Counter
from utilities import *
from data.common import *

# Get arguments from JSON
arguments = parse_arguments()
job_id, filename = arguments["job_id"], arguments["structure"]

# Open JSON
job = load_job(job_id)

try:
    # Get model
    model = get_model_for_job(filename)

    for family in get_structure_families():
        family = family.split("_")[0]
        # Update status
        save_job(job, status=f"Looking for {family} sites")

        # Find possible sites for this family
        possibles = list(model_to_family_inputs(model, family))

        # Convert possible sites to vectors
        dicts = [structure_site_to_sample(possible) for possible in possibles]
        remove = []
        for i, d in enumerate(dicts):
            if d is None: remove.append(i)
        possibles = [p for i, p in enumerate(possibles) if i not in remove]
        dicts = [d for i, d in enumerate(dicts) if i not in remove]
        vectors = [list(d.values()) for d in dicts]
        if not vectors: continue

        # Run vectors through models
        rf_model = joblib.load(f"predict/models/structure/{family}_100.joblib")
        predicted = rf_model.predict(vectors)
        probabilities = [p[1] for p in rf_model.predict_proba(vectors)]

        # Add sites to job object
        for site, positive, probability in zip(possibles, predicted, probabilities):
            l = job["sites"] if positive and probability > 0.95 else job["rejected_sites"]
            site = {
                "probability": probability, "family": family, "half": False,
                "residues": [{"name": res.name, "identifier": res.id} for res in site]
            }
            l.append(site)
            l.sort(key=lambda s: -s["probability"])
                
        # Save job
        save_job(job)
    

    save_job(job, status="complete")

except Exception as e:
    save_job(job, status="error")
    raise e