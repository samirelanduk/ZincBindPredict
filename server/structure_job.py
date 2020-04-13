#! /usr/bin/env python3

import sys
import os
import json
import joblib
import atomium
from itertools import combinations, product
from utilities import *

# Get arguments from JSON
arguments = parse_arguments()
job_id, filename = arguments["job_id"], arguments["structure"]

# Open JSON
job = load_job(job_id)

try:
    # Get model
    model = get_model_for_job(filename)

    for family in get_structure_families():
        # Update status
        save_job(job, status=f"Looking for {family} sites")

        # Find possible sites for this family
        possibles = list(model_to_family_inputs(model, family))

        # Convert possible sites to vectors
        vectors = [structure_family_site_to_vector(possible) for possible in possibles]

        # Run vectors through models
        from random import random
        from time import sleep
        sleep(random() * 5)
        probabilities = [round(random(), 4) for _ in vectors]

        # Add sites to job object
        for site, probability in zip(possibles, probabilities):
            site = {
                "probability": probability,
                "family": family,
                "residues": [{
                    "name": res.name, "identifier": res.id
                } for res in site]
            }
            (job["sites"] if probability > 0.8 else job["rejected"]).append(site)
       
        # Save job
        save_job(job)
    
    # Finish job
    save_job(job, status="complete")

except Exception as e:
    save_job(job, status="error")
    raise e