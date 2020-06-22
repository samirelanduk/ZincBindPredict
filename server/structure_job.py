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
        if arguments["use_families_models"] and (not arguments["find_half"]) and (family in arguments["families"] or arguments["families"] == []):
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
                l = job["sites"] if probability > 0.8 else job["rejected_sites"]
                site = {
                    "probability": probability, "family": family, "half": False,
                    "residues": [{"name": res.name, "identifier": res.id} for res in site]
                }
                l.append(site)
                l.sort(key=lambda s: -s["probability"])
        
            # Save job
            save_job(job)
    
    for half_family in get_structure_half_families():
        if arguments["use_families_models"] and arguments["find_half"]:
            # Update status
            save_job(job, status=f"Looking for {half_family} half-sites")

            # Find possible half sites for this family
            possibles = list(model_to_family_inputs(model, half_family))

            # Convert possible sites to vectors
            vectors = [structure_family_half_site_to_vector(possible) for possible in possibles]

            # Run vectors through models
            from random import random
            from time import sleep
            sleep(random() * 5)
            probabilities = [round(random(), 4) for _ in vectors]

            # Add sites to job object
            for site, probability in zip(possibles, probabilities):
                l = job["sites"] if probability > 0.8 else job["rejected_sites"]
                site = {
                    "probability": probability, "family": half_family, "half": True,
                    "residues": [{"name": res.name, "identifier": res.id} for res in site]
                }
                l.append(site)
                l.sort(key=lambda s: -s["probability"])
            
            # Save job
            save_job(job)
    
    if arguments["use_location_models"]:
        # Start looking through locations in model
        save_job(job, status=f"Looking at locations in structure")

        # Get possible centres of zinc binding
        possibles = list(get_model_locations(model))

        # Convert possible sites to vectors
        vectors = [structure_location_to_vector(possible, model) for possible in possibles]
        half_vectors = [structure_location_to_half_vector(possible, model) for possible in possibles]

        # Run vectors through model
        from random import random
        from time import sleep
        sleep(random() * 5)
        probabilities = [round(random(), 4) for _ in vectors]
        sleep(random() * 5)
        half_probabilities = [round(random(), 4) for _ in half_vectors]

        # Add sites to job object
        for site, probability, half_probability in zip(possibles, probabilities, half_probabilities):
            full = probability > 0.8
            half = not full and half_probability > 0.8
            l = job["locations"] if probability > 0.8 else job["rejected_locations"]
            site = {
                "probability": probability, "location": site, "half": half
            }
            l.append(site)
            l.sort(key=lambda s: -s["probability"])

    # Finish job
    save_job(job, status="complete")

except Exception as e:
    save_job(job, status="error")
    raise e