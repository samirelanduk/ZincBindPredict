#! /usr/bin/env python3

"""This script finds zinc binding sites in a sequence."""

import sys
import json
import joblib
import traceback
from itertools import combinations, product
from utilities import *

# Get arguments from JSON
arguments = parse_arguments()
job_id, sequence = arguments["job_id"], arguments["sequence"].upper()

# Open JSON
job = load_job(job_id)

try:
    for family in get_sequence_families():
        print(family)
        # Update status
        save_job(job, status=f"Looking for {family} sites")

        # Find possible sites for this family
        possibles = list(sequence_to_family_inputs(sequence, family))

        # Convert possible sites to vectors
        vectors = [sequence_site_to_vector(possible) for possible in possibles]

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
                "residues": site
            }
            (job["sites"] if probability > 0.8 else job["rejected"]).append(site)
            
        # Save job
        save_job(job)
    
    # Finish job
    save_job(job, status="complete")

except Exception as e:
    save_job(job, status="error")
    raise e