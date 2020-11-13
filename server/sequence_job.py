#! /usr/bin/env python3

"""This script finds zinc binding sites in a sequence."""

import sys
import json
import joblib
import traceback
from itertools import combinations, product
from collections import Counter
from utilities import *
from data.common import sequence_site_to_vector

from random import random
from time import sleep

def main():
    # Get arguments from JSON
    arguments = parse_arguments()
    job_id, sequence = arguments["job_id"], arguments["sequence"].upper()

    # Open JSON
    job = load_job(job_id)

    try:
        for family in get_sequence_families():
            if family in arguments.get("families", []) or not arguments.get("families"):
                # Update status
                save_job(job, status=f"Looking for {family} sites")

                # Find possible sites for this family
                possibles = list(sequence_to_family_inputs(sequence, family))

                # Convert possible sites to vectors
                dicts = [sequence_site_to_vector(possible) for possible in possibles]
                vectors = [list(d.values()) for d in dicts]
                if not vectors: continue
                max_length = max(len(v) for v in vectors)
                vectors = [v for v in vectors if len(v) == max_length]

                # Run vectors through models

                rf_model = joblib.load(f"predict/models/sequence/{family}-RF.joblib")
                predicted = rf_model.predict(vectors)
                probabilities = [p[1] for p in rf_model.predict_proba(vectors)]

                
                # Add sites to job object
                for site, positive, probability in zip(possibles, predicted, probabilities):
                    l = job["sites"] if positive and probability > 0.99 else job["rejected_sites"]
                    site = {
                        "probability": probability, "family": family, "residues": site
                    }
                    if positive and probability > 0.99: l.append(site)
                    l.sort(key=lambda s: -s["probability"])
                    
                # Save job
                save_job(job)
        
        # Finish job
        save_job(job, status="complete")

    except Exception as e:
        save_job(job, status="error")
        raise e


if __name__ == "__main__": main()