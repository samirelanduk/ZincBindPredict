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
            if family in arguments["families"] or arguments["families"] == []:
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
                knn_model = joblib.load(f"predict/models/sequence/{family}-KNN.joblib")
                knn_predicted = knn_model.predict(vectors)
                knn_probabilities = [p[1] for p in knn_model.predict_proba(vectors)]

                rf_model = joblib.load(f"predict/models/sequence/{family}-RF.joblib")
                rf_predicted = rf_model.predict(vectors)
                rf_probabilities = [p[1] for p in rf_model.predict_proba(vectors)]

                svm_model = joblib.load(f"predict/models/sequence/{family}-SVM.joblib")
                svm_predicted = svm_model.predict(vectors)
                svm_probabilities = [p[1] for p in svm_model.predict_proba(vectors)]

                predicted = [Counter(votes).most_common()[0][0] for votes in zip(knn_predicted, rf_predicted, svm_predicted)]
                probabilities = [sum(votes) / 3 for votes in zip(knn_probabilities, rf_probabilities, svm_probabilities)]

                # Add sites to job object
                for site, positive, probability in zip(possibles, predicted, probabilities):
                    l = job["sites"] if positive else job["rejected_sites"]
                    site = {
                        "probability": probability, "family": family, "residues": site
                    }
                    l.append(site)
                    l.sort(key=lambda s: -s["probability"])
                    
                # Save job
                save_job(job)
        
        # Finish job
        save_job(job, status="complete")

    except Exception as e:
        save_job(job, status="error")
        raise e


if __name__ == "__main__": main()