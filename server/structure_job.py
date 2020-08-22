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
        if arguments["use_families_models"] and (not arguments["find_half"]) and (family in arguments["families"] or arguments["families"] == []):
            # Update status
            save_job(job, status=f"Looking for {family} sites")

            # Find possible sites for this family
            possibles = list(model_to_family_inputs(model, family))

            # Convert possible sites to vectors
            dicts = [structure_family_site_to_vector(possible) for possible in possibles]
            vectors = [list(d.values()) for d in dicts]
            if not vectors: continue

            # Run vectors through models
            rf_model = joblib.load(f"predict/models/structure-families/{family}-RF.joblib")
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
    
    for half_family in get_structure_half_families():
        if arguments["use_families_models"] and arguments["find_half"]:
            # Update status
            save_job(job, status=f"Looking for {half_family} half-sites")

            # Find possible half sites for this family
            possibles = list(model_to_family_inputs(model, half_family))

            # Convert possible sites to vectors
            dicts = [structure_half_family_site_to_vector(possible) for possible in possibles]
            vectors = [list(d.values()) for d in dicts]
            if not vectors: continue

            # Run vectors through models
            knn_model = joblib.load(f"predict/models/structure-half-families/{half_family}-KNN.joblib")
            knn_predicted = knn_model.predict(vectors)
            knn_probabilities = [p[1] for p in knn_model.predict_proba(vectors)]

            rf_model = joblib.load(f"predict/models/structure-half-families/{half_family}-RF.joblib")
            rf_predicted = rf_model.predict(vectors)
            rf_probabilities = [p[1] for p in rf_model.predict_proba(vectors)]

            svm_model = joblib.load(f"predict/models/structure-half-families/{half_family}-SVM.joblib")
            svm_predicted = svm_model.predict(vectors)
            svm_probabilities = [p[1] for p in svm_model.predict_proba(vectors)]

            predicted = [Counter(votes).most_common()[0][0] for votes in zip(knn_predicted, rf_predicted, svm_predicted)]
            probabilities = [sum(votes) / 3 for votes in zip(knn_probabilities, rf_probabilities, svm_probabilities)]

            # Add sites to job object
            for site, positive, probability in zip(possibles, predicted, probabilities):
                l = job["sites"] if positive else job["rejected_sites"]
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
        locations = model_to_grid(model)

        # Convert possible locations to vectors
        vector_func = half_location_to_vector if arguments["find_half"] else location_to_vector
        dicts = [vector_func(location, model) for location in locations]
        vectors = [list(d.values()) for d in dicts]

        # Run vectors through models
        dataset = "half" if arguments["find_half"] else "full"
        knn_model = joblib.load(f"predict/models/locations/{dataset}-KNN.joblib")
        knn_predicted = knn_model.predict(vectors)
        knn_probabilities = [p[1] for p in knn_model.predict_proba(vectors)]

        rf_model = joblib.load(f"predict/models/locations/{dataset}-RF.joblib")
        rf_predicted = rf_model.predict(vectors)
        rf_probabilities = [p[1] for p in rf_model.predict_proba(vectors)]

        svm_model = joblib.load(f"predict/models/locations/{dataset}-SVM.joblib")
        svm_predicted = svm_model.predict(vectors)
        svm_probabilities = [p[1] for p in svm_model.predict_proba(vectors)]

        predicted = [Counter(votes).most_common()[0][0] for votes in zip(knn_predicted, rf_predicted, svm_predicted)]
        probabilities = [sum(votes) / 3 for votes in zip(knn_probabilities, rf_probabilities, svm_probabilities)]

        # Add sites to job object
        for site, positive, probability in zip(locations, predicted, probabilities):
            l = job["locations"] if positive else job["rejected_sites"]
            site = {
                "probability": probability, "location": site, "half": arguments["find_half"]
            }
            l.append(site)
            l.sort(key=lambda s: -s["probability"])

    # Finish job
    save_job(job, status="complete")

except Exception as e:
    save_job(job, status="error")
    raise e