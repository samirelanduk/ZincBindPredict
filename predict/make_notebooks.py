#! /usr/bin/env python3

import os
import json

for category in ["structure", "sequence"]:
    # Get models
    models = [m for m in os.listdir(f"predict/models/{category}")
     if m.endswith("joblib")]
    
    # Get notebook template
    with open(f"predict/notebooks/template.json") as f:
        template = f.read()
    
    # Make notebooks
    for model in models:
        name = model.replace("Classifier", "")[:-7]
        with open(f"predict/notebooks/{category}-{name}.ipynb", "w") as f:
            f.write(template.replace("_MODEL_", name).replace("_CATEGORY_", category).replace("_PATH_", model))

 


    '''# What families?
    families = [f[:-4] for f in os.listdir(f"data/csv/{category}") if f.endswith("csv")]

    # Get notebook template
    with open(f"data/notebooks/{category}-template.json") as f:
        template = f.read()

    # Make notebooks
    for family in families:
        with open(f"data/notebooks/{family}-{category}.ipynb", "w") as f:
            f.write(template.replace("X1", family))'''