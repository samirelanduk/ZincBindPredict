#! /usr/bin/env python3

import os
import json

for category in ["structure", "sequence"]:
    # What families?
    families = [f[:-4] for f in os.listdir(f"data/csv/{category}") if f.endswith("csv")]

    # Get notebook template
    with open(f"data/notebooks/{category}-template.json") as f:
        template = f.read()

    # Make notebooks
    for family in families:
        with open(f"data/notebooks/{family}-{category}.ipynb", "w") as f:
            f.write(template.replace("X1", family))