import os
import json

# What families?
structure_families = [f[:-4] for f in os.listdir("data/csv/structure") if f.endswith("csv")]

# Get notebook template
with open("data/notebooks/structure-template.json") as f:
    template = f.read()


for family in structure_families:
    with open(f"data/notebooks/{family}-structure.ipynb", "w") as f:
        f.write(template.replace("X1", family))


