import os
import json

# What families?
families = [f[:-4] for f in os.listdir("data/csv") if f.endswith("csv")]

# Get notebook template
with open("data/notebooks/template.json") as f:
    template = f.read()


for family in families:
    with open(f"data/notebooks/{family}.ipynb", "w") as f:
        f.write(template.replace("X1", family))


