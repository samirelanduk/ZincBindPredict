#! /usr/bin/env python3

import sys
sys.path.append("../zincbindpredict")
import kirjava
import atomium
from tqdm import tqdm
from common import structure_family_site_to_vector
from utilities import *
from common import location_to_vector, model_to_grid

API_URL = "https://api.zincbind.net/"

PDBS = """{
    pdbs(resolution__lt: 2) { edges { node { 
        id assembly zincsites { edges { node { id metals {
            edges { node { x y z } }
          } } } }
    } } }
}"""

# Get PDBs with zinc sites in
pdbs = fetch_data(API_URL, PDBS, {})
print(f"There are {len(pdbs)} PDBs")


# Go through each PDB
for pdb in tqdm(pdbs):

    # Build model
    try:
        model = atomium.fetch(pdb["id"]).generate_assembly(pdb["assembly"])
    except: continue
    model.optimise_distances()
    
    # Determine points to sample
    grid = model_to_grid(model)

    # Get positive locations
    positive_locations = set()
    for site in pdb["zincsites"]:
        for metal in site["metals"]:
            positive_locations.add((
                round(metal["x"]), round(metal["y"]), round(metal["z"])
            ))
    
    # Get negative locations
    negative_locations = set()
    while len(negative_locations) < len(positive_locations) * 10:
        point = random.choice(grid)
        if point not in positive_locations:
            negative_locations.add(point)

    # Get vectors for PDB
    positive_samples = [location_to_vector(point, model)
     for point in positive_locations]
    negative_samples = [location_to_vector(point, model)
     for point in negative_locations]
    
    # Save CSV
    save_csv(
        positive_samples, negative_samples, "full",
        os.path.join("data", "csv", "locations")
    )