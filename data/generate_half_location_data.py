#! /usr/bin/env python3

import sys
sys.path.append("../zincbindpredict")
import kirjava
import atomium
from tqdm import tqdm
from common import structure_family_site_to_vector
from utilities import *
from common import half_location_to_vector, model_to_grid

API_URL = "https://api.zincbind.net/"

PDBS = """{
    pdbs(resolution__lt: 2) { edges { node { 
        id assembly zincsites { edges { node {
            id metals { edges { node { x y z } } } chainInteractions {
                edges { node { id chain { atomiumId } } }
            }
        } } }
    } } }
}"""

# Get PDBs with zinc sites in
pdbs = fetch_data(API_URL, PDBS, {})
pdbs = [pdb for pdb in pdbs if len([
    site for site in pdb["zincsites"] if
    len(site["chainInteractions"]) > 1 and len(site["metals"]) == 1
])]
print(f"There are {len(pdbs)} PDBs with half-sites")

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
        if len(site["chainInteractions"]) > 1 and len(site["metals"]) == 1:
            metal = site["metals"][0]
            for interaction in site["chainInteractions"]:
                positive_locations.add((
                    (round(metal["x"]), round(metal["y"]), round(metal["z"])),
                    interaction["chain"]["atomiumId"]
                ))
    
    # Get negative locations
    negative_locations = set()
    while len(negative_locations) < len(positive_locations) * 10:
        point = random.choice(grid)
        if point not in positive_locations:
            list(model.atoms_in_sphere(point, 4))[0].chain.id
            negative_locations.add((
                point, list(model.atoms_in_sphere(point, 4))[0].chain.id
            ))

    # Get vectors for PDB
    positive_samples = [half_location_to_vector(point, chain, model)
     for point, chain in positive_locations]
    negative_samples = [half_location_to_vector(point, chain, model)
     for point, chain in negative_locations]
    
    # Save CSV
    save_csv(
        positive_samples, negative_samples, "half",
        os.path.join("data", "csv", "locations")
    )