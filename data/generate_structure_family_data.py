#! /usr/bin/env python3

import kirjava
import atomium
from tqdm import tqdm
from common import structure_family_site_to_vector
from utilities import *

API_URL = "https://api.zincbind.net/"

FAMILY_SITES_QUERY = """query familySites($family: String) {
    zincsites(family: $family, pdb__resolution__lt: 2) { edges { node { 
        id pdb { id assembly } residues(primary: true) { edges { node {
            id chainSignature atomiumId atoms(name: "CA") {
                edges { node { x y z } }
            }
        } } }
    } } }
}"""

# What families should be used?
with open("data/families.dat") as f: families = f.read().splitlines()

for family in families:
    # Download all binding sites for this family
    print(f"Fetching {family} data...")
    family_sites = fetch_data(API_URL, FAMILY_SITES_QUERY, {"family": family})
    
    # Keep track of some things
    positive_samples = []
    negative_samples = []
    current_pdb = None
    current_model = None
    positive_sites_from_model = []

    # Go through each site
    for site in tqdm(family_sites):
        # Is this a new PDB?
        if site["pdb"]["id"] != current_pdb:
            # Make negative samples before moving on
            if current_pdb is not None:
                negative_samples += create_negative_samples_for_model(
                    current_model, family, positive_sites_from_model
                )

            # Get model for this PDB
            pdb = atomium.fetch(site["pdb"]["id"])
            try:
                current_model = pdb.generate_assembly(site["pdb"]["assembly"])
            except ValueError: continue
            current_pdb = site["pdb"]["id"]
            current_model.optimise_distances()
            positive_sites_from_model = []

        # Get residues within model
        site_residues = [r for r in site["residues"] if r["chainSignature"]]
        try:
            at_residues = get_residues_from_model(current_model, site_residues)
        except: continue

        # Add vector to positives
        positive_samples.append(structure_family_site_to_vector(at_residues))

        # Add residues to the list of positive residue combinations for this pdb
        positive_sites_from_model.append(set(at_residues))

    # Make negative samples for final PDB
    negative_samples += create_negative_samples_for_model(
        current_model, family, positive_sites_from_model
    )
    
    # Create CSV file for family
    save_csv(
        positive_samples, negative_samples, family,
        os.path.join("data", "csv", "structure-families")
    )


    
