#! /usr/bin/env python3

import sys
import kirjava
import atomium
from atomium.data import CODES
from tqdm import tqdm
from common import structure_half_family_site_to_vector
from utilities import *

API_URL = "https://api.zincbind.net/"

ALL_SITES = """{
    zincsites(pdb__resolution__lt: 2) { edges { node {
        id chainInteractions { edges { node { id sequence } } }
    } } }
}"""

SITE = """query site($id: String!) { zincsite(id: $id) {
    pdb { id assembly } residues(primary: true) { edges { node {
        chainIdentifier name atomiumId atoms(name: "CA") {
            edges { node { x y z } }
        }
    }}}
} }"""

# What families should be used?
with open("data/half-families.dat") as f: families = f.read().splitlines()
for arg in sys.argv:
    if arg.startswith("--limit="):
        families = [f for f in families if f in arg[8:].split(",")]
    if arg.startswith("--exclude="):
        families = [f for f in families if f not in arg[10:].split(",")]
families = {family: [] for family in families}

# Get all half site IDs
print("Collating half-sites...")
sites = fetch_data(API_URL, ALL_SITES, {})
split_sites = [site for site in sites if len(site["chainInteractions"]) > 1]
for site in split_sites:
    for interaction in site["chainInteractions"]:
        family = chars_to_family([
            char for char in interaction["sequence"] if char.isupper()
        ])
        if family in families: families[family].append(site["id"])

# Go through each family
for family, sites in families.items():
    print(f"{family} half-sites...")

    # Keep track of some things
    positive_samples = []
    negative_samples = []
    current_pdb = None
    current_model = None
    positive_sites_from_model = []


    # Go through each site for that family
    for site_id in tqdm(sorted(set(sites))):

        # Download site
        site = fetch_data(API_URL, SITE, {"id": site_id})

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

        # Check each chain in site to see if it contains a relevant half site
        chains = set(residue["chainIdentifier"] for residue in site["residues"])
        for chain in chains:
            residues = [r for r in site["residues"] if r["chainIdentifier"] == chain]
            residue_family = chars_to_family([
                CODES.get(r["name"], "X") for r in residues
            ])
            if residue_family == family:

                # Got a half-site
                try:
                    at_residues = get_residues_from_model(current_model, residues)
                except: continue

                # Add vector to positives
                positive_samples.append(structure_half_family_site_to_vector(at_residues))

                # Add residues to the list of positive residue combinations for this pdb
                positive_sites_from_model.append(set(at_residues))
        
    # Make negative samples for final PDB
    negative_samples += create_negative_samples_for_model(
        current_model, family, positive_sites_from_model
    )

    # Create CSV file for family
    save_csv(
        positive_samples, negative_samples, family,
        os.path.join("data", "csv", "structure-half-families")
    )