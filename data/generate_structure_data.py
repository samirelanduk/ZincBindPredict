#! /usr/bin/env python3

import kirjava
import atomium
from tqdm import tqdm
import sys
sys.path.append("../zincbindpredict")
from utilities import *

API_URL = "https://api.zincbind.net/"

QUERY = """
query ZincSites($family: String) { pdbs(resolution__lt: 2) { edges { node {
    id assembly zincsites(family: $family) { edges { node { id residues {
        edges { node { atomiumId chainSignature atoms(name: "CA") {
            edges { node { x y z } }
        } } }
    } } } }
} } } }
"""

with open("data/families.dat") as f:
    FAMILIES = f.read().splitlines()

NEGATIVES = 100

def main():
    for family in FAMILIES:
        print(f"Getting data for {family} binding sites...")
        pdbs = fetch_data(API_URL, QUERY, variables={"family": family})
        pdbs = [pdb for pdb in pdbs if pdb["zincsites"]]
        update_data_file(family, "structure")
        site_count = sum([len(pdb["zincsites"]) for pdb in pdbs])
        print(f"Got {len(pdbs)} PDBs with {site_count} relevant binding sites")

        # Go through each PDB and get the relevant model
        for i, pdb in enumerate(tqdm(pdbs)):
            model = atomium.fetch(pdb["id"]).generate_assembly(pdb["assembly"])
            samples, seen_residue_ids = [], []

            # Get the positive cases
            for site in pdb["zincsites"]:
                # Get the actual atomium residues
                residues = [r for r in site["residues"] if r["chainSignature"]]
                try:
                    at_residues = get_residues_from_model(model, residues)
                except: continue

                # Make note of these residues for later
                seen_residue_ids.append(set([res.id for res in at_residues]))

                # Save sample
                sample = residues_to_sample(at_residues, site["id"])
                if sample: samples.append({**sample, **{"positive": 1}})
        
            # Get some negative cases
            negative_combos = model_to_residue_combos(
             model, family, ignore=seen_residue_ids
            )
            for combo in negative_combos:
                sample = residues_to_sample(combo, f"{pdb['id']}-N")
                if sample: samples.append({**sample, **{"positive": -1}})
            
            # Save samples to CSV
            update_data_file(family, "structure", samples)

        print()


if __name__ == "__main__":
    print()
    main()
    print()