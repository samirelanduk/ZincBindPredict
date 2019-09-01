#! /usr/bin/env python3

import kirjava
import atomium
from tqdm import tqdm
import sys
sys.path.append("../zincbindpredict")
from utilities import *

API_URL = "https://api.zincbind.net/"

QUERY = """
{ chainClusters { edges { node {
    id chains(first: 1) { edges { node { id sequence chainInteractions {
        edges { node { id sequence site { family } } }
    }}}}
}}}}
"""

with open("data/families.dat") as f:
    FAMILIES = f.read().splitlines()

NEGATIVES = 100

def main():
    for family in FAMILIES:
        # Get data
        print(f"Getting data for {family} binding sites...")
        chain_clusters = fetch_data(API_URL, QUERY, variables={})
        chains = [c["chains"][0] for c in chain_clusters]
        for chain in chains:
            chain["chainInteractions"] = [i for i in chain["chainInteractions"]
             if i["site"]["family"] == family and
              residue_sequence_count(i["sequence"]) == family_count(family)]
        chains = [chain for chain in chains if chain["chainInteractions"]]
        print(f"There are {len(chains)} usable, unique {family} chains")
        update_data_file(family, "sequence")
        
        # Go through each chain
        for chain in tqdm(chains):
            # Get the positive cases
            samples = []
            for i, interaction in enumerate(chain["chainInteractions"], start=1):
                sample = sequence_to_sample(
                 interaction["sequence"], f"{chain['id']}-{i}"
                )
                if sample: samples.append({**sample, **{"positive": 1}})
            
            # Get the negative cases
            exclude = [tuple([
             n for n, char in enumerate(i["sequence"]) if char.isupper()
            ]) for i in chain["chainInteractions"]]
            negative_combos = sequence_to_residue_combos(
             chain["sequence"], family, limit=len(samples) * NEGATIVES, ignore=exclude
            )
            for combo in negative_combos:
                sample = sequence_to_sample(combo, f"{chain['id']}-N")
                if sample: samples.append({**sample, **{"positive": -1}})

            # Save samples to CSV
            update_data_file(family, "sequence", samples)

if __name__ == "__main__":
    print()
    main()
    print()