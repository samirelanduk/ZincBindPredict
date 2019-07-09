#! /usr/bin/env python3

import kirjava
import atomium
from tqdm import tqdm
import sys
sys.path.append("../zincbindpredict")
from core.utilities import sequence_to_sample, sequence_to_residue_combos, split_family

API_URL = "https://api.zincbind.net/"

QUERY = """
query ZincSites($family: String) { zincsites(family: $family) { edges { node {
    id chainInteractions { edges { node { sequence chain { id sequence } } } }
} } } }
"""

def main():
    client = kirjava.Client(API_URL)
    for family in ["H3", "C4", "C2H2", "C3H1", "D1H2", "E1H2"]:
        print(f"Getting data for {family} binding sites")
        print("Fetching data from server...")
        results = client.execute(QUERY, variables={"family": family})
        with open(f"data/{family}-seq.csv", "w") as f: f.write("")

        # Remove sites with more than one chain
        sites = [edge["node"] for edge in results["data"]["zincsites"]["edges"] 
         if len(edge["node"]["chainInteractions"]["edges"]) == 1]
        
        # Remove sites with insufficient upper case residues
        sites = [site for site in sites if len([c for c in site["chainInteractions"]["edges"][0]["node"]["sequence"] if c.isupper()]) == sum([int(s[1:]) for s in split_family(family)])  ]
        
        # Organise into chains and chain interactions
        chains = {}
        for site in sites:
            interaction = site["chainInteractions"]["edges"][0]["node"]
            chain = interaction["chain"]
            if chain["id"] not in chains:
                chains[chain["id"]] = {"sequence": chain["sequence"], "sites": {}} 
            chains[chain["id"]]["sites"][site["id"]] = interaction["sequence"]
        sites = sum([len(chain["sites"]) for chain in chains.values()])
        print(f"Found {len(chains)} chains with {sites} relevant binding sites")
        
        # Go through each chain
        for i, chain_pair in enumerate(tqdm(chains.items())):
            chain_id, chain = chain_pair
            samples = []

            # Get the positive cases
            for site_id, site in chain["sites"].items():
                sample = sequence_to_sample(site, site_id)
                sample["positive"] = 1
                samples.append(sample)
            
            # Get the negative cases
            for sequence in sequence_to_residue_combos(chain["sequence"], family, limit=len(samples) * 100):
                sample = sequence_to_sample(sequence, f"{chain_id}-N")
                sample["positive"] = -1
                samples.append(sample)

            # Save samples to CSV
            csv_lines = [",".join(samples[0].keys()) + "\n"] if i == 0 else []
            for sample in samples:
                csv_lines.append(",".join([str(v) for v in sample.values()]) + "\n")
            with open(f"data/{family}-seq.csv", "a") as f:
                f.writelines(csv_lines)

        print()


if __name__ == "__main__":
    print()
    main()
    print()