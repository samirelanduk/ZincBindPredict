#! /usr/bin/env python3

import kirjava
import atomium
from tqdm import tqdm
import sys
sys.path.append("../zincbindpredict")
from core.utilities import sequence_to_sample, sequence_to_residue_combos

API_URL = "https://api.zincbind.net/"

QUERY = """
query ZincSites($family: String) {
  zincsites(family: $family) {
    edges {
      node {
        id
        chainInteractions {
          edges {
            node {
              sequence
              chain {
                id
                sequence
              }
            }
          }
        }
      }
    }
  }
}
"""

def main():
    client = kirjava.Client(API_URL)
    for family in ["H3", "C4", "C2H2", "C3H1", "D1H2", "E1H2"]:
        print(f"Getting data for {family} binding sites")
        print("Fetching data from server...")
        results = client.execute(QUERY, variables={"family": family})
        with open(f"data/{family}-seq.csv", "w") as f: f.write("")

        # Organise data
        chains = {}
        for site_edge in results["data"]["zincsites"]["edges"]:
            site = site_edge["node"]
            if len(site["chainInteractions"]["edges"]) > 1: continue
            interaction = site["chainInteractions"]["edges"][0]["node"]
            chain = interaction["chain"]
            if chain["id"] not in chains:
                chains[chain["id"]] = {
                 "sequence": chain["sequence"], "sites": {}
                }
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
            for indices in sequence_to_residue_combos(chain["sequence"], family, limit=len(samples) * 100):
                pass



            # Save samples to CSV
            csv_lines = [",".join(samples[0].keys()) + "\n"] if i == 0 else []
            for sample in samples:
                csv_lines.append(",".join([str(v) for v in sample.values()]) + "\n")
            with open(f"data/{family}-seq.csv", "a") as f:
                f.writelines(csv_lines)

        '''

        # Go through each PDB
        for i, pdb in enumerate(tqdm(pdbs)):
            

            # Get the positive cases
            for edge in pdb["zincsites"]["edges"]:
                site = edge["node"]
                site["residues"] = [e["node"] for e in site["residues"]["edges"]
                 if e["node"]["chainSignature"]]
                
                # Get atomium residues
                try:
                    residues = [[r for r in model.residues(id=residue["atomiumId"]) 
                    if r.atom(name="CA").location == (
                    residue["atoms"]["edges"][0]["node"]["x"],
                    residue["atoms"]["edges"][0]["node"]["y"],
                    residue["atoms"]["edges"][0]["node"]["z"]
                    )
                    ][0] for residue in site["residues"]]
                except: continue
                
                # Make note of these residues for later
                seen_ids.append(set([res.id for res in residues]))

                # Save sample
                sample = residues_to_sample(residues, site["id"])
                if sample:
                    sample["positive"] = 1
                    samples.append(sample)

            # Get the negative cases
            for combo in model_to_residue_combos(model, family, limit=len(samples) * 100):
                ids = set([res.id for res in combo])
                if ids not in seen_ids:
                    sample = residues_to_sample(combo, f"{pdb['id']}-N")
                    if sample:
                        sample["positive"] = -1
                        samples.append(sample)
                
            '''
        print()


if __name__ == "__main__":
    print()
    main()
    print()