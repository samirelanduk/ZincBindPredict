#! /usr/bin/env python3

import kirjava
import atomium
from tqdm import tqdm
import sys
sys.path.append("../zincbindpredict")
from core.utilities import residues_to_sample, model_to_residue_combos

API_URL = "http://localhost:8030/"

QUERY = """
query ZincSites($family: String) { pdbs(resolution__lt: 2) { edges { node {
    id assembly zincsites(family: $family) { edges { node { id residues {
        edges { node { atomiumId chainSignature atoms(name: "CA") {
            edges { node { x y z } }
        } } }
    } } } }
} } } }
"""

def main():
    client = kirjava.Client(API_URL)
    for family in ["H3", "C4", "C2H2", "C3H1", "D1H2", "E1H2"]:
        print(f"Getting data for {family} binding sites")
        print("Fetching data from server...")
        results = client.execute(QUERY, variables={"family": family})

        # Remove PDBs with no relevant binding sites
        pdbs = [edge["node"] for edge in results["data"]["pdbs"]["edges"]
         if edge["node"]["zincsites"]["edges"]]
        sites = sum([len(pdb["zincsites"]["edges"]) for pdb in pdbs])
        print(f"Found {len(pdbs)} PDBs with {sites} relevant binding sites")

        # Create list of data points that will be built up
        samples = []

        # Go through each PDB and get the relevant model
        for pdb in tqdm(pdbs[:2]):
            model = atomium.fetch(pdb["id"]).generate_assembly(pdb["assembly"])
            seen_ids = []

            # Get the positive cases
            for edge in pdb["zincsites"]["edges"]:
                site = edge["node"]
                site["residues"] = [e["node"] for e in site["residues"]["edges"]
                 if e["node"]["chainSignature"]]
                
                # Get atomium residues
                residues = [[r for r in model.residues(id=residue["atomiumId"]) 
                 if r.atom(name="CA").location == (
                  residue["atoms"]["edges"][0]["node"]["x"],
                  residue["atoms"]["edges"][0]["node"]["y"],
                  residue["atoms"]["edges"][0]["node"]["z"]
                 )
                ][0] for residue in site["residues"]]
                
                # Make note of these residues for later
                seen_ids.append(set([res.id for res in residues]))

                # Save sample
                sample = residues_to_sample(residues, site["id"])
                if sample:
                    sample["positive"] = 1
                    samples.append(sample)

            # Get the negative cases
            for combo in model_to_residue_combos(model, family):
                ids = set([res.id for res in combo])
                if ids not in seen_ids:
                    sample = residues_to_sample(combo, f"{pdb['id']}-N")
                    if sample:
                        sample["positive"] = -1
                        samples.append(sample)
                
        # Save samples to CSV
        csv_lines = [",".join(samples[0].keys()) + "\n"]
        for sample in samples:
            csv_lines.append(",".join([str(v) for v in sample.values()]) + "\n")
        with open(f"data/{family}.csv", "w") as f:
            f.writelines(csv_lines)
        print()


if __name__ == "__main__":
    print()
    main()
    print()