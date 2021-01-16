#! /usr/bin/env python3

import sys
import kirjava
import atomium
from tqdm import tqdm
from common import structure_site_to_sample
from utilities import *

API_URL = "https://api.zincbind.net/"

FAMILY_SITES_QUERY = """query familySites($family: String) {
    zincsites(family: $family, pdb__resolution__lt: 2) { edges { node { 
        id pdb { id assembly } chainInteractions { edges { node { sequence } } }
        residues(primary: true) { edges { node {
            id chainSignature atomiumId atoms(name: "CA") {
                edges { node { x y z } }
            }
        } } }
    } } }
}"""

# Get options
families = parse_data_args(sys.argv, clustering=False)

# Produce dataset for all families
for family in families:
    print(f"Fetching {family} data...")

    # Get binding sites
    sites = fetch_data(API_URL, FAMILY_SITES_QUERY, {"family": family})
    
    # Create dataset
    positives, negatives, failures = [], [], []
    with tqdm(total=len(sites) * 2) as pbar:
        
        # Positives
        for site in sites:
            try:
                site = get_atomium_site(site)
                sample = structure_site_to_sample(site)
                if sample["ca_max"] <= 30: positives.append(sample)
            except Exception: failures.append(site["id"])
            pbar.update()
        
        # Negatives
        all_codes = get_all_pdb_codes()
        while len(negatives) != len(positives):
            code = random.choice(all_codes)
            try:
                site = get_random_atomium_site(code, family)
                if site:
                    sample = structure_site_to_sample(site)
                    if sample["ca_max"] <= 30:
                        negatives.append(sample)
                        pbar.update()
            except Exception: pass

    # Report any sites which could not be turned into a sample
    if len(failures): print("Could not process:", ", ".join(failures))
    
    # Create CSV file for family
    save_csv(
        positives, negatives, family, None,
        os.path.join("data", "csv", "structure")
    )