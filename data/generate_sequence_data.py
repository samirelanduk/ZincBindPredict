#! /usr/bin/env python3

import sys
sys.path.append("../zincbindpredict")
import kirjava
import random
from tqdm import tqdm
import subprocess
import re
from common import sequence_site_to_sample
from utilities import *

API_URL = "https://api.zincbind.net/"

FAMILY_CHAINS_QUERY = """query familySites($family: String) {
    zincsites(family: $family) { edges { node { 
        id chainInteractions { edges { node { sequence } } }
    } } }
}"""

# Get options
families, clustering = parse_data_args(sys.argv)

# Produce dataset for all families
for family in families:
    print(f"Fetching {family} data...")

    # Get sequence per binding site
    sites = fetch_data(API_URL, FAMILY_CHAINS_QUERY, {"family": family})
    sequences = [site["chainInteractions"][0]["sequence"] for site in sites if\
        len(site["chainInteractions"]) == 1]
    sequences = [s for s in sequences if sequence_contains_family(s, family)]

    # Cluster them if necessary
    if clustering: sequences = cluster_sequences(sequences, clustering)

    # Create dataset
    positives, negatives = [], []
    with tqdm(total=len(sequences) * 2) as pbar:

        # Positives
        for sequence in sequences:
            positives.append(sequence_site_to_sample(sequence))
            pbar.update()

        # Negatives
        all_sequences = get_all_uniprot_sequences()
        while len(negatives) != len(positives):
            sequence = random.choice(all_sequences).lower()
            site = get_random_sequence_site(sequence, family)
            if site:
                negatives.append(sequence_site_to_sample(site))
                pbar.update()
    
    # Create CSV file for family
    save_csv(
        positives, negatives, family, clustering,
        os.path.join("data", "csv", "sequence")
    )