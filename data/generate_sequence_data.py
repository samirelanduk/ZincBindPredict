#! /usr/bin/env python3

import sys
sys.path.append("../zincbindpredict")
import kirjava
import random
from tqdm import tqdm
from common import sequence_site_to_vector
from utilities import *

API_URL = "https://api.zincbind.net/"

ALL_CHAINS_QUERY = """{chainClusters { edges { node { id chains(first: 1) {
    edges { node { sequence } }
} } } } }"""

FAMILY_CHAINS_QUERY = """query familySites($family: String) {
    zincsites(family: $family) { edges { node { 
        id chainInteractions { edges { node { sequence } } }
    } } }
}"""

# Download all unique chains - these will be used for generating negatives
clusters = fetch_data(API_URL, ALL_CHAINS_QUERY, {})
unique_sequences = [cluster["chains"][0]["sequence"] for cluster in clusters]
print(f"Using {len(unique_sequences)} unique sequences for negative samples")

# What families should be used?
with open("data/families.dat") as f: families = f.read().splitlines()
for arg in sys.argv:
    if arg.startswith("--limit="):
        families = [f for f in families if f in arg[8:].split(",")]
    if arg.startswith("--exclude="):
        families = [f for f in families if f not in arg[10:].split(",")]

for family in families:
    # Download all binding sites for this family
    print(f"Fetching {family} data...")
    family_sites = fetch_data(API_URL, FAMILY_CHAINS_QUERY, {"family": family})
    res_count = sum([int(c) for c in family if c.isdigit()])
    
    # Get one sequence for each site, and only use sites on one chain
    family_sequences = [
        site["chainInteractions"][0]["sequence"] for site in family_sites
         if len(site["chainInteractions"]) == 1
    ]

    # How many sequence sites are there and how many negatives should there be?
    positive_count = len(family_sequences)
    negative_count = positive_count * 10
    with tqdm(total=positive_count + negative_count) as pbar:

        # Get positive samples for them
        positive_samples = []
        for sequence in family_sequences:
            if len([c for c in sequence if c.isupper()]) == res_count:
                positive_samples.append(sequence_site_to_vector(sequence))
                pbar.update()

        # Get negative samples for this family - 10 per positive sample
        negative_samples = []
        while len(negative_samples) < negative_count:
            # Pick a random unique_chain
            unique_chain = random.sample(unique_sequences, 1)[0]

            # Pick a random site within it
            site = random_sequence_family_input(unique_chain, family)
            if not site: continue

            # If it's not actually a positive sample, add it
            if site not in family_sequences:
                negative_samples.append(sequence_site_to_vector(site))
                pbar.update()
    
    # Create CSV file for family
    save_csv(
        positive_samples, negative_samples, family,
        os.path.join("data", "csv", "sequence")
    )