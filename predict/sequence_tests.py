"""
For sequence in sequences:
	Work out single sites of these families
	Work out single sites of other families
	Work out half sites
	Pass to models
	List predicted sites
	Make three confusion matrices
"""

import sys
sys.path.append("../zincbindpredict")
from data.utilities import fetch_data
from data.common import sequence_site_to_vector
from server.utilities import sequence_to_family_inputs
from collections import Counter
import joblib
from tqdm import tqdm

QUERY = """{ chains(first: 500) { edges { node { id sequence chainInteractions {
	edges { node { id sequence site { family } } }
} } } } }"""

with open("data/families.dat") as f:
	FAMILIES = f.read().splitlines()

models = {f: joblib.load(f"predict/models/sequence/{f}-RF.joblib") for f in FAMILIES}

chains = fetch_data("http://localhost:7000", QUERY, {})

true_positives = 0
false_negatives = 0
false_positives = 0
true_negatives = 0

family_true_positives = 0
family_false_negatives = 0
family_false_positives = 0
family_true_negatives = 0

for chain in tqdm(chains):
	sequence = chain["sequence"]
	summary = sequence[:5] + "..." + sequence[-5:]
	
	interactions = chain["chainInteractions"]

	sites = []
	family_sites = []
	for interaction in interactions:
		counter = Counter([c for c in interaction["sequence"] if c.isupper()])
		counter = sorted(counter.most_common(), key=lambda x: x[0])
		family = "".join(["".join([c, str(i)]) for c, i in counter])
		if family == interaction["site"]["family"]:
			sites.append(interaction["sequence"])
			if family in FAMILIES:
				family_sites.append(interaction["sequence"])

	predicted_sites = []
	for family in FAMILIES:
		model = models[family]
		possibles = list(sequence_to_family_inputs(sequence, family))

		# Convert possible sites to vectors
		dicts = [sequence_site_to_vector(possible) for possible in possibles]
		vectors = [list(d.values()) for d in dicts]
		if not vectors: continue
		max_length = max(len(v) for v in vectors)
		vectors = [v for v in vectors if len(v) == max_length]

		predicted = model.predict(vectors)
		probabilities = [p[1] for p in model.predict_proba(vectors)]
		
		for index, site in enumerate(possibles):
			if predicted[index] and probabilities[index] > 0.99:
				predicted_sites.append(site)
	
	for site in sites:
		if site in predicted_sites:
			true_positives += 1
		else:
			false_negatives += 1
	for site in predicted_sites:
		if site not in sites:
			false_positives += 1
	for site in possibles:
		if site not in predicted_sites and site not in sites:
			true_negatives += 1
	
	for site in family_sites:
		if site in predicted_sites:
			family_true_positives += 1
		else:
			family_false_negatives += 1
	for site in predicted_sites:
		if site not in family_sites:
			family_false_positives += 1
	for site in possibles:
		if site not in predicted_sites and site not in family_sites:
			family_true_negatives += 1


print("TP", true_positives)	
print("FP", false_positives)	
print("FN", false_negatives)	
print("TN", true_negatives)
print("Recall", true_positives / (true_positives + false_negatives))
print("Precision", true_positives / (true_positives + false_positives))
print()
print("family_TP", family_true_positives)	
print("family_FP", family_false_positives)	
print("family_FN", family_false_negatives)	
print("family_TN", family_true_negatives)
print("family_Recall", family_true_positives / (family_true_positives + family_false_negatives))
print("family_Precision", family_true_positives / (family_true_positives + family_false_positives))
