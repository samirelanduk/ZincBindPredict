"""
For structure in structures:
	Work out single sites of these families
	Work out single sites of other families
	Work out half sites
	Pass to models
	List predicted sites
	Make three confusion matrices
"""

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
import atomium
from data.utilities import fetch_data
from data.common import structure_family_site_to_vector
from server.utilities import model_to_family_inputs
from collections import Counter
import joblib
from tqdm import tqdm

QUERY = """{ pdbs(first: 500) { edges { node { id assembly zincsites {
	edges { node { id family chainInteractions { count }
	residues(primary: true) { edges { node { 
		id atomiumId
	} } } } }
} } } } }"""

with open("data/families.dat") as f:
	FAMILIES = f.read().splitlines()

models = {f: joblib.load(f"predict/models/structure-families/{f}-RF.joblib") for f in FAMILIES}

pdbs = fetch_data("https://api.zincbind.net", QUERY, {})

true_positives = 0
false_negatives = 0
false_positives = 0
true_negatives = 0

family_true_positives = 0
family_false_negatives = 0
family_false_positives = 0
family_true_negatives = 0

for pdb in tqdm(pdbs):
	try:
		model = atomium.fetch(pdb["id"]).generate_assembly(pdb["assembly"])
	except: continue
	
	try:
		sites = []
		family_sites = []
		if len(pdb["zincsites"]) > 10: continue
		for site in pdb["zincsites"]:
			if site["chainInteractions"]["count"] == 1:
				s = frozenset([r["atomiumId"] for r in site["residues"]])
				sites.append(s)
				if site["family"] in FAMILIES:
					family_sites.append(s)
		
		predicted_sites = []
		for family in FAMILIES:
			rf_model = models[family]
			possibles = list(model_to_family_inputs(model, family))

			# Convert possible sites to vectors
			dicts = [structure_family_site_to_vector(possible) for possible in possibles]
			remove = []
			for i, d in enumerate(dicts):
				if d is None: remove.append(i)
			possibles = [p for i, p in enumerate(possibles) if i not in remove]
			dicts = [d for i, d in enumerate(dicts) if i not in remove]
			vectors = [list(d.values()) for d in dicts]
			if not vectors: continue
			max_length = max(len(v) for v in vectors)
			vectors = [v for v in vectors if len(v) == max_length]

			predicted = rf_model.predict(vectors)
			probabilities = [p[1] for p in rf_model.predict_proba(vectors)]
			
			for index, site in enumerate(possibles):
				if predicted[index] and probabilities[index] > 0.99:
					predicted_sites.append(frozenset([r.id for r in site]))
		
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
	except: break

		
	

	'''

	
	
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
			family_true_negatives += 1'''


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
