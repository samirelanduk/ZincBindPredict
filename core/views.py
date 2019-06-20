from itertools import combinations
import atomium
from django.http import JsonResponse

def split_family(family):
    subfamilies, subfamily = [], ""
    for char in family:
        if char.isalpha() and subfamily:
            subfamilies.append(subfamily)
            subfamily = ""
        subfamily += char
    subfamilies.append(subfamily)
    return subfamilies


def model_to_residue_combos(model, family):
    subfamilies = split_family(family)
    residue_combos = []
    for subfamily in subfamilies:
        residues = model.residues(code=subfamily[0])
        combos = list(combinations(residues, int(subfamily[1:])))
        residue_combos.append(combos)
    while len(residue_combos) > 1:
        residue_combos.insert(0, [(x,y) for x in residue_combos[0] for y in residue_combos[1]])
        residue_combos.pop(1)
        residue_combos.pop(1)
    return residue_combos[0]


def predict(request):
    f = request.FILES["data"]
    pdb = atomium.utilities.parse_string(f.read().decode(), str(f))
    model = pdb.model
    for family in ["H3", "C4", "C2H2M2"]:
        combos = model_to_residues(model, family)
    return JsonResponse({"structure": pdb.code})