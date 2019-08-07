import atomium
from django.http import JsonResponse
from data.utilities import model_to_residue_combos, residues_to_sample

def predict(request):
    f = request.FILES["data"]
    pdb = atomium.utilities.parse_string(f.read().decode(), str(f))
    model = pdb.model
    for family in ["H3", "C4", "C2H2M2"]:
        combos = model_to_residue_combos(model, family)
        for combo in combos:
            residues_to_sample(combo)
    return JsonResponse({"structure": pdb.code, "sites-found": 0})