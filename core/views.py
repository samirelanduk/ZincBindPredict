import atomium
from django.http import JsonResponse

def predict(request):
    f = request.FILES["data"]
    pdb = atomium.utilities.parse_string(f.read().decode(), str(f))
    return JsonResponse({"structure": pdb.code})