import atomium
from django.http import JsonResponse

PATHS = {
 "/predict/structure?code=XXXX": "Find zinc binding sites in PDB structure"
}

def root(request):
    return JsonResponse({
     "/": "root", **PATHS
    })


def predict(request):
    return JsonResponse(PATHS)


def predict_structure(request):
    pass


def predict_sequence(request):
    pass