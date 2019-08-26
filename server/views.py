import atomium
import os
import time
import json
from datetime import datetime
from subprocess import Popen
import requests
from django.http import JsonResponse
from django.conf import settings
from .utilities import model_to_residue_combos, residues_to_sample

PATHS = {
 "/predict/structure?code=XXXX": "Find zinc binding sites in PDB structure"
}

def root(request):
    return JsonResponse({
     "/": "root", **PATHS
    }, json_dumps_params={"indent": 4})


def predict(request):
    return JsonResponse(PATHS, json_dumps_params={"indent": 4})


def predict_structure(request):
    job_id = int(time.time() * 1000)
    os.mkdir(f"server/jobs/{job_id}")

    if "code" in request.GET:
        code = request.GET["code"]
        url = f"https://mmtf.rcsb.org/v1.0/full/{code}"
        response = requests.get(url)
        if response.status_code != 200:
            return JsonResponse({
             "error": f"{code} does not seem to be a valid PDB"
            }, status=422, json_dumps_params={"indent": 4})
        
        
        with open(f"server/jobs/{job_id}/{code}.mmtf", "wb") as f:
            f.write(response.content)
    
    elif "file" in request.FILES:

        uploaded_file = request.FILES["file"]
        for chunk in uploaded_file.chunks():
            with open(f"server/jobs/{job_id}/{uploaded_file.name}", "ab") as f:
                f.write(chunk)


    else:
        return JsonResponse({
         "error": "You must provide a PDB code or a structure file"
        }, status=422, json_dumps_params={"indent": 4})
    
    p = Popen(["server/job.py", str(job_id)])
    s = "s" if settings.ALLOWED_HOSTS else ""
    return JsonResponse({
     "job": f"http{s}://{request.META['HTTP_HOST']}{request.path}{job_id}",
     "expires": datetime.fromtimestamp(
      job_id / 1000 + (settings.JOB_EXPIRATION * 60 * 60 * 24)
     ).strftime("%-d %B %Y, %H:%M UTC")
    }, json_dumps_params={"indent": 2})


def job(request, id):
    with open(f"server/jobs/{id}/job.json") as f:
        data = json.load(f)
    return JsonResponse(data, json_dumps_params={"indent": 2})


def predict_sequence(request):
    pass