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
    })


def predict(request):
    return JsonResponse(PATHS)


def predict_structure(request):
    if "code" not in request.GET:
        return JsonResponse({
         "error": "You must provide a PDB code"
        }, status=422)
    code = request.GET["code"]
    url = f"https://mmtf.rcsb.org/v1.0/full/{code}"
    response = requests.get(url)
    if response.status_code != 200:
        return JsonResponse({
         "error": f"{code} does not seem to be a valid PDB"
        }, status=422)
    job_id = int(time.time() * 1000)
    os.mkdir(f"server/jobs/{job_id}")
    with open(f"server/jobs/{job_id}/{code}.mmtf", "wb") as f:
        f.write(response.content)
    p = Popen(["server/job.py", str(job_id)])
    s = "s" if settings.ALLOWED_HOSTS else ""
    return JsonResponse({
     "job": f"http{s}://{request.META['HTTP_HOST']}{request.path}{job_id}",
     "expires": datetime.fromtimestamp(
      job_id / 1000 + (settings.JOB_EXPIRATION * 60 * 60 * 24)
     ).strftime("%-d %B %Y, %H:%M UTC")
    })


def job(request, id):
    with open(f"server/jobs/{id}/job.json") as f:
        data = json.load(f)
    return JsonResponse(data)


def predict_sequence(request):
    pass