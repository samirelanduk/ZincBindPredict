import atomium
import os
import time
import json
from datetime import datetime
from subprocess import Popen
import requests
from django.http import JsonResponse
from django.conf import settings
from .utilities import *

PATHS = {
 "/structure?code=XXXX": "Find zinc binding sites in PDB structure",
 "/structure/{ID}/": "The results of a structure job",
 "/sequence?sequence=XXXXXXXXXX": "Find zinc binding sites in peptide sequence",
 "/sequence/{ID}/": "The results of a sequence job"
}

def root(request):
    return JsonResponse({
     "/": "root", **PATHS
    }, json_dumps_params={"indent": 4})


def predict_structure(request):
    job_id = int(time.time() * 1000)
    os.mkdir(f"server/jobs/{job_id}")
    if "code" in request.GET:
        if not save_pdb_code(request.GET["code"], job_id):
            return error_json("That does not seem to be a valid PDB code")
    elif "file" in request.FILES:
        save_uploaded_file(request.FILES["file"], job_id)
    else:
        return error_json(f"You must provide a PDB code or a structure file")
    p = Popen(["server/structure_job.py", str(job_id)])
    s = "s" if settings.ALLOWED_HOSTS else ""
    return JsonResponse({
     "job": f"http{s}://{request.META['HTTP_HOST']}{request.path}{job_id}",
     "expires": datetime.fromtimestamp(
      job_id / 1000 + (settings.JOB_EXPIRATION * 60 * 60 * 24)
     ).strftime("%-d %B %Y, %H:%M UTC")
    }, json_dumps_params={"indent": 2})


def predict_sequence(request):
    job_id = int(time.time() * 1000)
    os.mkdir(f"server/jobs/{job_id}")
    if "sequence" not in request.GET:
        return error_json(f"You must provide a sequence")
    p = Popen(["server/sequence_job.py", str(job_id), request.GET["sequence"]])
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