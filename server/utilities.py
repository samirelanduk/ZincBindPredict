import sys
sys.path.append(".")
import requests
import os
import json
from django.http import JsonResponse
from data.utilities import model_to_residue_combos, residues_to_sample
from data.utilities import sequence_to_residue_combos, sequence_to_sample

def get_job_location(id):
    return f"server{os.path.sep}jobs{os.path.sep}{id}.json"


def save_job(job):
    with open(get_job_location(job["id"]), "w") as f:
        json.dump(job, f)


def initialise_job(id, type, protein):
    return {
        "id": id,
        "status": "initialising",
        "type": type,
        "protein": protein,
        "sites": [], "rejected": []
    }









def save_pdb_code(code, job_id):
    url = f"https://mmtf.rcsb.org/v1.0/full/{code}"
    response = requests.get(url)
    if response.status_code != 200: return None
    with open(f"server/jobs/{job_id}/{code}.mmtf", "wb") as f:
        f.write(response.content)
    return response


def save_uploaded_file(uploaded_file, job_id):
    for chunk in uploaded_file.chunks():
        with open(f"server/jobs/{job_id}/{uploaded_file.name}", "ab") as f:
            f.write(chunk)


def error_json(error):
    return JsonResponse({"error": error}, status=422, json_dumps_params={"indent": 4})


def get_job_structure_file(job_id):
    filename = None
    for f in os.listdir(f"server/jobs/{job_id}"):
        if f[-4:] in ["mmtf", ".pdb", ".cif"]:
            filename = f
            break
    return filename


def write_job(job):
    with open(f"server/jobs/{job['job_id']}/job.json", "w") as f:
        json.dump(job, f)