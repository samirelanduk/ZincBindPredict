import sys
sys.path.append(".")
import os
import json
from data.utilities import model_to_residue_combos, residues_to_sample

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