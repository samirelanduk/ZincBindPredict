import sys
import os

# What is the job ID
if len(sys.argv) < 2:
    print("Please provide Job ID")
    sys.exit(1)
job_id = sys.argv[1]

# Make sure there is a folder for this job
if not os.path.exists(f"server/jobs/{job_id}"):
    print("There's no job with that ID")
    sys.exit(1)
