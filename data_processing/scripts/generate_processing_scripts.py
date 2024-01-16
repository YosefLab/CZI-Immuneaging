## This script generates .sh files with commands for running the current processing jobs in the job queue on S3 (one script for running the library processing jobs in the queue and one for the sample processing jobs).
## Example: python generate_processing_scripts.py <aws_credentials_file> /data/yosef2/scratch/immuneaging/processing_results /data/yosef2/users/erahmani/projects/immune_aging/code /data/yosef2/scratch/immuneaging/jobs_sh

import sys
import os
import json
import subprocess
from pathlib import Path

s3_access_file = sys.argv[1]
output_destination = sys.argv[2]
code_path = sys.argv[3]
output_path = sys.argv[4]

sys.path.append("code_path")
from utils import *
set_access_keys(s3_access_file)

os.system("mkdir -p " + output_destination)
jobs_queue_destination = os.path.join(output_destination,"job_queue")
os.system("mkdir -p " + jobs_queue_destination)
jobs_run_destination = os.path.join(output_destination,"job_run")
os.system("mkdir -p " + jobs_run_destination)

code_version = subprocess.check_output(
    ["git", "describe", "--always"], cwd=Path(code_path).resolve().parent).strip().decode()
print(code_version)

outfiles = []

for job_type in ("process_library", "process_sample"):
    sh_cmds = {}
    ls_cmd = 'aws s3 ls s3://immuneaging/job_queue/{}/'.format(job_type)
    aws_response = os.popen(ls_cmd).read()
    if len(aws_response):
        sync_cmd = 'aws s3 cp --recursive --no-progress s3://immuneaging/job_queue/{}/ {}'.format(job_type, jobs_queue_destination)
        print("Running aws command: " + sync_cmd)
        print("output:\n" + os.popen(sync_cmd).read())
    else:
        print("No jobs of type " + job_type)
        next
    for i in aws_response.split("\n"):
        if "configs" in i:
            configs_filename = i.split(" ")[-1]
            queue_filename = os.path.join(jobs_queue_destination,configs_filename)
            with open(queue_filename, "r") as f:
                configs = json.loads(f.read())
            run_filename = os.path.join(jobs_run_destination, configs_filename)
            configs["code_path"] = code_path
            configs["output_destination"] = output_destination
            configs["s3_access_file"] = s3_access_file
            configs["code_version"] = code_version
            if configs["donor"] not in sh_cmds:
                sh_cmds[configs["donor"]] = []
            # save configs file
            with open(run_filename, 'w') as f:
                json.dump(configs, f)
            # add commands for running the current job
            sh_cmds[configs["donor"]].append('{}'.format(run_filename))
    # move the original file on S3 from the queue to the running folder
    mv_cmd = 'aws s3 mv --recursive s3://immuneaging/job_queue/{0}/ s3://immuneaging/job_queue/{0}.running/ --exclude "*" --include "{0}*"'.format(job_type)
    os.system(mv_cmd)
    # save the commands in sh_cmds, for each donor separately
    if len(sh_cmds):
        for i in sh_cmds:
            if not os.path.exists(os.path.join(output_path,"{}".format(job_type))):
                os.mkdir(os.path.join(output_path,"{}".format(job_type)))
            outfiles.append(os.path.join(output_path,"{}/{}.{}_jobs.sh".format(job_type,i,job_type)))
            with open(outfiles[-1], 'w') as f:
                for cmd in sh_cmds[i]:
                    f.write(cmd+"\n")

if (len(outfiles)):
    print("Generated the following files:")
    print(outfiles)
else:
    print("There are no jobs to run. Generated no files.")