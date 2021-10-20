## This script can be used to digest a set of logs from the proess_sample.py script for a given donor,seq_run pair.
## Run as follows: python digest_process_sample_logs.py <donor_id> <seq_run> <logs_location> <version> <working_dir> <s3_access_file>

import sys
import os
import math
import traceback

import utils
from logger import RichLogger, is_alertable_log_line

def get_process_sample_log_file_name(sample_id: str, library_type: str, version: str):
    prefix = "{}_{}".format(sample_id, library_type)
    return "process_sample.{}.{}.log".format(prefix, version)

def remove_logs(logs_dir: str):
    if os.path.isdir(logs_dir):
        os.system("rm {}/*".format(logs_dir))

donor_id = sys.argv[1]
seq_run = sys.argv[2]
logs_location = sys.argv[3] # must be either "aws" or the absolute path to the logs location on the local disk
version = sys.argv[4] # must be either the version to use (e.g. "v1") or the latest version for each sample ("latest")
working_dir = sys.argv[5] # must be either "" or the absolute path to the local directory where we will place the logs downloaded from aws - only used if logs_location is "aws"
s3_access_file = sys.argv[6] # must be either "" or the absolute path to the aws credentials file - only used if logs_location is "aws"

# sanity check command input
assert(logs_location=="aws" or os.path.isdir(logs_location))
assert(working_dir == "" if logs_location!="aws" else os.path.isdir(working_dir))
assert(s3_access_file == "" if logs_location!="aws" else os.path.isfile(s3_access_file))

# set aws credentials
if logs_location == "aws":
    utils.set_access_keys(s3_access_file)

# gather the full set of sample id's we are interested in
samples = utils.read_immune_aging_sheet("Samples", quiet=True)
indices = samples["Donor ID"] == donor_id
sample_ids = samples[indices]["Sample_ID"]

library_type = "GEX" # we assume this for now, if this changes we can have it be passed as a command line param
downloaded_from_aws = logs_location == "aws"

try:
    # if the logs location is aws, download all log files to the local working directory
    logger = RichLogger()
    if logs_location == "aws":
        logs_location = working_dir
        for sample_id in sample_ids:
            prefix = "{}_{}".format(sample_id, library_type)
            filename = get_process_sample_log_file_name(sample_id, library_type, version)
            sync_cmd = 'aws s3 sync --no-progress s3://immuneaging/processed_samples/{}/{} {} --exclude "*" --include {}'.format(
                prefix, version, working_dir, filename)
            logger.add_to_log("syncing {}...".format(filename))
            logger.add_to_log("sync_cmd: {}".format(sync_cmd))
            logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

    if version == "latest":
        # TODO issue #15
        raise NotImplementedError

    # for each sample id, parse its process_sample logs and report any noteworthy log events
    log_lines_to_print = {}
    for sample_id in sample_ids:
        filename = get_process_sample_log_file_name(sample_id, library_type, version)
        filepath = os.path.join(logs_location, filename)
        if not os.path.isfile(filepath):
            logger.add_to_log("File not found: {}. Skipping.".format(filepath))
            continue
        lines = []
        with open(filepath, 'r') as f:
            lines = f.readlines()
        for line in lines:
            # for now, we only print logs above a hard coded severity threshold (see `is_alertable_log_line`),
            # but we can extend this script to take a user-provided threshold 
            if is_alertable_log_line(line):
                if filepath not in log_lines_to_print:
                    log_lines_to_print[filepath] = [line]
                else:
                    log_lines_to_print[filepath].append(line)

    if len(log_lines_to_print) == 0:
        logger.add_to_log("No relevant log lines were found.", "info")
    else:
        logger.add_to_log("Found the following relevant log lines.", "info")
        first_item = True
        for key,value in log_lines_to_print.items():
            if not first_item:
                # draw a separator line between files
                width = os.get_terminal_size().columns / 5
                print(" " * math.floor(width)  + "\u2014" * 3 * math.floor(width) + " " * math.floor(width) + "\n")
            file_name = key.split("/")[-1]
            print(file_name + ":\n")
            for line in value:
                print("\t" + line)
            first_item = False

    # clean up the logs we downloaded from aws if any
    if downloaded_from_aws:
        remove_logs(logs_location)
except Exception as err:
    logger.add_to_log("Execution failed with the following error:\n{}".format(traceback.format_exc()), "critical")
    # clean up the logs we downloaded from aws if any
    if downloaded_from_aws:
        remove_logs(logs_location)
    sys.exit()
