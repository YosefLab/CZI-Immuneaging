## This script can be used to digest a set of logs from the proess_sample.py script for a given donor,seq_run pair.
## Run as follows: TODO

import sys
import os
import math
import utils
import logger

donor_id = sys.argv[1]
seq_run = sys.argv[2]
logs_location = sys.argv[3] # must be one of "aws" or path to the logs location on the local disk
version = sys.argv[4] # must be either the version to use (e.g. "v1") or the latest version for each sample ("latest")
working_dir = sys.argv[5] # must be either "" or the path to the local directory where we will temporarily place the logs downloaded from aws - only used if logs_location is "aws"

# sanity check command input
assert(logs_location=="aws" or os.path.isdir(logs_location))
assert(os.path.isdir(working_dir))

# gather the full set of sample id's we are interested in
samples = utils.read_immune_aging_sheet("Samples")
indices = samples["Donor ID"] == donor_id
sample_ids = samples[indices]["Sample_ID"]

# if the logs location is aws, download all log files to the local working directory
if logs_location == "aws":
    # TODO download logs from aws to a local working directory
    raise NotImplementedError

if version == "latest":
    # TODO
    raise NotImplementedError

# for each sample id, parse its process_sample logs and report any noteworthy log events
log_lines_to_print = {}
for sample_id in sample_ids:
    filename = os.path.join(logs_location, "process_sample.{}.{}.log".format(sample_id, version))
    lines = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        # for now, we only print logs above a hard coded level threshold (ERROR),
        # but we can extend this script to take a user-provided threshold 
        if logger.is_alertable_log_line(line):
            if filename not in log_lines_to_print:
                log_lines_to_print[filename] = [line]
            else:
                log_lines_to_print[filename].append(line)

rich_logger = logger.RichLogger()
if len(log_lines_to_print) == 0:
    rich_logger.add_to_log("No relevant log lines were found.", "info")
else:
    rich_logger.add_to_log("Found the following relevant log lines.", "info")
    first_item = True
    for item in log_lines_to_print:
        if not first_item:
            # draw a separator line between files
            width = os.get_terminal_size().columns / 5
            print(" " * math.floor(width)  + "\u2014" * 3 * math.floor(width) + " " * math.floor(width))
        print(item.key + ":\n")
        for line in item.value:
            print("\t\t" + line)
        first_item = False
