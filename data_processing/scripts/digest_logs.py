# This script can be used to digest a set of logs from different processing scripts (currently supported: process_sample logs and process_library logs) for a given donor,seq_run pair. Various actions can be requested, for example: print noteworthy log lines, get a digest in csv, etc.
# - For sample processing logs, run as follows:
#   python digest_logs.py <action> sample <donor_id> <seq_run> <logs_location> <version> <working_dir> <s3_access_file>
# - For library processing logs, run as follows:
#   python digest_logs.py <action> library <donor_id> <seq_run> <logs_location> <version> <working_dir> <s3_access_file>

import sys
import os
import traceback
from abc import ABC, abstractmethod
from typing import List, Dict
import pandas as pd
import csv
from parse import *
import logging
import io

import utils
from logger import RichLogger

logging.getLogger('parse').setLevel(logging.WARNING)

class BaseDigestClass(ABC):
    def __init__(self, args: List[str]):
        self._ingest_and_sanity_check_input(args)
        self.library_type = "GEX" # we assume this for now, if this changes we can have it be passed as a command line param
        self.downloaded_from_aws = self.logs_location == "aws"

        # set aws credentials
        if self.logs_location == "aws":
            utils.set_access_keys(self.s3_access_file)

    def _ingest_and_sanity_check_input(self, args):
        assert(len(args) == 9)
        self.donor_id = args[3]
        self.seq_run = args[4]
        self.logs_location = args[5] # must be either "aws" or the absolute path to the logs location on the local disk
        self.version = args[6] # must be either the version to use (e.g. "v1") or the latest version for each sample ("latest")
        self.working_dir = args[7] # must be either "" or the absolute path to the local directory where we will place the logs downloaded from aws - only used if logs_location is "aws"
        self.s3_access_file = args[8] # must be either "" or the absolute path to the aws credentials file - only used if logs_location is "aws"

        assert(self.logs_location == "aws" or os.path.isdir(self.logs_location))
        assert(self.working_dir == "" if self.logs_location != "aws" else os.path.isdir(self.working_dir))
        assert(self.s3_access_file == "" if self.logs_location != "aws" else os.path.isfile(self.s3_access_file))

    def _get_all_samples(self) -> pd.DataFrame:
        samples = utils.read_immune_aging_sheet("Samples", quiet=True)
        indices = samples["Donor ID"] == self.donor_id
        return samples[indices]

    @abstractmethod
    def _get_object_ids(self):
        """
        Gathers the full set of object id's we are interested in. These id's can then be used to get the corresponding
        log file name (e.g. using ``_get_log_file_name``) to download logs from.
        """
        pass

    @abstractmethod
    def _get_object_prefix(self, object_id: str):
        pass
        
    @abstractmethod
    def _get_log_file_name(self, object_id: str):
        pass

    @abstractmethod
    def _get_aws_dir_name(self):
        pass

    @staticmethod
    def _remove_logs(logs_dir: str):
        if os.path.isdir(logs_dir):
            os.system("rm {}/*".format(logs_dir))

    @staticmethod
    def _is_alertable_log_line(line: str) -> bool:
        return "WARNING" in line or "ERROR" in line or "CRITICAL" in line

    def _get_log_lines(self) -> Dict[str, List]:
        files_to_lines = {}
        logger = RichLogger()
        try:
            object_ids = self._get_object_ids()
            
            # if the logs location is aws, download all log files to the local working directory
            if self.logs_location == "aws":
                # set this first so that if we hit an exception mid-way through the loop below
                # we can still attempt to clean up the downloaded logs
                self.logs_location = self.working_dir
                for object_id in object_ids:
                    prefix = self._get_object_prefix(object_id)
                    filename = self._get_log_file_name(object_id)
                    aws_dir_name = self._get_aws_dir_name()
                    sync_cmd = 'aws s3 sync --no-progress s3://immuneaging/{}/{}/{} {} --exclude "*" --include {}'.format(aws_dir_name, prefix, self.version, self.working_dir, filename)
                    logger.add_to_log("syncing {}...".format(filename))
                    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
                    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

            if self.version == "latest":
                # TODO issue #15
                raise NotImplementedError

            # for each object id, get all its log lines and add it to the dict
            for object_id in object_ids:
                filename = self._get_log_file_name(object_id)
                filepath = os.path.join(self.logs_location, filename)
                if not os.path.isfile(filepath):
                    logger.add_to_log("File not found: {}. Skipping.".format(filepath))
                    continue
                lines = []
                with open(filepath, 'r') as f:
                    lines = f.readlines()
                for line in lines:
                    if filepath not in files_to_lines:
                        files_to_lines[filepath] = [line]
                    else:
                        files_to_lines[filepath].append(line)

            # clean up the logs we downloaded from aws if any
            if self.downloaded_from_aws:
                self._remove_logs(self.logs_location)
        except Exception:
            logger.add_to_log("Execution failed with the following error:\n{}".format(traceback.format_exc()), "critical")
            # clean up the logs we downloaded from aws if any
            if self.downloaded_from_aws:
                self._remove_logs(self.logs_location)
            sys.exit()
        return files_to_lines

    def print_digest(self):
        logger = RichLogger()
        try:
            files_to_lines = self._get_log_lines()
            
            # for each file, parse its logs and report any noteworthy log events
            log_lines_to_print = {}
            for filepath,lines in files_to_lines.items():
                for line in lines:
                    # for now, we only print logs above a hard coded severity threshold (see `_is_alertable_log_line`),
                    # but we can extend this script to take a user-provided threshold 
                    if self._is_alertable_log_line(line):
                        if filepath not in log_lines_to_print:
                            log_lines_to_print[filepath] = [line]
                        else:
                            log_lines_to_print[filepath].append(line)

            # print digested log lines
            if len(log_lines_to_print) == 0:
                logger.add_to_log("No relevant log lines were found.", "info")
            else:
                logger.add_to_log("Found the following relevant log lines.", "info")
                first_item = True
                for key,value in log_lines_to_print.items():
                    if not first_item:
                        # draw a separator line between files
                        utils.draw_separator_line()
                    file_name = key.split("/")[-1]
                    print(file_name + ":\n")
                    for line in value:
                        print("\t" + line)
                    first_item = False
        except Exception:
            logger.add_to_log("Execution failed with the following error:\n{}".format(traceback.format_exc()), "critical")
            sys.exit()

    @abstractmethod
    def get_digest_csv(self):
        pass

class DigestSampleProcessingLogs(BaseDigestClass):
    def __init__(self, args: List[str]):
        super().__init__(args)

    def _get_object_ids(self):
        samples_df = self._get_all_samples()
        return samples_df["Sample_ID"]

    def _get_object_prefix(self, object_id: str):
        return "{}_{}".format(object_id, self.library_type)

    def _get_log_file_name(self, object_id: str):
        prefix = self._get_object_prefix(object_id)
        return "process_sample.{}.{}.log".format(prefix, self.version)

    def _get_object_id(self, log_file_name: str):
        # log file name is process_sample.{prefix}.{version}.log where prefix is {object_id}_{library_type}
        prefix = log_file_name.split(".")[1]
        return prefix.split("_")[0]

    def _get_aws_dir_name(self):
        return "processed_samples"

    def get_digest_csv(self):
        logger = RichLogger()
        try:
            files_to_lines = self._get_log_lines()

            # define the csv digest headers
            CSV_HEADER_SAMPLE_ID: str = "Sample ID"
            CSV_HEADER_FAILED: str = "Failed?"
            CSV_HEADER_WARNING: str = "Warning?"
            CSV_HEADER_FAILURE_REASON: str = "Failure Reason"
            CSV_HEADER_WARNING_REASON: str = "Warning Reason"
            CSV_HEADER_DOUBLETS: str = "% doublets"
            CSV_HEADER_AMBIENT_RNA: str = "% ambient RNA"
            CSV_HEADER_VDJ: str = "% vdj genes"
            CSV_HEADER_RBC: str = "% RBC"

            # for each file, parse its logs and add digest info to csv
            csv_rows = []
            for filepath,lines in files_to_lines.items():
                csv_row = {
                    CSV_HEADER_SAMPLE_ID: self._get_object_id(filepath),
                    CSV_HEADER_FAILED: "No",
                    CSV_HEADER_WARNING: "No",
                    CSV_HEADER_FAILURE_REASON: "",
                    CSV_HEADER_WARNING_REASON: "",
                    CSV_HEADER_DOUBLETS: -1,
                    CSV_HEADER_AMBIENT_RNA: -1,
                    CSV_HEADER_VDJ: -1,
                    CSV_HEADER_RBC: -1,
                }
                for line in lines:
                    if "ERROR" in line or "CRITICAL" in line:
                        csv_row[CSV_HEADER_FAILED] = "Yes"
                        csv_row[CSV_HEADER_FAILURE_REASON] = line.strip('"').strip()
                    if "WARNING" in line:
                        csv_row[CSV_HEADER_WARNING] = "Yes"
                        if csv_row[CSV_HEADER_WARNING_REASON] == "": # first one
                            csv_row[CSV_HEADER_WARNING_REASON] += line.strip('"').strip()
                        else:
                            csv_row[CSV_HEADER_WARNING_REASON] += " --- " + line.strip('"').strip()
                    # doublets
                    # TODO update env with parse
                    doublets_qc_parsed = search(utils.QC_STRING_DOUBLETS, line)
                    if doublets_qc_parsed:
                        csv_row[CSV_HEADER_DOUBLETS] = doublets_qc_parsed[1]
                    # ambient rna
                    ambient_rna_qc_parsed = search(utils.QC_STRING_AMBIENT_RNA, line)
                    if ambient_rna_qc_parsed:
                        csv_row[CSV_HEADER_AMBIENT_RNA] = ambient_rna_qc_parsed[1]
                    # vdj genes
                    vdj_qc_parsed = search(utils.QC_STRING_VDJ, line)
                    if vdj_qc_parsed:
                        csv_row[CSV_HEADER_VDJ] = vdj_qc_parsed[1]
                    # red blood cells
                    rbc_qc_parsed = search(utils.QC_STRING_RBC, line)
                    if rbc_qc_parsed:
                        csv_row[CSV_HEADER_RBC] = rbc_qc_parsed[1]
                csv_rows.append(csv_row)

            # write the csv
            csv_file = io.StringIO()
            field_names = [
                CSV_HEADER_SAMPLE_ID,
                CSV_HEADER_FAILED,
                CSV_HEADER_WARNING,
                CSV_HEADER_FAILURE_REASON,
                CSV_HEADER_WARNING_REASON,
                CSV_HEADER_DOUBLETS,
                CSV_HEADER_AMBIENT_RNA,
                CSV_HEADER_VDJ,
                CSV_HEADER_RBC,
            ]
            writer = csv.DictWriter(csv_file, fieldnames=field_names)
            writer.writeheader()
            writer.writerows(csv_rows)
            print(csv_file.getvalue())
            csv_file.close()
        except Exception:
            logger.add_to_log("Execution failed with the following error:\n{}".format(traceback.format_exc()), "critical")
            sys.exit()

class DigestLibraryProcessingLogs(BaseDigestClass):
    def __init__(self, args: List[str]):
        super().__init__(args)

    def _get_object_ids(self):
        object_ids = set()
        samples_df = self._get_all_samples()
        for i in samples_df[self.library_type + " lib"]:
            for j in i.split(","):
                object_ids.add(j)
        return object_ids

    def _get_object_prefix(self, object_id: str):
        return "{}_{}_{}_{}".format(self.donor_id, self.seq_run, self.library_type, object_id)

    def _get_log_file_name(self, object_id: str):
        prefix = self._get_object_prefix(object_id)
        return "process_library.{}.{}.log".format(prefix, self.version)

    def _get_aws_dir_name(self):
        return "processed_libraries"

    def get_digest_csv(self):
        raise NotImplementedError


process_type = sys.argv[2] # must be one of "sample" or "library"
action = sys.argv[1] # must be one of "sample" or "library"
assert(process_type in ["sample", "library"])
assert(action in ["print_digest", "get_csv"])
if process_type == "sample":
    digest_class = DigestSampleProcessingLogs(sys.argv)
else:
    digest_class = DigestLibraryProcessingLogs(sys.argv)
if action == "print_digest":
    digest_class.print_digest()
else:
    digest_class.get_digest_csv()
