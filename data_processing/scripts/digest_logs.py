# This script can be used to digest a set of logs from the proess_sample.py script or the process_library.py script for a given donor,seq_run pair.
# - For sample processing logs, run as follows:
#   python digest_processing_logs.py sample <donor_id> <seq_run> <logs_location> <version> <working_dir> <s3_access_file>
# - For library processing logs, run as follows:
#   python digest_processing_logs.py library <donor_id> <seq_run> <logs_location> <version> <working_dir> <s3_access_file>

import sys
import os
import traceback
from abc import ABC, abstractmethod
from typing import List
import pandas as pd

import utils
from logger import RichLogger

class BaseDigestClass(ABC):
    def __init__(self, args: List[str]):
        self._ingest_and_sanity_check_input(args)
        self.library_type = "GEX" # we assume this for now, if this changes we can have it be passed as a command line param
        self.downloaded_from_aws = self.logs_location == "aws"

        # set aws credentials
        if self.logs_location == "aws":
            utils.set_access_keys(self.s3_access_file)

    def _ingest_and_sanity_check_input(self, args):
        assert(len(args) == 8)
        self.donor_id = args[2]
        self.seq_run = args[3]
        self.logs_location = args[4] # must be either "aws" or the absolute path to the logs location on the local disk
        self.version = args[5] # must be either the version to use (e.g. "v1") or the latest version for each sample ("latest")
        self.working_dir = args[6] # must be either "" or the absolute path to the local directory where we will place the logs downloaded from aws - only used if logs_location is "aws"
        self.s3_access_file = args[7] # must be either "" or the absolute path to the aws credentials file - only used if logs_location is "aws"

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

    def digest_logs(self):
        logger = RichLogger()
        try:
            object_ids = self._get_object_ids()
            
            # if the logs location is aws, download all log files to the local working directory
            if self.logs_location == "aws":
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

            # for each object id, parse its logs and report any noteworthy log events
            log_lines_to_print = {}
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
                    # for now, we only print logs above a hard coded severity threshold (see `_is_alertable_log_line`),
                    # but we can extend this script to take a user-provided threshold 
                    if self._is_alertable_log_line(line):
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
                        utils.draw_separator_line()
                    file_name = key.split("/")[-1]
                    print(file_name + ":\n")
                    for line in value:
                        print("\t" + line)
                    first_item = False

            # clean up the logs we downloaded from aws if any
            if self.downloaded_from_aws:
                self._remove_logs(self.logs_location)
        except Exception as err:
            logger.add_to_log("Execution failed with the following error:\n{}".format(traceback.format_exc()), "critical")
            # clean up the logs we downloaded from aws if any
            if self.downloaded_from_aws:
                self._remove_logs(self.logs_location)
            sys.exit()


class DigestSampleProcessingLogs(BaseDigestClass):
    def __init__(self, args: List[str]):
        super().__init__(args)

    def _get_object_ids(self):
        samples_df = self._get_all_samples()
        object_ids = samples_df["Sample_ID"]
        return object_ids

    def _get_object_prefix(self, object_id: str):
        return "{}_{}".format(object_id, self.library_type)

    def _get_log_file_name(self, object_id: str):
        prefix = self._get_object_prefix(object_id)
        return "process_sample.{}.{}.log".format(prefix, self.version)

    def _get_aws_dir_name(self):
        return "processed_samples"


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
        return "{}_{}_{}_{}".format(self.donor, self.seq_run, self.library_type, object_id)

    def _get_log_file_name(self, object_id: str):
        prefix = self._get_object_prefix()
        return "process_library.{}.{}.log".format(prefix, self.version)

    def _get_aws_dir_name(self):
        return "processed_libraries"


process_type = sys.argv[1] # must be one of "sample" or "library"
assert(process_type in ["sample", "library"])
if process_type == "sample":
    digest_class = DigestSampleProcessingLogs(sys.argv)
else:
    digest_class = DigestLibraryProcessingLogs(sys.argv)
digest_class.digest_logs()
