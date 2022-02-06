## Run as follows: python gather_lib_alignment_qcs.py <code_path> <output_destination> <donor_id> <seq_run> <s3_access_file>

import re
import sys
import os
import json
from typing import List
from logger import RichLogger

code_path = sys.argv[1]
output_destination = sys.argv[2]
donor_id = sys.argv[3]
seq_run = sys.argv[4]
s3_access_file = sys.argv[5]

sys.path.append(code_path)
from utils import *

samples = read_immune_aging_sheet("Samples")
indices = (samples["Donor ID"] == donor_id) & (samples["Seq run"] == float(seq_run))

logger = RichLogger()

# TODO crawl through all donors and merge them in all in a single csv per lib type

def add_lib(lib_type: str, all_libs: set) -> None:
    if lib_type not in ["GEX", "BCR", "TCR"]:
        raise ValueError("Unsupported lib_type: {}. Must be one of: GEX, BCR, TCR".format(lib_type))
    column_name = "{} lib".format(lib_type)
    libs_all = samples[indices][column_name]
    failed_libs = set()
    for i in range(len(libs_all)):
        if libs_all.iloc[i] is np.nan:
            continue
        libs = libs_all.iloc[i].split(",")
        for lib in libs:
            # find the latest aligned_lib_version
            set_access_keys(s3_access_file)
            latest_version = -1
            ls_cmd = "aws s3 ls s3://immuneaging/aligned_libraries --recursive"
            ls  = os.popen(ls_cmd).read()
            if len(ls) != 0:
                filenames = ls.split("\n")
                if lib_type == "GEX":
                    # search for patterns of <lib id>.vX.h5ad. If there is a match, group
                    # one is "<lib id>.v" and group 2 is "X" (X can be any integer >0)
                    pattern = "({}\.v)(\d+)\.h5ad$".format(lib)
                else:
                    pattern = "({}_{}\.cellranger\.filtered_contig_annotations\.v)(\d+)\.csv".format(lib_type, lib)
                for filename in filenames:
                    m = re.search(pattern, filename)
                    if bool(m):
                        version = int(m[2])
                        if latest_version < version:
                            latest_version = version
            if latest_version == -1:
                failed_libs.add(lib)
                # skip
                continue
            aligned_lib_version = "v" + str(latest_version)
            all_libs.add((lib,lib_type,aligned_lib_version))
    for fl in failed_libs:
        logger.add_to_log("No aligned libraries found on AWS for lib id {} lib type {}. Skipping.".format(fl, lib_type), level="warning")

def combine_metrics_for_lib(libs: set, generic_lib_type: str, data_dir: str):
    lib_data_dir = os.path.join(data_dir, generic_lib_type)
    os.system("mkdir -p " + lib_data_dir)
    all_metrics_files = []
    for lib in libs:
        lib_id = lib[0]
        lib_type = lib[1]
        aligned_lib_version = lib[2]
        logger.add_to_log("Downloading metrics.csv file for lib id {}, lib type {} from S3...".format(lib_id, lib_type))
        metrics_csv_file_name = "{}_{}_{}_{}.cellranger.metrics_summary.csv".format(donor_id, seq_run, lib_type, lib_id)
        sync_cmd = 'aws s3 sync --no-progress s3://immuneaging/aligned_libraries/{}/{}_{}_{}_{}/ {} --exclude "*" --include {}'.format(
            aligned_lib_version, donor_id, seq_run, lib_type, lib_id, lib_data_dir, metrics_csv_file_name
        )
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
        metrics_csv_file = os.path.join(lib_data_dir, metrics_csv_file_name)
        if not os.path.isfile(metrics_csv_file):
            msg = "Failed to download file {} from S3.".format(metrics_csv_file)
            logger.add_to_log(msg, level="error")
            raise ValueError(msg)
        all_metrics_files.append((metrics_csv_file,lib_id,lib_type,aligned_lib_version))

    # create the combined csv metrics file
    all_dfs = []
    for f in all_metrics_files:
        df = pd.read_csv(f[0])
        df["Lib Id"] = f[1]
        df["Lib Type"] = f[2]
        df["Aligned Lib Version"] = f[3]
        # TODO del any df columns that we don't want, if we only desire to keep a subset
        all_dfs.append(df)
    combined_df = pd.concat(all_dfs, ignore_index=True)
    combined_metrics = os.path.join(data_dir, "{}_{}_all_{}_metrics.csv".format(donor_id,seq_run,generic_lib_type))
    with open(combined_metrics, 'w') as f:
        combined_df.to_csv(f)
    # upload the combined csv file to AWS
    logger.add_to_log("Uploading combined metrics file {} to S3...".format(combined_metrics.split("/")[-1]))
    sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/combined_lib_alignment_metrics/{}/ --exclude "*" --include {}'.format(data_dir, generic_lib_type, combined_metrics.split("/")[-1])
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

gex_libs = set()
ir_libs = set()
add_lib("GEX", gex_libs)
add_lib("BCR", ir_libs)
add_lib("TCR", ir_libs)

# create a new directory for the outputs
data_dir = os.path.join(output_destination, "_".join([donor_id, seq_run]))
os.system("mkdir -p " + data_dir)

combine_metrics_for_lib(gex_libs, "GEX", data_dir)
combine_metrics_for_lib(ir_libs, "IR", data_dir)