## Run as follows: python gather_lib_alignment_qcs.py <lib_type> <code_path> <output_destination> <donor_id> <seq_run> <s3_access_file>

import re
import sys
import os
import json
from typing import List
from logger import RichLogger

lib_type = sys.argv[1]
code_path = sys.argv[2]
output_destination = sys.argv[3]
donor_id = sys.argv[4]
seq_run = sys.argv[5]
s3_access_file = sys.argv[6]

sys.path.append(code_path)
from utils import *

samples = read_immune_aging_sheet("Samples")
indices = (samples["Donor ID"] == donor_id) & (samples["Seq run"] == float(seq_run))

def add_lib(lib_type: str, all_libs: set) -> None:
    if lib_type not in ["GEX", "BCR", "TCR"]:
        raise ValueError("Unsupported lib_type: {}. Must be one of: GEX, BCR, TCR".format(lib_type))
    column_name = "{} lib".format(lib_type)
    libs_all = samples[indices][column_name]
    gex_libs_all = samples[indices]["GEX lib"]
    failed_libs = set()
    for i in range(len(libs_all)):
        if libs_all.iloc[i] is np.nan:
            continue
        libs = libs_all.iloc[i].split(",")
        gex_libs = gex_libs_all.iloc[i].split(",")
        for j in range(len(libs)):
            lib = libs[j]
            corresponding_gex_lib = gex_libs[j]
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
            all_libs.add((lib,lib_type,aligned_lib_version,corresponding_gex_lib))
    logger = RichLogger()
    for fl in failed_libs:
        logger.add_to_log("No aligned libraries found on AWS for lib id {} lib type {}. Skipping.".format(fl, lib_type), level="warning")

gex_libs = set()
ir_libs = set()
add_lib("GEX", gex_libs)
add_lib("BCR", ir_libs)
add_lib("TCR", ir_libs)

for lib in all_libs:
    lib_id = lib[0]
    lib_type = lib[1]
    aligned_lib_version = lib[2]
    corresponding_gex_lib = lib[3]
    lib_configs = {
        "sandbox_mode": "False",
        "data_owner": "valehvpa",
        "code_path": "./",
        "output_destination": "./",
        "s3_access_file": "./",
        "donor": donor_id,
        "seq_run": seq_run,
        "library_type": lib_type,
        "library_id": lib_id,
        "corresponding_gex_lib": corresponding_gex_lib,
        "filter_cells_min_genes": 600,
        "filter_cells_min_umi": 1000,
        "filter_genes_min_cells": 0,
        "filter_cells_max_pct_counts_mt": 20,
        "filter_cells_min_pct_counts_ribo": 0,
        "genes_to_exclude": "MALAT1",
        "exclude_mito_genes": "True",
        "hashsolo_priors": "0.01,0.8,0.19",
        "hashsolo_number_of_noise_barcodes": 2,
        "aligned_library_configs_version": aligned_lib_version,
        "python_env_version": "immune_aging.py_env.v4",
        "r_setup_version": "immune_aging.R_setup.v2"
    }
    filename = os.path.join(output_destination,
        "process_library.{}.{}.{}.{}.configs.txt".format(donor_id,seq_run,lib_id, lib_type))
    with open(filename, 'w') as f:
        json.dump(lib_configs, f)