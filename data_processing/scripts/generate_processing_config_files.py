## This script can be used to generate config files for all libraries and all samples from a given donor; currently only GEX libraries are considered (i.e. ignores BCR/TCR libs).
## Note that 
## Run as follows: python generate_processing_config_files.py <code_path> <output_destination> <donor_id> <seq_run>
## Note that code_path, output_destination, and s3_access_file will be set later via generate_processing_scripts.py
## after generating the files run the following commands to upload to aws (after setting the aws credentials as env variables):
## aws s3 sync config_files s3://immuneaging/job_queue/process_library/ --exclude "*" --include "process_library*.configs.txt"
## aws s3 sync config_files s3://immuneaging/job_queue/process_sample/ --exclude "*" --include "process_sample.configs*.txt"
## assuming config_files is the directory containign the config files.
## Then, the script generate_processing_scripts.py can be used to create .sh files to execute the actual processing based on the config files.

import sys
import os
import json
import subprocess

code_path = sys.argv[1]
output_destination = sys.argv[2]
donor_id = sys.argv[3]
seq_run = sys.argv[4]

sys.path.append(code_path)
from utils import *

samples = read_immune_aging_sheet("Samples")
indices = samples["Donor ID"] == donor_id

# create config files for library processing
library_ids = set()
for i in samples[indices]["GEX lib"]:
    for j in i.split(","):
        library_ids.add(j)

for library_id in library_ids:
    lib_configs = {
        "sandbox_mode": "False",
        "data_owner": "erahmani",
        "code_path": "./",
        "output_destination": "./",
        "s3_access_file": "./",
        "donor": donor_id,
        "seq_run": seq_run,
        "library_type": "GEX",
        "library_id": library_id,
        "filter_cells_min_genes": 200,
        "filter_genes_min_cells": 0,
        "filter_cells_max_pct_counts_mt": 20,
        "filter_cells_min_pct_counts_ribo": 5,
        "genes_to_exclude": "MALAT1",
        "exclude_mito_genes": "True",
        "hashsolo_priors": "0.01,0.8,0.19",
        "hashsolo_number_of_noise_barcodes": 2,
        "aligned_library_configs_version": "v1",
        "python_env_version": "immune_aging.py_env.v2",
        "r_setup_version": "immune_aging.R_setup.v2"
    }
    filename = os.path.join(output_destination,
        "process_library.{}.{}.{}.configs.txt".format(donor_id,seq_run,library_id))
    with open(filename, 'w') as f:
                json.dump(lib_configs, f)


# create config files for sample processing

sample_ids = samples[indices]["Sample_ID"]

for sample_id in sample_ids:
    library_ids = [i for i in samples[samples["Sample_ID"] == sample_id]["GEX lib"].iloc[0].split(",")]
    processed_library_configs_version = ["v1" for i in range(len(library_ids))]
    sample_configs = {
        "sandbox_mode": "False",
        "data_owner": "erahmani",
        "code_path": "./",
        "output_destination": "./",
        "s3_access_file": "./",
        "donor": donor_id,
        "seq_run": seq_run,
        "library_type": "GEX",
        "sample_id": sample_id,
        "library_ids": ",".join(library_ids),
        "processed_library_configs_version": ",".join(processed_library_configs_version),
        "min_cells_per_library": 200,
        "filter_decontaminated_cells_min_genes": 100,
        "normalize_total_target_sum": 10000,
        "n_highly_variable_genes": 3000,
        "highly_variable_genes_flavor": "seurat_v3",
        "scvi_max_epochs": 400,
        "totalvi_max_epochs": 400,
        "solo_filter_genes_min_cells": 5,
        "solo_max_epochs": 400,
        "neighborhood_graph_n_neighbors": 15,
        "umap_min_dist": 0.5,
        "umap_spread": 1.0,
        "umap_n_components": 2,
        "python_env_version": "immune_aging.py_env.v2",
        "r_setup_version": "immune_aging.R_setup.v2"
    }
    filename = os.path.join(output_destination,"process_sample.configs.{}.txt".format(sample_id))
    with open(filename, 'w') as f:
        json.dump(sample_configs, f)
