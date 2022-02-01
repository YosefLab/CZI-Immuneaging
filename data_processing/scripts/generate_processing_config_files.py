## This script can be used to generate config files for all libraries and all samples from a given donor; currently only GEX libraries are considered (i.e. ignores BCR/TCR libs).
## Note that 
## Run as follows: python generate_processing_config_files.py <config_type> <code_path> <output_destination> <donor_id> <seq_run> <processed_lib_version>
## where config_type can be one of "library", "sample" or "all"
## and "processed_lib_version" is the version for all processed libraries that we use to process a sample (if config_type is "library" this is unused)
## Note that code_path, output_destination, and s3_access_file will be set later via generate_processing_scripts.py
## After generating the files run the following commands to upload to aws (after setting the aws credentials as env variables):
## aws s3 sync config_files s3://immuneaging/job_queue/process_library/ --exclude "*" --include "process_library*.configs.txt"
## aws s3 sync config_files s3://immuneaging/job_queue/process_sample/ --exclude "*" --include "process_sample.configs*.txt"
## assuming config_files is the directory containing the config files.
## Then, the script generate_processing_scripts.py can be used to create .sh files to execute the actual processing based on the config files.

import sys
import os
import json
from typing import List

config_type = sys.argv[1]
code_path = sys.argv[2]
output_destination = sys.argv[3]
donor_id = sys.argv[4]
seq_run = sys.argv[5]
processed_lib_version = sys.argv[6]

sys.path.append(code_path)
from utils import *

samples = read_immune_aging_sheet("Samples")
indices = (samples["Donor ID"] == donor_id) & (samples["Seq run"] == float(seq_run))

if config_type in ["library", "all"]:
    # create config files for library processing
    def add_lib(lib_type: str, all_libs: set) -> None:
        column_name = "{} lib".format(lib_type)
        libs_all = samples[indices][column_name]
        gex_libs_all = samples[indices]["GEX lib"]
        for i in range(len(libs_all)):
            libs = libs_all.iloc[i].split(",")
            gex_libs = gex_libs_all.iloc[i].split(",")
            for j in range(len(libs)):
                lib = libs[j]
                corresponding_gex_lib = gex_libs[j]
                aligned_lib_version = "v1" # TODO infer this from AWS since different lib types can have different aligned versions (github issue #51)
                all_libs.add((lib,lib_type,aligned_lib_version,corresponding_gex_lib))

    all_libs = set()
    add_lib("GEX", all_libs)
    add_lib("BCR", all_libs)
    add_lib("TCR", all_libs)

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
            "python_env_version": "immune_aging.py_env.v3",
            "r_setup_version": "immune_aging.R_setup.v2"
        }
        filename = os.path.join(output_destination,
            "process_library.{}.{}.{}.{}.configs.txt".format(donor_id,seq_run,lib_id, lib_type))
        with open(filename, 'w') as f:
            json.dump(lib_configs, f)

if config_type in ["sample", "all"]:
    # create config files for sample processing
    sample_ids = samples[indices]["Sample_ID"]
    organs = samples[indices]["Organ"]
    for i in range(len(sample_ids)):
        sample_id = sample_ids.iloc[i]
        is_jejunum = organs.iloc[i] in ["JEJ", "JEJEPI", "JEJLP"]

        def add_libs(lib_type: str, all_libs: List[str], all_lib_types: List[str], all_lib_versions: List[str]) -> None:
            lib_ids = [i for i in samples[samples["Sample_ID"] == sample_id]["{} lib".format(lib_type)].iloc[0].split(",")]
            all_libs += lib_ids
            all_lib_types += [lib_type] * len(lib_ids)
            all_lib_versions += [processed_lib_version] * len(lib_ids)

        all_libs = []
        all_lib_types = []
        all_lib_versions = []
        add_libs("GEX", all_libs, all_lib_types, all_lib_versions)
        add_libs("BCR", all_libs, all_lib_types, all_lib_versions)
        add_libs("TCR", all_libs, all_lib_types, all_lib_versions)

        sample_configs = {
            "sandbox_mode": "False",
            "data_owner": "valehvpa",
            "code_path": "./",
            "output_destination": "./",
            "s3_access_file": "./",
            "processed_libraries_dir": "",
            "donor": donor_id,
            "seq_run": seq_run,
            "sample_id": sample_id,
            "library_ids": ",".join(all_libs),
            "library_types": ",".join(all_lib_types),
            "processed_library_configs_version": ",".join(all_lib_versions),
            "min_cells_per_library": 50 if is_jejunum else 200, # jejunum samples are generally less enriched as they are less available to sequence
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
            "celltypist_model_urls": "https://celltypist.cog.sanger.ac.uk/models/v6/Immune_All_High.pkl,https://celltypist.cog.sanger.ac.uk/models/v6/Immune_All_Low.pkl",
            "rbc_model_url": "s3://immuneaging/unpublished_celltypist_models/RBC_model_CZI.pkl",
            "vdj_genes": "s3://immuneaging/vdj_genes/vdj_gene_list_v1.csv",
            "python_env_version": "immune_aging.py_env.v3",
            "r_setup_version": "immune_aging.R_setup.v2",
        }
        filename = os.path.join(output_destination,"process_sample.configs.{}.txt".format(sample_id))
        with open(filename, 'w') as f:
            json.dump(sample_configs, f)
