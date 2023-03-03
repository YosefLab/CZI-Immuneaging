## This script can be used to generate config files for all libraries and all samples from a given donor; currently only GEX libraries are considered (i.e. ignores BCR/TCR libs).
## Note that 
## Run as follows: python generate_processing_config_files.py <config_type> <code_path> <output_destination> <donor_id> <seq_run> <s3_access_file> <processed_gex_lib_version> <processed_ir_lib_version>
## where config_type can be one of "library", "sample" or "all"
## processed_._lib_version args are only required if config_type is "sample", otherwise you can pass any string 
## Note that code_path, output_destination, and s3_access_file will be set later via generate_processing_scripts.py
## After generating the files run the following commands to upload to aws (after setting the aws credentials as env variables):
## aws s3 sync config_files s3://immuneaging/job_queue/process_library/ --exclude "*" --include "process_library*.configs.txt"
## aws s3 sync config_files s3://immuneaging/job_queue/process_sample/ --exclude "*" --include "process_sample.configs*.txt"
## assuming config_files is the directory containing the config files.
## Then, the script generate_processing_scripts.py can be used to create .sh files to execute the actual processing based on the config files.

import re
import sys
import os
import json
from typing import List
from logger import RichLogger

config_type = sys.argv[1]
code_path = sys.argv[2]
output_destination = sys.argv[3]
donor_id = sys.argv[4]
seq_run = sys.argv[5]
s3_access_file = sys.argv[6]
processed_gex_lib_version = sys.argv[7]
processed_ir_lib_version = sys.argv[8]
sandbox_mode = sys.argv[9]

sys.path.append(code_path)
from utils import *

celltypist_model_urls = "https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_Low.pkl,https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_High.pkl"
rbc_model_url = "s3://immuneaging/unpublished_celltypist_models/RBC_model_CZI.pkl"

samples = read_immune_aging_sheet("Samples")
indices = (samples["Donor ID"] == donor_id) & (samples["Seq run"] == float(seq_run))

if config_type in ["library", "all"]:
    # create config files for library processing
    def add_lib(lib_type: str, all_libs: set) -> None:
        if lib_type not in ["GEX", "BCR", "TCR"]:
            raise ValueError("Unsupported lib_type: {}. Must be one of: GEX, BCR, TCR".format(lib_type))
        column_name = "{} lib".format(lib_type)
        libs_all = samples[indices][column_name]
        gex_libs_all = samples[indices]["GEX lib"]
        failed_libs = set()
        #Read from AWS, single call and fetching is faster.
        set_access_keys(s3_access_file)
        ls_cmd = "aws s3 ls s3://immuneaging/aligned_libraries --recursive"
        ls  = os.popen(ls_cmd).read()
        
        for i in range(len(libs_all)):
            if libs_all.iloc[i] is np.nan:
                continue
            libs = libs_all.iloc[i].split(",")
            gex_libs = gex_libs_all.iloc[i].split(",")
            for j in range(len(libs)):
                lib = libs[j]
                corresponding_gex_lib = gex_libs[j]
                # find the latest aligned_lib_version
                latest_version = -1
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
            "sandbox_mode": sandbox_mode,
            "data_owner": "cane11",
            "code_path": code_path,
            "output_destination": output_destination,
            "s3_access_file": s3_access_file,
            "donor": donor_id,
            "seq_run": seq_run,
            "library_type": lib_type,
            "library_id": lib_id,
            "corresponding_gex_lib": corresponding_gex_lib,
            "filter_cells_min_genes": 400,
            "filter_cells_min_umi": 800,
            "filter_genes_min_cells": 0,
            "filter_cells_max_pct_counts_mt": 20,
            "filter_cells_min_pct_counts_ribo": 0,
            "genes_to_exclude": "MALAT1",
            "exclude_mito_genes": "True",
            "hashsolo_priors": "0.05,0.7,0.25",
            "hashsolo_number_of_noise_barcodes": None,
            "aligned_library_configs_version": aligned_lib_version,
            "python_env_version": "immune_aging.py_env.v4",
            "pipeline_version": "qc_230227",
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
            comma_sep_libs = samples[samples["Sample_ID"] == sample_id]["{} lib".format(lib_type)].iloc[0]
            if comma_sep_libs is np.nan:
                return
            lib_ids = [i for i in comma_sep_libs.split(",")]
            all_libs += lib_ids
            all_lib_types += [lib_type] * len(lib_ids)
            processed_lib_version = processed_gex_lib_version if lib_type == "GEX" else processed_ir_lib_version
            all_lib_versions += [processed_lib_version] * len(lib_ids)

        all_libs = []
        all_lib_types = []
        all_lib_versions = []
        add_libs("GEX", all_libs, all_lib_types, all_lib_versions)
        add_libs("BCR", all_libs, all_lib_types, all_lib_versions)
        add_libs("TCR", all_libs, all_lib_types, all_lib_versions)

        sample_configs = {
            "sandbox_mode": sandbox_mode,
            "data_owner": "cane11",
            "code_path": code_path,
            "output_destination": output_destination,
            "s3_access_file": s3_access_file,
            "donor": donor_id,
            "seq_run": seq_run,
            "processed_libraries_dir": "",
            "donor": donor_id,
            "seq_run": seq_run,
            "sample_id": sample_id,
            "library_ids": ",".join(all_libs),
            "library_types": ",".join(all_lib_types),
            "processed_library_configs_version": ",".join(all_lib_versions),
            "min_cells_per_library": 50, # jejunum samples are generally less enriched as they are less available to sequence, otherwise completely removing skin sample.
            "filter_decontaminated_cells_min_genes": 30,
            "normalize_total_target_sum": 10000,
            "n_highly_variable_genes": 3000,
            "gene_likelihood": "nb",
            "highly_variable_genes_flavor": "seurat_v3",
            "scvi_max_epochs": None,
            "totalvi_max_epochs": None,
            "early_stopping": True,
            "reduce_lr_on_plateau": False,
            "n_epochs_kl_warmup": 30,
            "reduce_lr_on_plateau": False,
            "empirical_protein_background_prior": "False",
            "solo_filter_genes_min_cells": 30,
            "neighborhood_graph_n_neighbors": 15,
            "umap_min_dist": 0.5,
            "umap_spread": 1.0,
            "umap_n_components": 2,
            "celltypist_model_urls": celltypist_model_urls,
            "rbc_model_url": rbc_model_url,
            "vdj_genes": "s3://immuneaging/vdj_genes/vdj_gene_list_v1.csv",
            "python_env_version": "immune_aging.py_env.v4",
            "rscript": "/home/eecs/cergen/anaconda3/envs/new_decontx/bin/Rscript",
            "pipeline_version": "qc_230227",
            "percolation_score": {
                "doublet_probability" : {"score_key": "doublet_probability"},
                "doublet_hypothesis_probability": {"score_key": "doublet_hypothesis_probability"},
                "pct_counts_hb": {"score_key": "pct_counts_hb"},
                "double_ir": {"score_key": "double_ir", "threshold": "True"},
                "celltypist_predicted_labels.RBC_model_CZI": {"score_key": "celltypist_predicted_labels.RBC_model_CZI", "threshold": "RBC"},
                "total_counts" : {"score_key": "total_counts", "threshold": 2000, "exclude_high": False},
                "n_genes" : {"score_key": "n_genes", "threshold": 1200, "exclude_high": False},
                "n_proteins" : {"score_key": "n_proteins", "threshold": 200, "exclude_high": False}
            }
        }
        filename = os.path.join(output_destination,"process_sample.configs.{}.txt".format(sample_id))
        with open(filename, 'w') as f:
            json.dump(sample_configs, f)
