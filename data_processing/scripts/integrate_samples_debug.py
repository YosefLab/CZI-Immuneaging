
#########################################################
###### INITIALIZATIONS AND PREPARATIONS BEGIN HERE ######
#########################################################

import logging
import os
import json
import scanpy as sc
import numpy as np
import pandas as pd
import scvi
import hashlib
import traceback
import sys

from utils import *
from logger import SimpleLogger

logging.getLogger('numba').setLevel(logging.WARNING)

# This does two things:
# 1. Makes the logger look good in a log file
# 2. Changes a bit how torch pins memory when copying to GPU, which allows you to more easily run models in parallel with an estimated 1-5% time hit
scvi.settings.reset_logging_handler()
scvi.settings.dl_pin_memory_gpu_training = False

def do_all(configs_file: str):
    with open(configs_file) as f: 
        data = f.read()	

    configs = json.loads(data)
    sandbox_mode = configs["sandbox_mode"] == "True"
    sys.path.append(configs["code_path"])

    output_destination = configs["output_destination"]
    prefix = configs["output_prefix"]
    s3_access_file = configs["s3_access_file"]

    # sort the sample ids lexicographically (and accordingly the versions) in order to avoid generating a new version if the exact same set of samples were previously used but in a differnet order
    order = np.argsort(configs["sample_ids"].split(","))
    sample_ids = np.array(configs["sample_ids"].split(","))[order]
    processed_sample_configs_version = np.array(configs["processed_sample_configs_version"].split(","))[order]
    # the followings are required because we check if the integration of the requested samples already exist on aws
    configs["sample_ids"] = ",".join(sample_ids)
    configs["processed_sample_configs_version"] = ",".join(processed_sample_configs_version)

    assert(len(sample_ids)>1)

    VARIABLE_CONFIG_KEYS = ["data_owner","s3_access_file","code_path","output_destination"] # config changes only to these fields will not initialize a new configs version
    sc.settings.verbosity = 3   # verbosity: errors (0), warnings (1), info (2), hints (3)

    # apply the aws credentials to allow access though aws cli; make sure the user is authorized to run in non-sandbox mode if applicable
    s3_dict = set_access_keys(s3_access_file, return_dict = True)
    assert sandbox_mode or hashlib.md5(bytes(s3_dict["AWS_SECRET_ACCESS_KEY"], 'utf-8')).hexdigest() in AUTHORIZED_EXECUTERS, "You are not authorized to run this script in a non sanbox mode; please set sandbox_mode to True"
    set_access_keys(s3_access_file)

    # create a new directory for the data and outputs
    data_dir = os.path.join(output_destination, prefix)
    os.system("mkdir -p " + data_dir)

    # check for previous versions of integrated data
    s3_url = "s3://immuneaging/integrated_samples/{}_level".format(configs["integration_level"])
    is_new_version, version = get_configs_status(configs, s3_url + "/" + prefix, "integrate_samples.configs."+prefix,
        VARIABLE_CONFIG_KEYS, data_dir)
    output_configs_file = "integrate_samples.configs.{}.{}.txt".format(prefix,version)

    # set up logger
    logger_file = "integrate_samples.{}.{}.log".format(prefix,version)
    logger_file_path = os.path.join(data_dir, logger_file)
    if os.path.isfile(logger_file_path):
        os.remove(logger_file_path)

    output_h5ad_file = "{}.{}.h5ad".format(prefix, version)
    output_h5ad_model_file = "{}.{}.model_data.h5ad".format(prefix, version)

    logger = SimpleLogger(filename = logger_file_path)
    logger.add_to_log("Running integrate_samples.py...")
    logger.add_to_log("Starting time: {}".format(get_current_time()))

    logger.add_to_log("using the following configurations:\n{}".format(str(configs)))
    logger.add_to_log("Configs version: " + version)
    logger.add_to_log("New configs version: " + str(is_new_version))

    h5ad_file = "{}.{}.h5ad".format(prefix, version)
    if is_new_version:
        if not sandbox_mode:
            logger.add_to_log("Uploading new configs version to S3...")
            with open(os.path.join(data_dir,output_configs_file), 'w') as f:
                json.dump(configs, f)
            sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{} --exclude "*" --include {}'.format(
                data_dir, s3_url, prefix, version, output_configs_file)
            logger.add_to_log("sync_cmd: {}".format(sync_cmd))
            logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
    else:
        logger.add_to_log("Checking if h5ad file already exists on S3...")
        h5ad_file_exists = False
        logger_file_exists = False
        ls_cmd = "aws s3 ls {}/{}/{} --recursive".format(s3_url, prefix,version)
        files = os.popen(ls_cmd).read()
        logger.add_to_log("aws response: {}\n".format(files))
        for f in files.rstrip().split('\n'):
            if f.split('/')[-1] == h5ad_file:
                h5ad_file_exists = True
            if f.split('/')[-1] == logger_file:
                logger_file_exists = True
        if h5ad_file_exists and logger_file_exists:
            logger.add_to_log("The following h5ad file is already on S3: {}\nTerminating execution.".format(h5ad_file))
            sys.exit()
        if h5ad_file_exists and not logger_file_exists:
            logger.add_to_log("The following h5ad file is already on S3: {}\nhowever, log file is missing on S3; proceeding with execution.".format(h5ad_file))
        if not h5ad_file_exists:
            logger.add_to_log("The following h5ad file does not exist on S3: {}".format(h5ad_file))

    h5ad_files = []
    logger.add_to_log("Downloading h5ad files of processed samples from S3...")
    for j in range(len(sample_ids)):
        sample_id = sample_ids[j]
        sample_version = processed_sample_configs_version[j]
        sample_h5ad_file = "{}_GEX.processed.{}.h5ad".format(sample_id,sample_version)
        sample_h5ad_path = os.path.join(data_dir,sample_h5ad_file)
        h5ad_files.append(sample_h5ad_path)
        sync_cmd = 'aws s3 sync --no-progress s3://immuneaging/processed_samples/{}_GEX/{}/ {} --exclude "*" --include {}'.format(
            sample_id,sample_version,data_dir,sample_h5ad_file)
        logger.add_to_log("syncing {}...".format(sample_h5ad_file))
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        aws_response = os.popen(sync_cmd).read()
        logger.add_to_log("aws response: {}\n".format(aws_response))
        if not os.path.exists(sample_h5ad_path):
            logger.add_to_log("h5ad file does not exist on aws for sample {}. Terminating execution.".format(sample_id))
            sys.exit()

    ############################################
    ###### SAMPLE INTEGRATION BEGINS HERE ######
    ############################################

    alerts = []
    logger.add_to_log("Reading h5ad files of processed samples...")
    adata_dict = {}
    for j in range(len(h5ad_files)):
        h5ad_file = h5ad_files[j]
        sample_id = sample_ids[j]
        adata_dict[sample_id] = sc.read_h5ad(h5ad_file)
        adata_dict[sample_id].obs["sample_id"] = sample_id
        for field in ('X_pca', 'X_scVI', 'X_totalVI', 'X_umap_pca', 'X_umap_scvi', 'X_umap_totalvi'):
            if field in adata_dict[sample_id].obsm:
                del adata_dict[sample_id].obsm[field]

    # get the names of all proteins and control proteins across all samples (in case the samples were collected using more than one protein panel)
    proteins = set()
    for j in range(len(sample_ids)):
        if "protein_expression" in adata_dict[sample_ids[j]].obsm:
            proteins.update(adata_dict[sample_ids[j]].obsm["protein_expression"].columns)

    proteins_ctrl = set()
    for j in range(len(sample_ids)):
        if "protein_expression_Ctrl" in adata_dict[sample_ids[j]].obsm:
            proteins_ctrl.update(adata_dict[sample_ids[j]].obsm["protein_expression_Ctrl"].columns)

    # each sample should include all proteins; missing values are set as NaN
    if len(proteins) > 0 or len(proteins_ctrl) > 0:
        proteins = [i for i in proteins]
        proteins_ctrl = [i for i in proteins_ctrl]
        for j in range(len(sample_ids)):
            df = pd.DataFrame(columns = proteins, index = adata_dict[sample_ids[j]].obs.index)
            if "protein_expression" in adata_dict[sample_ids[j]].obsm:
                df[adata_dict[sample_ids[j]].obsm["protein_expression"].columns] = adata_dict[sample_ids[j]].obsm["protein_expression"].copy()
            adata_dict[sample_ids[j]].obsm["protein_expression"] = df
            df_ctrl = pd.DataFrame(columns = proteins_ctrl, index = adata_dict[sample_ids[j]].obs.index)
            if "protein_expression_Ctrl" in adata_dict[sample_ids[j]].obsm:
                df_ctrl[adata_dict[sample_ids[j]].obsm["protein_expression_Ctrl"].columns] = adata_dict[sample_ids[j]].obsm["protein_expression_Ctrl"].copy()
            adata_dict[sample_ids[j]].obsm["protein_expression_Ctrl"] = df_ctrl

    logger.add_to_log("Concatenating all cells from all samples...")
    adata = adata_dict[sample_ids[0]]
    if len(sample_ids) > 1:
        adata = adata.concatenate([adata_dict[sample_ids[j]] for j in range(1,len(sample_ids))], join="outer")

    # protein QC
    if "protein_levels_max_sds" in configs:
        protein_levels_max_sds = configs["protein_levels_max_sds"]
    else:
        protein_levels_max_sds = None

    if (len(proteins) > 0) and protein_levels_max_sds is not None:
        logger.add_to_log("Running protein QC...")
        n_cells_before = adata.obsm["protein_expression"].shape[0]
        n_proteins_before = adata.obsm["protein_expression"].shape[1]
        # a function for excluding cells or proteins based on extreme values (does not excludes cells that have only missing values)
        def exclude_outliers(x, num_sds):
            # a value is considered as an outlier if it is more extreme that the mean plus (or minus) num_sds times the standard deviation
            is_outlier = np.logical_and(x >= (np.nanmean(x)-np.nanstd(x)*num_sds), x <= (np.nanmean(x)+np.nanstd(x)*num_sds))
            # do not consider missing values as outliers
            return np.logical_or(x.isna(), is_outlier)
        # (1) remove cells that demonstrate extremely high (or low) protein library size; normalize library size by the total number of non NaN proteins available for the cell.
        normalized_cell_lib_size = adata.obsm["protein_expression"].sum(axis=1)/(adata.obsm["protein_expression"].shape[1]-adata.obsm["protein_expression"].isnull().sum(axis=1))
        keep = exclude_outliers(normalized_cell_lib_size, protein_levels_max_sds)
        logger.add_to_log("Removing {} cells with extreme protein values (normalized library size more extreme than {} standard deviations)...".format(np.sum(~keep), protein_levels_max_sds))
        adata = adata[keep,].copy()
        # (2) remove proteins that have extreme library size (across all cells; consider only cells that have non missing values and normalize by the number of cells with no missing values for the protein)
        normalized_protein_lib_sizes = adata.obsm["protein_expression"].sum()/(adata.obsm["protein_expression"].shape[0]-adata.obsm["protein_expression"].isnull().sum(axis=0))
        keep = exclude_outliers(normalized_protein_lib_sizes, protein_levels_max_sds)
        logger.add_to_log("Removing {} proteins with extreme total number of reads across cells (normalized levels - by the number of cells with no missing values for the protein - more extreme than {} standard deviations)...".format(
            np.sum(~keep), protein_levels_max_sds))
        adata.obsm["protein_expression"] = adata.obsm["protein_expression"][adata.obsm["protein_expression"].columns[keep.values]]
        # (3) remove proteins that were left with zero reads
        non_zero_proteins = adata.obsm["protein_expression"].sum() > 0
        num_zero_proteins = np.sum(~non_zero_proteins)
        if num_zero_proteins:
            logger.add_to_log("Removing {} proteins with zero reads across all cells...".format(num_zero_proteins))
        adata.obsm["protein_expression"] = adata.obsm["protein_expression"][adata.obsm["protein_expression"].columns[non_zero_proteins.values]]
        # summarize
        logger.add_to_log("Protein QC summary: a total of {} cells ({}%) and {} proteins ({}%) were filtered out owing to extreme values.".format(
            n_cells_before-adata.n_obs, round(100*(n_cells_before-adata.n_obs)/n_cells_before,2),
            n_proteins_before-adata.obsm["protein_expression"].shape[1], round(100*(n_proteins_before-adata.obsm["protein_expression"].shape[1])/n_proteins_before,2)))

    if "is_solo_singlet" in adata.obs:
        del adata.obs["is_solo_singlet"]

    if "Classification" in adata.obs:
        del adata.obs["Classification"]

    logger.add_to_log("A total of {} cells and {} genes are available after merge.".format(adata.n_obs, adata.n_vars))

    adata.obs["age"] = adata.obs["age"].astype(str)
    adata.obs["BMI"] = adata.obs["BMI"].astype(str)
    adata.obs["height"] = adata.obs["height"].astype(str)

    return adata