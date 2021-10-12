## For speed, this script should be executed on s130 which has a GPU.
import sys
integrate_samples_script = sys.argv[0]
configs_file = sys.argv[1]

#########################################################
###### INITIALIZATIONS AND PREPARATIONS BEGIN HERE ######
#########################################################

import logging
import os
import json
import anndata
import scanpy as sc
import numpy as np
import pandas as pd
import scvi
import subprocess
import hashlib
import time

# This does two things:
# 1. Makes the logger look good in a log file
# 2. Changes a bit how torch pins memory when copying to GPU, which allows you to more easily run models in parallel with an estimated 1-5% time hit
scvi.settings.reset_logging_handler()
scvi.settings.dl_pin_memory_gpu_training = False

with open(configs_file) as f: 
    data = f.read()	

configs = json.loads(data)
sandbox_mode = configs["sandbox_mode"] == "True"
sys.path.append(configs["code_path"])
from utils import *

timestamp = get_current_time()

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
    
AUTHORIZED_EXECUTERS = ["b750bd0287811e901c88dc328187e25f", "1c75133ab6a1fc3ed9233d3fe40b3d73"] # md5 checksums of the AWS_SECRET_ACCESS_KEY value of those that are authorized to run process_library and upload outputs to the server; note that individuals with upload permission to aws can bypass that by changing the code - this is just designed to alert users that they should only use sandbox mode.
VARIABLE_CONFIG_KEYS = ["berkeley_user","s3_access_file","code_path","output_destination"] # config changes only to these fields will not initialize a new configs version
LOGGER_LEVEL = logging.DEBUG
sc.settings.verbosity = 3   # verbosity: errors (0), warnings (1), info (2), hints (3)

# apply the aws credentials to allow access though aws cli; make sure the user is authorized to run in non-sandbox mode if applicable
s3_dict = set_access_keys(s3_access_file, return_dict = True)
assert sandbox_mode or hashlib.md5(bytes(s3_dict["AWS_SECRET_ACCESS_KEY"], 'utf-8')).hexdigest() in AUTHORIZED_EXECUTERS, "You are not authorized to run this script in a non sanbox mode; please set sandbox_mode to True"
set_access_keys(s3_access_file)

# create a new directory for the data and outputs
data_dir = os.path.join(output_destination, prefix)
os.system("mkdir -p " + data_dir)

# check for previous versions of integrated data
s3_path = "s3://immuneaging/integrated_samples/{}".format(prefix)
is_new_version, version = get_configs_status(configs, s3_path, "integrate_samples.configs."+prefix,
    VARIABLE_CONFIG_KEYS, data_dir)
output_configs_file = "integrate_samples.configs.{}.{}.txt".format(prefix,version)

# set up logger
logger_file = "integrate_samples.{}.{}.log".format(prefix,version)
logger_file_path = os.path.join(data_dir, logger_file)
if os.path.isfile(logger_file_path):
	os.remove(logger_file_path)

output_h5ad_file = "{}.{}.h5ad".format(prefix, version)
output_h5ad_model_file = "{}.{}.model_data.h5ad".format(prefix, version)

start_logger(level = LOGGER_LEVEL, filename = logger_file_path)
add_to_log("Running integrate_samples.py...")
add_to_log("Starting time: {}".format(timestamp))
with open(integrate_samples_script, "r") as f:
    add_to_log("process_sample.py md5 checksum: {}\n".format(hashlib.md5(bytes(f.read(), 'utf-8')).hexdigest()))

add_to_log("using the following configurations:\n{}".format(str(configs)))
add_to_log("Configs version: " + version)
add_to_log("New configs version: " + str(is_new_version))

h5ad_file = "{}.{}.h5ad".format(prefix, version)
if is_new_version:
    if not sandbox_mode:
        add_to_log("Uploading new configs version to S3...")
        with open(os.path.join(data_dir,output_configs_file), 'w') as f:
            json.dump(configs, f)
        sync_cmd = 'aws s3 sync {} s3://immuneaging/integrated_samples/{}/{} --exclude "*" --include {}'.format(
            data_dir, prefix, version, output_configs_file)
        add_to_log("sync_cmd: {}".format(sync_cmd))
        add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
else:
    add_to_log("Checking if h5ad file already exists on S3...")
    h5ad_file_exists = False
    logger_file_exists = False
    ls_cmd = "aws s3 ls s3://immuneaging/integrated_samples/{}/{} --recursive".format(prefix,version)
    files = os.popen(ls_cmd).read()
    add_to_log("aws response: {}\n".format(files))
    for f in files.rstrip().split('\n'):
        if f.split('/')[-1] == h5ad_file:
            h5ad_file_exists = True
        if f.split('/')[-1] == logger_file:
            logger_file_exists = True
    if h5ad_file_exists and logger_file_exists:
        add_to_log("The following h5ad file is already on S3: {}\nTerminating execution.".format(h5ad_file))
        sys.exit()
    if h5ad_file_exists and not logger_file_exists:
        add_to_log("The following h5ad file is already on S3: {}\nhowever, log file is missing on S3; proceeding with execution.".format(h5ad_file))
    if not h5ad_file_exists:
        add_to_log("The following h5ad file does not exist on S3: {}".format(h5ad_file))

h5ad_files = []
add_to_log("Downloading h5ad files of processed samples from S3...")
for j in range(len(sample_ids)):
    sample_id = sample_ids[j]
    sample_version = processed_sample_configs_version[j]
    sample_h5ad_file = "{}_GEX.processed.{}.h5ad".format(sample_id,sample_version)
    sample_h5ad_path = os.path.join(data_dir,sample_h5ad_file)
    h5ad_files.append(sample_h5ad_path)
    sync_cmd = 'aws s3 sync s3://immuneaging/processed_samples/{}_GEX/{}/ {} --exclude "*" --include {}'.format(
        sample_id,sample_version,data_dir,sample_h5ad_file)
    add_to_log("syncing {}...".format(sample_h5ad_file))
    add_to_log("sync_cmd: {}".format(sync_cmd))
    aws_response = os.popen(sync_cmd).read()
    add_to_log("aws response: {}\n".format(aws_response))
    if not os.path.exists(sample_h5ad_path):
        add_to_log("h5ad file does not exist on aws for sample {}. Terminating execution.".format(sample_id))
        sys.exit()

############################################
###### SAMPLE INTEGRATION BEGINS HERE ######
############################################

alerts = []
add_to_log("Reading h5ad files of processed samples...")
adata_dict = {}
for j in range(len(h5ad_files)):
    h5ad_file = h5ad_files[j]
    sample_id = sample_ids[j]
    adata_dict[sample_id] = sc.read_h5ad(h5ad_file)
    adata_dict[sample_id].obs["sample_id"] = sample_id

add_to_log("Concatenating all cells from all samples...".format(sample_id))
adata = adata_dict[sample_ids[0]]
if len(sample_ids) > 1:
    adata = adata.concatenate([adata_dict[sample_ids[j]] for j in range(1,len(sample_ids))], join="outer")

if "is_solo_singlet" in adata.obs:
    del adata.obs["is_solo_singlet"]

if "Classification" in adata.obs:
    del adata.obs["Classification"]

add_to_log("A total of {} cells and {} genes are available after merge.".format(adata.n_obs, adata.n_vars))

# set the train size; this train size was justified by an experiment that is described here https://yoseflab.github.io/scvi-tools-reproducibility/runtime_analysis_reproducibility/
if 0.1 * adata.n_obs < 20000:
    train_size = 0.9
else:
    train_size = 1-(20000/adata.n_obs)
configs["train_size"] = train_size

try:
    is_cite = "protein_expression" in adata.obsm ## TODO make sure that if there are only hashtag antibodies then protein_expression is not in obsm
    assert is_cite ## TODO remove after adding the case of no cite data
    if is_cite:
        add_to_log("Detected Antibody Capture features.")    
    rna = adata.copy()
    batch_key = "batch"
    add_to_log("Detecting highly variable genes...")
    sc.pp.highly_variable_genes(rna, n_top_genes=configs["n_highly_variable_genes"], subset=True,
        flavor=configs["highly_variable_genes_flavor"], batch_key=batch_key, span = 1.0)
    # use the decontaminated counts as the input for scvi and totalvi
    rna.X = np.round(rna.layers["decontaminated_counts"])
    if is_cite:
        # if there is no protein information for some of the cells set them to zero (instead of NaN)
        rna.obsm["protein_expression"] = rna.obsm["protein_expression"].fillna(0)
        key = "X_totalVI_integrated"
        rna, totalvi_model, totalvi_model_file = run_model(rna, configs, batch_key, "protein_expression", "totalvi", prefix, version, data_dir, key)
        add_to_log("Calculate neighbors graph and UMAP...")
        neighbors_key = "totalvi_integrated_neighbors"
        sc.pp.neighbors(rna, n_neighbors=configs["neighborhood_graph_n_neighbors"],
            use_rep=key, key_added=neighbors_key) 
        rna.obsm["X_umap_totalvi_integrated"] = sc.tl.umap(rna, min_dist=configs["umap_min_dist"], spread=float(configs["umap_spread"]),
            n_components=configs["umap_n_components"], neighbors_key=neighbors_key, copy=True).obsm["X_umap"]
        adata.obsm.update(rna.obsm)
except Exception as err:
    add_to_log("Execution failed with the following error:\n{}".format(err))
    add_to_log("Terminating execution prematurely.")
    if not sandbox_mode:
        # upload log to S3
        sync_cmd = 'aws s3 sync {} s3://immuneaging/integrated_samples/{}/{}/ --exclude "*" --include {}'.format( \
            data_dir, prefix, version, logger_file)
        os.system(sync_cmd)
    print(err)
    sys.exit()

add_to_log("Saving h5ad files...")
adata.write(os.path.join(data_dir,output_h5ad_file))
# TODO: save also the version of the data that was used for the model fitting (after feature selection etc.); can be useful for reference-based integration
# rna.write(os.path.join(data_dir,"{}.totalvi_model".format(prefix),output_h5ad_model_file))

###############################################################
###### OUTPUT UPLOAD TO S3 - ONLY IF NOT IN SANDBOX MODE ######
###############################################################

if not sandbox_mode:
    add_to_log("Uploading h5ad file to S3...")
    sync_cmd = 'aws s3 sync {} s3://immuneaging/integrated_samples/{}/{} --exclude "*" --include {}'.format(
        data_dir, prefix, version, output_h5ad_file)
    add_to_log("sync_cmd: {}".format(sync_cmd))
    add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

add_to_log("Execution of integrate_samples.py is complete.")

add_to_log("Number of cells: {}, number of genes: {}.".format(adata.n_obs, adata.n_vars))

logging.shutdown()
if not sandbox_mode:
    # Uploading log file to S3.
    sync_cmd = 'aws s3 sync {} s3://immuneaging/integrated_samples/{}/{} --exclude "*" --include {}'.format(
        data_dir, prefix, version, logger_file)
    os.system(sync_cmd)
