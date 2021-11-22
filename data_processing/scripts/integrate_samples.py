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
import scanpy as sc
import numpy as np
import pandas as pd
import scvi
import hashlib
import traceback

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
from logger import SimpleLogger

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

logger = SimpleLogger(filename = logger_file_path)
logger.add_to_log("Running integrate_samples.py...")
logger.add_to_log("Starting time: {}".format(get_current_time()))
with open(integrate_samples_script, "r") as f:
    logger.add_to_log("integrate_samples.py md5 checksum: {}\n".format(hashlib.md5(bytes(f.read(), 'utf-8')).hexdigest()))

logger.add_to_log("using the following configurations:\n{}".format(str(configs)))
logger.add_to_log("Configs version: " + version)
logger.add_to_log("New configs version: " + str(is_new_version))

s3_url = "s3://immuneaging/integrated_samples/{}_level".format(configs["integration_level"])
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
        adata_dict[sample_ids[j]].obsm["protein_expression"] = df_ctrl

logger.add_to_log("Concatenating all cells from all samples...")
adata = adata_dict[sample_ids[0]]
if len(sample_ids) > 1:
    adata = adata.concatenate([adata_dict[sample_ids[j]] for j in range(1,len(sample_ids))], join="outer")

# protein QC
if len(proteins) > 0:
    logger.add_to_log("Running protein QC...")
    num_sds = 5
    # (1) remove proteins that have very low library size (across all cells; consider only cells that have non NaN values)
    # get library size per cell, normalized by num of non-nan proteins
    lib_sizes = adata.obsm["protein_expression"].sum(min_count = 1)/(adata.obsm["protein_expression"].shape[1]-adata.obsm["protein_expression"].isnull.sum())
    sd = np.nanstd(lib_sizes)
    mean = np.nanmean(lib_sizes)
    keep1 = np.logical_or(lib_sizes.isna(), np.logical_and(lib_sizes >= (mean-sd*num_sds), lib_sizes <= (mean+sd*num_sds)))
	# (2) remove cells that demonstrate extremely high values in any given protein
    sds = np.nanstd(adata.obsm["protein_expression"], axis=1)
    means = np.nanmean(adata.obsm["protein_expression"], axis=1)
    keep2 = np.logical_or(adata.obsm["protein_expression"].isnull(),
        np.logical_and(adata.obsm["protein_expression"] >= (mean-sd*num_sds), adata.obsm["protein_expression"] <= (mean+sd*num_sds)))
    keep = np.logical_and(keep1, keep2)
    adata = adata[keep,]
    logger.add_to_log("{} cells with outlier protein levels were filtered out, leaving a total of {} cells ({}% removed)".format(
        len(keep)-adata.n_obs, adata.n_obs, 100*(len(keep)-adata.n_obs)/len(keep)))

if "is_solo_singlet" in adata.obs:
    del adata.obs["is_solo_singlet"]

if "Classification" in adata.obs:
    del adata.obs["Classification"]

logger.add_to_log("A total of {} cells and {} genes are available after merge.".format(adata.n_obs, adata.n_vars))

# set the train size; this train size was justified by an experiment that is described here https://yoseflab.github.io/scvi-tools-reproducibility/runtime_analysis_reproducibility/
if 0.1 * adata.n_obs < 20000:
    train_size = 0.9
else:
    train_size = 1-(20000/adata.n_obs)

configs["train_size"] = train_size

try:
    is_cite = "protein_expression" in adata.obsm
    if is_cite:
        logger.add_to_log("Detected Antibody Capture features.")
    rna = adata.copy()
    batch_key = "batch"
    logger.add_to_log("Filtering out vdj genes...")
    filter_vdj_genes(rna, configs["vdj_genes"], data_dir, logger)
    logger.add_to_log("Detecting highly variable genes...")
    sc.pp.highly_variable_genes(rna, n_top_genes=configs["n_highly_variable_genes"], subset=True,
        flavor=configs["highly_variable_genes_flavor"], batch_key=batch_key, span = 1.0)
    # use the decontaminated counts as the input for scvi and totalvi
    rna.X = np.round(rna.layers["decontaminated_counts"])
    # scvi
    key = "X_scvi_integrated"
    _, scvi_model_file = run_model(rna, configs, batch_key, None, "scvi", prefix, version, data_dir, logger, key)
    logger.add_to_log("Calculate neighbors graph and UMAP based on scvi components...")
    neighbors_key = "scvi_integrated_neighbors"
    sc.pp.neighbors(rna, n_neighbors=configs["neighborhood_graph_n_neighbors"], use_rep=key, key_added=neighbors_key) 
    rna.obsm["X_umap_scvi_integrated"] = sc.tl.umap(rna, min_dist=configs["umap_min_dist"], spread=float(configs["umap_spread"]),
            n_components=configs["umap_n_components"], neighbors_key=neighbors_key, copy=True).obsm["X_umap"]
    if is_cite:
        # totalVI
        key = "X_totalVI_integrated"
        # if there is no protein information for some of the cells set them to zero (instead of NaN)
        rna.obsm["protein_expression"] = rna.obsm["protein_expression"].fillna(0)
        # there are known spurious failures with totalVI (such as "invalid parameter loc/scale")
        # so we try a few times then carry on with the rest of the script as we can still mine the
        # rest of the data regardless of CITE info
        retry_count = 4
        try:
            _, totalvi_model_file = run_model(rna, configs, batch_key, "protein_expression", "totalvi", prefix, version, data_dir, logger, latent_key=key, max_retry_count=retry_count)
            logger.add_to_log("Calculate neighbors graph and UMAP based on totalVI components...")
            neighbors_key = "totalvi_integrated_neighbors"
            sc.pp.neighbors(rna, n_neighbors=configs["neighborhood_graph_n_neighbors"],use_rep=key, key_added=neighbors_key) 
            rna.obsm["X_umap_totalvi_integrated"] = sc.tl.umap(rna, min_dist=configs["umap_min_dist"], spread=float(configs["umap_spread"]),
                n_components=configs["umap_n_components"], neighbors_key=neighbors_key, copy=True).obsm["X_umap"]
        except Exception as err:
            logger.add_to_log("Execution of totalVI failed with the following error (latest) with retry count {}: {}. Moving on...".format(retry_count, err), "warning")
            is_cite = False
    # pca
    logger.add_to_log("Calculating PCA...")
    sc.pp.pca(rna)
    logger.add_to_log("Calculating neighborhood graph and UMAP based on PCA...")
    key = "pca_neighbors"
    sc.pp.neighbors(rna, n_neighbors=configs["neighborhood_graph_n_neighbors"],
        use_rep="X_pca", key_added=key) 
    rna.obsm["X_umap_pca"] = sc.tl.umap(rna, min_dist=configs["umap_min_dist"], spread=float(configs["umap_spread"]),
        n_components=configs["umap_n_components"], neighbors_key=key, copy=True).obsm["X_umap"]
    # update the adata with the components of the dim reductions and umap coordinates
    adata.obsm.update(rna.obsm)
    # save the identify of the most variable genes used
    adata.var["is_highly_variable_gene"] = adata.var.index.isin(rna.var.index)
except Exception as err:
    logger.add_to_log("Execution failed with the following error: {}.\n{}".format(err, traceback.format_exc()), "critical")
    logger.add_to_log("Terminating execution prematurely.")
    if not sandbox_mode:
        # upload log to S3
        sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{}/ --exclude "*" --include {}'.format( \
            data_dir, s3_url, prefix, version, logger_file)
        os.system(sync_cmd)
    print(err)
    sys.exit()

logger.add_to_log("Saving h5ad files...")
adata.obs["age"] = adata.obs["age"].astype(str)
adata.obs["BMI"] = adata.obs["BMI"].astype(str)
adata.obs["height"] = adata.obs["height"].astype(str)
adata.write(os.path.join(data_dir,output_h5ad_file))

###############################################################
###### OUTPUT UPLOAD TO S3 - ONLY IF NOT IN SANDBOX MODE ######
###############################################################

if not sandbox_mode:
    logger.add_to_log("Uploading h5ad file to S3...")
    sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{} --exclude "*" --include {}'.format(
        data_dir, s3_url, prefix, version, output_h5ad_file)
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
    logger.add_to_log("Uploading model files (a single .zip file for each model) to S3...")
    sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{}/ --exclude "*" --include {}'.format(
        data_dir, s3_url, prefix, version, scvi_model_file)
    if is_cite:
        sync_cmd += ' --include {}'.format(totalvi_model_file)
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))    

logger.add_to_log("Execution of integrate_samples.py is complete.")

logger.add_to_log("Number of cells: {}, number of genes: {}.".format(adata.n_obs, adata.n_vars))

logging.shutdown()
if not sandbox_mode:
    # Uploading log file to S3.
    sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{} --exclude "*" --include {}'.format(
        data_dir, s3_url, prefix, version, logger_file)
    os.system(sync_cmd)
