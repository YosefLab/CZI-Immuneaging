import sys
import logging
import os
import json
import anndata
import scanpy as sc
import numpy as np
import pandas as pd
import subprocess
import hashlib
import time

from logger import SimpleLogger
from utils import *

process_lib_script = sys.argv[0]
configs_file = sys.argv[1]

#########################################################
###### INITIALIZATIONS AND PREPARATIONS BEGIN HERE ######
#########################################################

with open(configs_file) as f: 
    data = f.read()	

configs = json.loads(data)

sys.path.append(configs["code_path"])

AUTHORIZED_EXECUTERS = ["b750bd0287811e901c88dc328187e25f"] # md5 checksums of the AWS_SECRET_ACCESS_KEY value of those that are authorized to run process_library and upload outputs to the server; not that individuals with upload permission to aws can bypass that by changing the code - this is just designed to alert users that they should only use sandbox mode.
VARIABLE_CONFIG_KEYS = ["data_owner","s3_access_file","code_path","output_destination"] # config changes only to these fields will not initialize a new configs version
sc.settings.verbosity = 3   # verbosity: errors (0), warnings (1), info (2), hints (3)

timestamp = get_current_time()
sandbox_mode = configs["sandbox_mode"] == "True"

# apply the aws credentials to allow access though aws cli; make sure the user is authorized to run in non-sandbox mode if applicable
s3_dict = set_access_keys(configs["s3_access_file"], return_dict = True)
assert hashlib.md5(bytes(s3_dict["AWS_SECRET_ACCESS_KEY"], 'utf-8')).hexdigest() in AUTHORIZED_EXECUTERS, "You are not authorized to run this script in a non sanbox mode; please set sandbox_mode to True"
set_access_keys(configs["s3_access_file"])

# create a new directory for the data and outputs
data_dir = os.path.join(configs["output_destination"], "_".join([configs["donor"], configs["seq_run"]]))
os.system("mkdir -p " + data_dir)
prefix = "{}_{}_{}_{}".format(configs["donor"], configs["seq_run"], configs["library_type"], configs["library_id"])
data_dir = os.path.join(data_dir, prefix)
os.system("mkdir -p " + data_dir)

# get the latest processed version of the library
s3_path = "s3://immuneaging/processed_libraries/" + prefix
is_new_version, version = get_configs_status(configs, s3_path, "process_library.configs." + prefix,
    VARIABLE_CONFIG_KEYS, data_dir)
output_configs_file = "process_library.configs.{}.{}.txt".format(prefix,version)

# set up logger
logger_file = os.path.join("process_library.{}.{}.log".format(prefix,version))
logger_file_path = os.path.join(data_dir, logger_file)
if os.path.isfile(logger_file_path):
	os.remove(logger_file_path)

logger = SimpleLogger(filename = logger_file_path)
logger.add_to_log("Running process_library.py...")
logger.add_to_log("Starting time: {}".format(timestamp))
with open(process_lib_script, "r") as f:
    logger.add_to_log("process_library.py md5 checksum: {}\n".format(hashlib.md5(bytes(f.read(), 'utf-8')).hexdigest()))

logger.add_to_log("using the following configurations:\n{}".format(str(configs)))
logger.add_to_log("Configs version: " + version)
logger.add_to_log("New configs version: " + str(is_new_version))

h5ad_file = "{}.processed.{}.h5ad".format(prefix, version)
if is_new_version:
    logger.add_to_log("Uploading new configs version to S3...")
    cp_cmd = "cp {} {}".format(configs_file, os.path.join(data_dir,output_configs_file))
    os.system(cp_cmd)
    sync_cmd = 'aws s3 sync {} s3://immuneaging/processed_libraries/{}/{}/ --exclude "*" --include {}'.format(
        data_dir, prefix, version, output_configs_file)
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
else:
    logger.add_to_log("Checking if h5ad file already exists on S3...")
    h5ad_file_exists = False
    logger_file_exists = False
    ls_cmd = "aws s3 ls s3://immuneaging/processed_libraries/{}/{} --recursive".format(prefix,version)
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

summary = ["\n{0}\nExecution summary\n{0}".format("="*25)]

############################################
###### LIBRARY PROCESSING BEGINS HERE ######
############################################

logger.add_to_log("Downloading h5ad file of aligned library from S3...")
aligned_h5ad_file = "{}_{}.{}.{}.h5ad".format(configs["donor"], configs["seq_run"],
    configs["library_id"], configs["aligned_library_configs_version"])
sync_cmd = 'aws s3 sync s3://immuneaging/aligned_libraries/{}/{}_{}_{}_{}/ {} --exclude "*" --include {}'.format(
    configs["aligned_library_configs_version"], configs["donor"], configs["seq_run"], configs["library_type"],
    configs["library_id"], data_dir, aligned_h5ad_file)
logger.add_to_log("sync_cmd: {}".format(sync_cmd))
logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

logger.add_to_log("Reading aligned h5ad file...")
aligned_h5ad_file = os.path.join(data_dir, aligned_h5ad_file)
adata = sc.read_h5ad(aligned_h5ad_file)

summary.append("Started with a total of {} cells and {} genes.".format(
    adata.n_obs, adata.n_vars))

logger.add_to_log("Applying basic filters...")

n_cells_before = adata.n_obs
sc.pp.filter_cells(adata, min_genes=configs["filter_cells_min_genes"])
logger.add_to_log("Filtered out {} cells that have less than {} genes expressed.".format(n_cells_before-adata.n_obs, configs["filter_cells_min_genes"]))
n_genes_before = adata.n_vars
sc.pp.filter_genes(adata, min_cells=configs["filter_genes_min_cells"])
logger.add_to_log("Filtered out {} genes that are detected in less than {} cells.".format(n_genes_before-adata.n_vars, configs["filter_genes_min_cells"]))

adata.var['mt'] = adata.var_names.str.startswith('MT-') # mitochondrial genes
adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL")) # ribosomal genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)
n_cells_before = adata.n_obs
adata = adata[adata.obs['pct_counts_mt'] <= configs["filter_cells_max_pct_counts_mt"], :]
logger.add_to_log("Filtered out {} cells with more than {}\% counts coming from mitochondrial genes.".format(n_cells_before-adata.n_obs, configs["filter_cells_max_pct_counts_mt"]))
n_cells_before = adata.n_obs
adata = adata[adata.obs['pct_counts_ribo'] >= configs["filter_cells_min_pct_counts_ribo"], :]
logger.add_to_log("Filtered out {} cells with less than {}\% counts coming from ribosomal genes.".format(n_cells_before-adata.n_obs, configs["filter_cells_min_pct_counts_ribo"]))

genes_to_exclude = np.zeros((adata.n_vars,), dtype=bool)
if configs["genes_to_exclude"] != "None":
    for gene in configs["genes_to_exclude"].split(','):
        genes_to_exclude = np.add(genes_to_exclude,adata.var_names.str.startswith(gene))

if configs["exclude_mito_genes"] == "True":
    genes_to_exclude = np.add(genes_to_exclude,adata.var_names.str.startswith('MT-'))

genes_to_exclude_names = adata.var_names[np.where(genes_to_exclude)]
genes_to_keep = np.invert(genes_to_exclude)
n_genes_before = adata.n_vars
adata = adata[:,genes_to_keep]
logger.add_to_log("Filtered out the following {} genes: {}".format(n_genes_before-adata.n_vars, ", ".join(genes_to_exclude_names)))

cell_hashing = [i for i in adata.var_names[np.where(adata.var_names.str.startswith(configs["donor"]+"-"))]]
if len(cell_hashing)>1:
    logger.add_to_log("Demultiplexing is needed; using hashsolo...")
    # add the HTO counts to .obs
    adata.obs[cell_hashing] = adata[:,cell_hashing].X.toarray()
    hashsolo_priors = [float(i) for i in configs["hashsolo_priors"].split(',')]
    sc.external.pp.hashsolo(adata, cell_hashing_columns = cell_hashing, priors = hashsolo_priors, inplace = True,
        number_of_noise_barcodes = min(len(cell_hashing)-1,configs["hashsolo_number_of_noise_barcodes"]))
    num_doublets = sum(adata.obs["Classification"] == "Doublet")
    logger.add_to_log("Removing {:.2f}% of the droplets ({} droplets out of {}) called by hashsolo as doublets...".format(100*num_doublets/adata.n_obs, num_doublets, adata.n_obs))
    adata = adata[adata.obs["Classification"] != "Doublet"]
    logger.add_to_log("Adding the library ID to the cell barcode name (will allow to distinguish between differnet cells with the same barcode when integrating differnet libraries that were used to collect the same samples)...")
    adata.obs_names = adata.obs_names + "_" + configs["library_id"]
    
logger.add_to_log("Saving h5ad file...")
adata.write(os.path.join(data_dir,h5ad_file))

if not sandbox_mode:
    logger.add_to_log("Uploading h5ad file to S3...")
    sync_cmd = 'aws s3 sync {} s3://immuneaging/processed_libraries/{}/{}/ --exclude "*" --include {}'.format(
        data_dir, prefix, version, h5ad_file)
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

logger.add_to_log("Execution of process_library.py is complete.")

summary.append("Final number of cells: {}, final number of genes: {}.".format(
    adata.n_obs, adata.n_vars))
for i in summary:
    logger.add_to_log(i)

logging.shutdown()
if not sandbox_mode:
    sync_cmd = 'aws s3 sync {} s3://immuneaging/processed_libraries/{}/{}/ --exclude "*" --include {}'.format(
        data_dir, prefix, version, logger_file)
    os.system(sync_cmd)
    #logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    #logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
