import sys
process_sample_script = sys.argv[0]
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
import zipfile

with open(configs_file) as f: 
    data = f.read()	

configs = json.loads(data)
sandbox_mode = configs["sandbox_mode"] == "True"
output_destination = configs["output_destination"]
donor = configs["donor"]
seq_run = configs["seq_run"]
sample_id = configs["sample_id"]
library_type = configs["library_type"]

sys.path.append(configs["code_path"])
from utils import *

AUTHORIZED_EXECUTERS = ["b750bd0287811e901c88dc328187e25f"] # md5 checksums of the AWS_SECRET_ACCESS_KEY value of those that are authorized to run process_library and upload outputs to the server; not that individuals with upload permission to aws can bypass that by changing the code - this is just designed to alert users that they should only use sandbox mode.
VARIABLE_CONFIG_KEYS = ["data_owner","s3_access_file","code_path","output_destination"] # config changes only to these fields will not initialize a new configs version
LOGGER_LEVEL = logging.DEBUG
sc.settings.verbosity = 3   # verbosity: errors (0), warnings (1), info (2), hints (3)

# a map between fields in the Donors sheet of the Google Spreadsheet to metadata fields
DONORS_FIELDS = {"Donor ID": "donor_id",
    "Site (UK/ NY)": "site",
    "DCD/DBD": "DCD/DBD",
    "Age (years)": "age",
    "Sex": "sex",
    "ethnicity/race": "ethnicity/race",
    "cause of death": "death_cause",
    "mech of injury": "mech_injury",
    "height (cm)": "height",
    "BMI (kg/m^2)": "BMI",
    "lipase level": "lipase_level",
    "blood sugar (mg/dL)": "blood_sugar",
    "Period of time in relation to smoking": "smoking_period",
    "smoker (pack-years)": "smoking_packs_year",
    "EBV status": "EBV",
    "CMV status": "CMV"}

# a map between fields in the Samples sheet of the Google Spreadsheet to metadata fields
SAMPLES_FIELDS = {"Sample_ID": "sample_id",
    "Seq run": "seq_run",
    "Fresh/frozen": "Fresh/frozen",
    "Cell type": "sample_cell_type",
    "Sorting": "sorting",
    "Stimulation": "stimulation"}

timestamp = get_current_time()

# apply the aws credentials to allow access though aws cli; make sure the user is authorized to run in non-sandbox mode if applicable
s3_dict = set_access_keys(configs["s3_access_file"], return_dict = True)
assert hashlib.md5(bytes(s3_dict["AWS_SECRET_ACCESS_KEY"], 'utf-8')).hexdigest() in AUTHORIZED_EXECUTERS, "You are not authorized to run this script in a non sanbox mode; please set sandbox_mode to True"
set_access_keys(configs["s3_access_file"])

# create a new directory for the data and outputs
data_dir = os.path.join(output_destination, "_".join([donor, seq_run]))
os.system("mkdir -p " + data_dir)

# check for previous versions of the processed sample
prefix = "{}_{}".format(sample_id, library_type)
s3_path = "s3://immuneaging/processed_samples/" + prefix
is_new_version, version = get_configs_status(configs, s3_path, "process_sample.configs." + prefix,
    VARIABLE_CONFIG_KEYS, data_dir)
output_configs_file = "process_sample.configs.{}.{}.txt".format(prefix,version)

# set up logger
logger_file = "process_sample.{}.{}.log".format(prefix,version)
logger_file_path = os.path.join(data_dir, logger_file)
if os.path.isfile(logger_file_path):
	os.remove(logger_file_path)

start_logger(level = LOGGER_LEVEL, filename = logger_file_path)
add_to_log("Running process_sample.py...")
add_to_log("Starting time: {}".format(timestamp))
with open(process_sample_script, "r") as f:
    add_to_log("process_sample.py md5 checksum: {}\n".format(hashlib.md5(bytes(f.read(), 'utf-8')).hexdigest()))

add_to_log("using the following configurations:\n{}".format(str(configs)))
add_to_log("Configs version: " + version)
add_to_log("New configs version: " + str(is_new_version))

h5ad_file = "{}.processed.{}.h5ad".format(prefix, version)
if is_new_version:
    add_to_log("Uploading new configs version to S3...")
    cp_cmd = "cp {} {}".format(configs_file, os.path.join(data_dir,output_configs_file))
    os.system(cp_cmd)
    sync_cmd = 'aws s3 sync {} s3://immuneaging/processed_samples/{}/{}/ --exclude "*" --include {}'.format(
        data_dir, prefix, version, output_configs_file)
    add_to_log("sync_cmd: {}".format(sync_cmd))
    add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
else:
    add_to_log("Checking if h5ad file already exists on S3...")
    h5ad_file_exists = False
    logger_file_exists = False
    ls_cmd = "aws s3 ls s3://immuneaging/processed_samples/{}/{} --recursive".format(prefix,version)
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

library_versions = configs["processed_library_configs_version"].split(',')
library_ids = configs["library_ids"].split(',')
if not sandbox_mode:
    add_to_log("Downloading h5ad files of processed libraries from S3...")
    for j in range(len(library_ids)):
        library_id = library_ids[j]
        library_version = library_versions[j]
        lib_h5ad_file = "{}_{}_{}_{}.processed.{}.h5ad".format(donor, seq_run,
            library_type, library_id, library_version)    
        sync_cmd = 'aws s3 sync s3://immuneaging/processed_libraries/{}_{}_{}_{}/{}/ {} --exclude "*" --include {}'.format(
            donor, seq_run, library_type, library_id, library_version, data_dir, lib_h5ad_file)
        add_to_log("syncing {}...".format(lib_h5ad_file))
        add_to_log("sync_cmd: {}".format(sync_cmd))
        add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

############################################
###### SAMPLE PROCESSING BEGINS HERE #######
############################################

add_to_log("Reading h5ad files of processed libraries...")
adata_dict = {}
for j in range(len(library_ids)):
    library_id = library_ids[j]
    library_version = library_versions[j]
    lib_h5ad_file = os.path.join(data_dir, "{}_{}_{}_{}.processed.{}.h5ad".format(donor, seq_run,
        library_type, library_id, library_version))
    adata_dict[library_id] = sc.read_h5ad(lib_h5ad_file)
    adata_dict[library_id].obs["library_id"] = library_id

add_to_log("Concatenating all cells of sample {} from all libraries...".format(sample_id))
adata = adata_dict[library_ids[0]][adata_dict[library_ids[0]].obs["Classification"] == sample_id]
adata = adata.concatenate([adata_dict[library_ids[j]][adata_dict[library_ids[j]].obs["Classification"] == sample_id] for j in range(1,len(library_ids))])

add_to_log("Adding metadata...")

add_to_log("Downloading the Donors sheet from the Google spreadsheets...")
donors = read_immune_aging_sheet("Donors")
donor_index = donors["Donor ID"] == donor
add_to_log("Adding donor-level metadata...")
for k in DONORS_FIELDS.keys():
    adata.obs[DONORS_FIELDS[k]] = donors[k][donor_index].values[0]

add_to_log("Downloading the Samples sheet from the Google spreadsheets...")
samples = read_immune_aging_sheet("Samples")
add_to_log("Adding sample-level metadata...")
sample_index = samples["Sample_ID"] == sample_id
for k in SAMPLES_FIELDS.keys():
    adata.obs[SAMPLES_FIELDS[k]] = samples[k][sample_index].values[0]

if library_type == "GEX":
    adata.obs["GEX_chem"] = samples["GEX chem"][sample_index].values[0]
    adata.obs["HTO_chem"] = samples["HTO chem"][sample_index].values[0]
    adata.obs["ADT_chem"] = samples["CITE chem"][sample_index].values[0]

is_hto = len(library_ids)>1
is_cite = "feature_types" in adata.var and "Antibody Capture" in np.unique(adata.var.feature_types)
if is_cite:
    add_to_log("Detected Antibody Capture features.")    
    add_to_log("Normalizing CITE data separately...")
    rna = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()    
else:
    rna = adata.copy()

add_to_log("Detecting highly variable genes...")
sc.pp.highly_variable_genes(rna, n_top_genes=configs["n_highly_variable_genes"], subset=True,
    flavor=configs["highly_variable_genes_flavor"], batch_key="batch")

add_to_log("Setting up scvi...")
scvi.data.setup_anndata(rna, batch_key="batch")

add_to_log("Training scvi model...")
model = scvi.model.SCVI(rna)
model.train(check_val_every_n_epoch=1000, max_epochs=configs["scvi_max_epochs"])

add_to_log("Saving scvi latent representation...")
latent = model.get_latent_representation()
rna.obsm["X_scVI"] = latent

scvi_model_file = "{}.processed.{}.scvi_model.zip".format(prefix, version)
add_to_log("Saving the scvi model into {}...".format(scvi_model_file))
scvi_model_file_path = os.path.join(data_dir, scvi_model_file)
scvi_model_dir_path = os.path.join(data_dir,"{}.scvi_model/".format(prefix))
if os.path.isdir(scvi_model_dir_path):
    os.system("rm -r " + scvi_model_dir_path)

model.save(scvi_model_dir_path)
zipf = zipfile.ZipFile(scvi_model_file_path, 'w', zipfile.ZIP_DEFLATED)
zipdir(scvi_model_dir_path, zipf)
zipf.close()

add_to_log("Running solo for detecting doublets...")
if is_hto:
    batches = np.unique(rna.obs["batch"])
    add_to_log("Running solo on the following batches separately: {}".format(batches))
    is_solo_singlet = np.ones((rna.n_obs,), dtype=bool)
    for batch in batches:
        add_to_log("Running solo on batch {}...".format(batch))
        solo_batch = scvi.external.SOLO.from_scvi_model(model, restrict_to_batch=batch)
        solo_batch.train(max_epochs=configs["solo_max_epochs"])
        is_solo_singlet[(rna.obs["batch"] == batch).values] = solo_batch.predict(soft=False) == "singlet"
    rna.obs["is_solo_singlet"] = is_solo_singlet
else:
    add_to_log("Running solo...")
    solo = scvi.external.SOLO.from_scvi_model(model)
    solo_batch.train(max_epochs=configs["solo_max_epochs"])
    is_solo_singlet = solo.predict(soft=False) == "singlet"

add_to_log("Removing doublets...")
n_obs_before = rna.n_obs
rna = rna[is_solo_singlet,]
add_to_log("Removed {} estimated doublets; {} droplets remained.".format(n_obs_before-rna.n_obs,rna.n_obs))

add_to_log("Normalizing RNA...")
rna.layers["counts"] = rna.X.copy()
sc.pp.normalize_total(rna, target_sum=configs["normalize_total_target_sum"])
sc.pp.log1p(rna)
sc.pp.scale(rna)
rna.raw = rna

add_to_log("Calculating PCA...")
sc.pp.pca(rna)

add_to_log("Calculating neighborhood graph and UMAP based on PCA...")
key = "pca_neighbors"
rna.neighbors = sc.pp.neighbors(rna, n_neighbors=configs["neighborhood_graph_n_neighbors"],
    use_rep="X_pca", key_added=key) 
rna.obsm["X_umap_pca"] = sc.tl.umap(rna, min_dist=configs["umap_min_dist"], spread=float(configs["umap_spread"]),
    n_components=configs["umap_n_components"], neighbors_key=key, copy=True).obsm["X_umap"]

add_to_log("Calculating neighborhood graph and UMAP based on SCVI components...")
key = "scvi_neighbors"
rna.neighbors = sc.pp.neighbors(rna, n_neighbors=configs["neighborhood_graph_n_neighbors"],
    use_rep="X_scVI", key_added=key) 
rna.obsm["X_umap_scvi"] = sc.tl.umap(rna, min_dist=configs["umap_min_dist"], spread=float(configs["umap_spread"]),
    n_components=configs["umap_n_components"], neighbors_key=key, copy=True).obsm["X_umap"]


add_to_log("Gathering data...")
adata = adata[is_solo_singlet,]
adata.obsm.update(rna.obsm)
adata.obs[rna.obs.columns] = rna.obs
# save raw counts
adata.layers["counts"] = adata.X.copy()
adata.raw = adata
# normalize RNA
if is_cite:
    # normalize RNA and protein separately
    rna = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
    sc.pp.normalize_total(rna, target_sum=configs["normalize_total_target_sum"])
    sc.pp.log1p(rna)
    protein = adata[:, adata.var["feature_types"] == "Antibody Capture"].copy()
    sc.pp.normalize_total(protein, target_sum=configs["normalize_total_target_sum"])
    sc.pp.log1p(protein)
else:
    rna = adata.copy()

adata.X[:, (adata.var["feature_types"] == "Gene Expression").values] = rna.X
if is_cite:
    adata.X[:, (adata.var["feature_types"] == "Antibody Capture").values] = protein.X

add_to_log("Saving h5ad file...")
adata.write(os.path.join(data_dir,h5ad_file))


###############################################################
###### OUTPUT UPLOAD TO S3 - ONLY IF NOT IN SANDBOX MODE ######
###############################################################

if not sandbox_mode:
    add_to_log("Uploading h5ad file to S3...")
    sync_cmd = 'aws s3 sync {} s3://immuneaging/processed_libraries/{}/{}/ --exclude "*" --include {}'.format(
        data_dir, prefix, version, h5ad_file)
    add_to_log("sync_cmd: {}".format(sync_cmd))
    add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
    add_to_log("Uploading scvi model files (a single .zip file) to S3...")
    sync_cmd = 'aws s3 sync {} s3://immuneaging/processed_libraries/{}/{}/ --exclude "*" --include {}'.format(
        data_dir, prefix, version, scvi_model_file)
    add_to_log("sync_cmd: {}".format(sync_cmd))
    add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
    add_to_log("Execution of process_sample.py is complete.")
    # Uploading log file to S3...
    sync_cmd = 'aws s3 sync {} s3://immuneaging/processed_libraries/{}/{}/ --exclude "*" --include {}'.format(
        data_dir, prefix, version, logger_file)
    add_to_log("sync_cmd: {}".format(sync_cmd))
    add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
