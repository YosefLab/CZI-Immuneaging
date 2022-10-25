## For speed, this script should be executed on s130 which has a GPU and high RAM
import shutil
import sys
integrate_using_scanvi_script = sys.argv[0]
configs_file = sys.argv[1]

#########################################################
###### INITIALIZATIONS AND PREPARATIONS BEGIN HERE ######
#########################################################

import logging
import os
import json
import scanpy as sc
import pandas as pd
import scvi
import hashlib
import traceback
import celltypist
import urllib.request
import zipfile

logging.getLogger('numba').setLevel(logging.WARNING)

# This does two things:
# 1. Makes the logger look good in a log file
# 2. Changes a bit how torch pins memory when copying to GPU, which allows you to more easily run models in parallel with an estimated 1-5% time hit
scvi.settings.reset_logging_handler()
scvi.settings.dl_pin_memory_gpu_training = False

with open(configs_file) as f: 
    data = f.read()

configs = json.loads(data)
sandbox_mode = configs["sandbox_mode"] == "True"
integration_level = configs["integration_level"]
tissue_integration = integration_level == "tissue"

sys.path.append(configs["code_path"])
from utils import *
from vdj_utils import *
from logger import SimpleLogger

output_destination = configs["output_destination"]
output_prefix = configs["output_prefix"]
s3_access_file = configs["s3_access_file"]

VARIABLE_CONFIG_KEYS = ["data_owner","s3_access_file","code_path","output_destination"] # config changes only to these fields will not initialize a new configs version
sc.settings.verbosity = 3   # verbosity: errors (0), warnings (1), info (2), hints (3)

# apply the aws credentials to allow access though aws cli; make sure the user is authorized to run in non-sandbox mode if applicable
s3_dict = set_access_keys(s3_access_file, return_dict = True)
assert sandbox_mode or hashlib.md5(bytes(s3_dict["AWS_SECRET_ACCESS_KEY"], 'utf-8')).hexdigest() in AUTHORIZED_EXECUTERS, "You are not authorized to run this script in a non sandbox mode; please set sandbox_mode to True"
set_access_keys(s3_access_file)

# create a new directory for the data and outputs
data_dir = os.path.join(output_destination, output_prefix)
os.system("mkdir -p " + data_dir)

# check for previous versions of integrated data
s3_url = "s3://immuneaging/scanvi_integrated_samples/{}_level".format(integration_level)
is_new_version, version = get_configs_status(configs, f'{s3_url}/{output_prefix}', f'integrate_using_scanvi.configs.{output_prefix}', VARIABLE_CONFIG_KEYS, data_dir)
output_configs_file = f"integrate_using_scanvi.configs.{output_prefix}.{version}.txt"

# set up logger
logger_file = "integrate_using_scanvi.{}.{}.log".format(output_prefix,version)
logger_file_path = os.path.join(data_dir, logger_file)
if os.path.isfile(logger_file_path):
    os.remove(logger_file_path)

output_h5ad_file = "{}.{}.h5ad".format(output_prefix, version)
output_h5ad_file_unstim = "{}.unstim.{}.h5ad".format(output_prefix, version)
output_h5ad_file_stim = "{}.stim.{}.h5ad".format(output_prefix, version)

logger = SimpleLogger(filename = logger_file_path)
logger.add_to_log("Running integrate_using_scanvi.py...")
logger.add_to_log("Starting time: {}".format(get_current_time()))
with open(integrate_using_scanvi_script, "r") as f:
    logger.add_to_log("integrate_using_scanvi.py md5 checksum: {}\n".format(hashlib.md5(bytes(f.read(), 'utf-8')).hexdigest()))

logger.add_to_log("using the following configurations:\n{}".format(str(configs)))
logger.add_to_log("Configs version: " + version)
logger.add_to_log("New configs version: " + str(is_new_version))

if is_new_version:
    if not sandbox_mode:
        logger.add_to_log("Uploading new configs version to S3...")
        with open(os.path.join(data_dir,output_configs_file), 'w') as f:
            json.dump(configs, f)
        sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{} --exclude "*" --include {}'.format(data_dir, s3_url, output_prefix, version, output_configs_file)
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
else:
    logger.add_to_log("Checking if h5ad file already exists on S3...")
    h5ad_file_exists = False
    logger_file_exists = False
    ls_cmd = "aws s3 ls {}/{}/{} --recursive".format(s3_url, output_prefix, version)
    files = os.popen(ls_cmd).read()
    logger.add_to_log("aws response: {}\n".format(files))
    for f in files.rstrip().split('\n'):
        if f.split('/')[-1] == output_h5ad_file:
            h5ad_file_exists = True
        if f.split('/')[-1] == logger_file:
            logger_file_exists = True
    if h5ad_file_exists and logger_file_exists:
        logger.add_to_log("The following h5ad file is already on S3: {}\nTerminating execution.".format(output_h5ad_file))
        sys.exit()
    if h5ad_file_exists and not logger_file_exists:
        logger.add_to_log("The following h5ad file is already on S3: {}\nhowever, log file is missing on S3; proceeding with execution.".format(output_h5ad_file))
    if not h5ad_file_exists:
        logger.add_to_log("The following h5ad file does not exist on S3: {}".format(output_h5ad_file))

samples = read_immune_aging_sheet("Samples")

logger.add_to_log("Downloading model and data files of integrated tissues/compartments from S3...")
integrated_object_version = configs["latest_integrated_object_version"]
# for each of stim, unstim and stim_unstim samples we run a separate integration
stim_integrated_files = []
unstim_integrated_files = []
stim_unstim_integrated_files = []

# example integrated files for the T compartment:
# T.stim.v4.h5ad
# T.stim.v4.scvi_model_batch_key_donor_id.zip
# T.stim.v4.scvi_model_batch_key_donor_id+tissue.zip
#
# example integrated files for the BLO tissue:
# BLO.stim.v4.h5ad
# BLO.stim.v4.scvi_model.zip

def append_stim_unsim(stim_status: str, file: str):
    if stim_status == "stim":
        stim_integrated_files.append(file)
    elif stim_status == "unstim":
        unstim_integrated_files.append(file)
    else:
        stim_unstim_integrated_files.append(file)

integration_modes = ["stim", "unstim", "stim_unstim"]
batch_keys = configs["batch_key"].split(",")
for mode in integration_modes:
    mode_prefix = f"{mode}." if mode != "stim_unstim" else ""
    # add h5ad
    integrated_h5ad_file = f"{output_prefix}.{mode_prefix}{integrated_object_version}.h5ad" # e.g. T.stim.v4.h5ad or BLO.stim.v4.h5ad
    append_stim_unsim(mode, integrated_h5ad_file)
    # add model(s)
    if tissue_integration:
        integrated_model_file = f"{output_prefix}.{mode_prefix}{integrated_object_version}.scvi_model.zip" # e.g. BLO.stim.v4.scvi_model.zip
        append_stim_unsim(mode, integrated_model_file)
    else:
        for bk in batch_keys:
            integrated_model_file = f"{output_prefix}.{mode_prefix}{integrated_object_version}.scvi_model_batch_key_{bk}.zip" # e.g. T.stim.v4.scvi_model_batch_key_donor_id.zip
            append_stim_unsim(mode, integrated_model_file)

# sync down everything
all_files = stim_integrated_files + unstim_integrated_files + stim_unstim_integrated_files
for file in all_files:
    s3_url_integrated = "s3://immuneaging/integrated_samples/{}_level".format(integration_level)
    sync_cmd = f'aws s3 sync --no-progress {s3_url_integrated}/{output_prefix}/{integrated_object_version} {data_dir} --exclude "*" --include {file}'
    logger.add_to_log("syncing {}...".format(file))
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    aws_response = os.popen(sync_cmd).read()
    logger.add_to_log("aws response: {}\n".format(aws_response))
    if not os.path.exists(os.path.join(data_dir, file)):
        if ".stim." in file:
            logger.add_to_log(f"file {file} not found on aws. Will skip stim-only integration", level="warning")
            if "stim" in integration_modes:
                integration_modes.remove("stim")
        elif ".unstim." in file:
            logger.add_to_log(f"file {file} not found on aws. Will skip unstim-only integration", level="warning")
            if "unstim" in integration_modes:
                integration_modes.remove("unstim")
        else:
            logger.add_to_log(f"file {file} not found on aws. Terminating execution.", level="error")
            sys.exit()

ann_version = configs["latest_annotated_object_version"]
logger.add_to_log(f"Downloading the latest annotation csv file (version: {ann_version}) from S3...")
annotation_csv_file = f"annotated_barcodes.csv"
aws_sync(f"s3://immuneaging/annotated_objects/{ann_version}", data_dir, annotation_csv_file, logger)
annotation_csv_path = os.path.join(data_dir, annotation_csv_file)
if not os.path.exists(annotation_csv_path):
    logger.add_to_log(f"annotated_barcodes.csv file does not exist on aws. Terminating execution.")
    sys.exit()

############################################
###### SAMPLE INTEGRATION BEGINS HERE ######
############################################

def get_h5ad_and_models(files: List[str]):
    h5ad_files = [f for f in files if f.endswith(".h5ad")]
    assert len(h5ad_files) == 1
    h5ad_file = h5ad_files[0]
    model_files = [f for f in files if not f.endswith(".h5ad")]
    assert len(model_files) > 0
    return h5ad_file, model_files

# run the following integration pipeline up to three times - once per integration mode (stim, unstim and stim_unstim)
valid_libs = get_all_libs("GEX")
for integration_mode in integration_modes:
    logger.add_to_log("Running {} integration pipeline...".format(integration_mode))
    if integration_mode == "stim_unstim":
        h5ad_file, _ = get_h5ad_and_models(stim_unstim_integrated_files)
        mode_suffix = ""
    elif integration_mode == "unstim":
        h5ad_file, _ = get_h5ad_and_models(unstim_integrated_files)
        mode_suffix = ".unstim"
    else:
        h5ad_file, _ = get_h5ad_and_models(stim_integrated_files)
        mode_suffix = ".stim"

    logger.add_to_log(f"Reading integrated h5ad file ({h5ad_file})...")
    adata = anndata.read_h5ad(os.path.join(data_dir, h5ad_file))
    prefix = output_prefix + mode_suffix

    logger.add_to_log(f"Adding any available manual labels to the adata...")
    annotations = pd.read_csv(annotation_csv_path, index_col="barcode")
    labels_key = configs["labels_key"]
    unlabeled_category = configs["unlabeled_category"]
    add_annotations_to_adata(adata, labels_key, unlabeled_category, annotations, valid_libs)

    n_labeled = adata[adata.obs[labels_key] != unlabeled_category].n_obs
    pct_labeled = 100*(n_labeled / adata.n_obs)
    logger.add_to_log(f"A total of {adata.n_obs} cells (and {adata.n_vars} genes) are available in the dataset of which {n_labeled} (i.e. {pct_labeled:.2f}%) have manual labels.")
    logger.add_to_log(f"Manual label distribution:\n{adata.obs[labels_key].value_counts()}")
    if n_labeled == 0:
        logger.add_to_log(f"0 cells have a label, so we cannot run scanvi (and no point in doing so). Moving on...")
        continue
    
    try:
        # just log if we have proteins though we wont be running totalvi
        if "protein_expression" in adata.obsm:
            logger.add_to_log("Detected Antibody Capture features.")
        # iterate over batch keys
        scanvi_model_files = {}
        for batch_key in batch_keys:
            logger.add_to_log("Running for batch_key {}...".format(batch_key))
            if tissue_integration:
                model_file = f"{prefix}.{integrated_object_version}.scvi_model.zip"
                data_file = f"{prefix}.{integrated_object_version}.scvi_model.data.h5ad"
            else:
                model_file = f"{prefix}.{integrated_object_version}.scvi_model_batch_key_{batch_key}.zip"
                data_file = f"{prefix}.{integrated_object_version}.scvi_model_batch_key_{batch_key}.data.h5ad"
            # unzip the model file
            logger.add_to_log(f"Unzipping model file {model_file}...".format(batch_key))
            model_path = os.path.join(data_dir, model_file)
            shutil.unpack_archive(model_path, data_dir)
            model_path = model_path.replace(f"{integrated_object_version}.", "")[:-4] # remove the .zip
            # # rename the adata file to adata.h5ad so that the model load call can pick it up
            # old_name = os.path.join(model_path, data_file)
            # new_name = os.path.join(model_path, "adata.h5ad")
            # shutil.move(old_name, new_name)

            logger.add_to_log("Loading the pre-trained scvi model and creating a scanvi model from it...")
            vae_data = anndata.read_h5ad(os.path.join(model_path, data_file))
            add_annotations_to_adata(vae_data, labels_key, unlabeled_category, annotations, valid_libs)
            vae = scvi.model.SCVI.load(model_path, adata=vae_data)
            scvi.model.SCANVI.setup_anndata(
                vae.adata,
                batch_key=batch_key,
                labels_key=labels_key,
            )
            lvae = scvi.model.SCANVI.from_scvi_model(
                vae,
                unlabeled_category,
                adata=vae.adata,
                n_latent=configs["n_latent"],
            )

            logger.add_to_log("Training scanvi model...")
            lvae.train(n_samples_per_label=configs["n_samples_per_label"])
            if "elbo_train" in lvae.history_:
                logger.add_to_log(f'Number of scanvi training epochs: {len(lvae.history_["elbo_train"])}...')
            logger.add_to_log("Saving scanvi latent representation...")
            latent = lvae.get_latent_representation()
            latent_key = f"X_scanvi_integrated_batch_key_{batch_key}"
            adata.obsm[latent_key] = latent
            model_file = "{}.{}.scanvi_model_batch_key_{}.zip".format(prefix, version, batch_key)
            logger.add_to_log("Saving the model into {}...".format(model_file))
            model_file_path = os.path.join(data_dir, model_file)
            model_dir_path = os.path.join(data_dir,"{}.scanvi_model_batch_key_{}/".format(prefix, batch_key))
            if os.path.isdir(model_dir_path):
                os.system("rm -r " + model_dir_path)
            lvae.save(model_dir_path)
            # save the data used for fitting the scanvi model separately because we need to call `write_anndata_with_object_cols`
            logger.add_to_log("Saving the data used for fitting the model...")
            data_file = "{}.{}.scanvi_model_batch_key_{}.data.h5ad".format(prefix, version, batch_key)
            write_anndata_with_object_cols(lvae.adata, model_dir_path, data_file)
            # zip the dir with all the model outputs
            zipf = zipfile.ZipFile(model_file_path, 'w', zipfile.ZIP_DEFLATED)
            zipdir(model_dir_path, zipf)
            zipf.close()
            scanvi_model_files[batch_key] = model_file
            # done with training scanvi
            logger.add_to_log("Calculate neighbors graph and UMAP based on scanvi components...")
            neighbors_key = f"scanvi_integrated_neighbors_batch_key_{batch_key}"
            sc.pp.neighbors(adata, n_neighbors=configs["neighborhood_graph_n_neighbors"], use_rep=latent_key, key_added=neighbors_key) 
            adata.obsm[f"X_umap_scvi_integrated_batch_key_{batch_key}"] = sc.tl.umap(
                adata,
                min_dist=configs["umap_min_dist"],
                spread=float(configs["umap_spread"]),
                n_components=configs["umap_n_components"],
                neighbors_key=neighbors_key,
                copy=True
            ).obsm["X_umap"]
    except Exception as err:
        logger.add_to_log("Execution failed with the following error: {}.\n{}".format(err, traceback.format_exc()), "critical")
        logger.add_to_log("Terminating execution prematurely.")
        if not sandbox_mode:
            # upload log to S3
            sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{}/ --exclude "*" --include {}'.format( \
                data_dir, s3_url, output_prefix, version, logger_file)
            os.system(sync_cmd)
        print(err)
        sys.exit()
    logger.add_to_log("Using CellTypist for annotations...")
    leiden_resolutions = [float(j) for j in configs["leiden_resolutions"].split(",")]
    celltypist_model_urls = configs["celltypist_model_urls"].split(",")
    celltypist_dotplot_min_frac = configs["celltypist_dotplot_min_frac"]
    logger.add_to_log("Downloading CellTypist models...")
    # download celltypist models
    celltypist_model_paths = []
    for celltypist_model_url in celltypist_model_urls:
        celltypist_model_path = os.path.join(data_dir, celltypist_model_url.split("/")[-1])
        urllib.request.urlretrieve(celltypist_model_url, celltypist_model_path)
        celltypist_model_paths.append(celltypist_model_path)
    dotplot_paths = []
    for batch_key in batch_keys:
        dotplot_paths += annotate(
            adata,
            model_paths = celltypist_model_paths,
            model_urls = celltypist_model_urls,
            components_key = f"X_scanvi_integrated_batch_key_{batch_key}",
            neighbors_key = f"neighbors_scanvi",
            n_neighbors = configs["neighborhood_graph_n_neighbors"],
            resolutions = leiden_resolutions,
            model_name = f"scanvi_batch_key_{batch_key}" + mode_suffix,
            dotplot_min_frac = celltypist_dotplot_min_frac,
            logger = logger,
            save_all_outputs = True
        )
    dotplot_dirname = "dotplots" + mode_suffix
    dotplot_dir = os.path.join(data_dir, dotplot_dirname)
    os.system("rm -r -f {}".format(dotplot_dir))
    os.system("mkdir {}".format(dotplot_dir))
    for dotplot_path in dotplot_paths:
        os.system("mv {} {}".format(dotplot_path, dotplot_dir))
    dotplots_zipfile = "{}.{}.celltypist_dotplots.zip".format(prefix, version)
    zipf = zipfile.ZipFile(os.path.join(data_dir,dotplots_zipfile), 'w', zipfile.ZIP_DEFLATED)
    zipdir(dotplot_dir, zipf)
    zipf.close()
    logger.add_to_log("Saving h5ad files...")
    output_h5ad_file = "{}.{}.h5ad".format(prefix, version)
    write_anndata_with_object_cols(adata, data_dir, output_h5ad_file)
    # OUTPUT UPLOAD TO S3 - ONLY IF NOT IN SANDBOX MODE
    if not sandbox_mode:
        logger.add_to_log("Uploading h5ad file to S3...")
        sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{} --exclude "*" --include {}'.format(
            data_dir, s3_url, output_prefix, version, output_h5ad_file)
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
        logger.add_to_log("Uploading model files (a single .zip file for each model) and CellTypist dot plots to S3...")
        inclusions = [f"--include {file}" for file in list(scanvi_model_files.values())]
        sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{}/ --exclude "*" {}'.format(
            data_dir, s3_url, output_prefix, version, " ".join(inclusions))
        sync_cmd += ' --include {}'.format(dotplots_zipfile)         
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
    logger.add_to_log("Number of cells: {}, number of genes: {}.".format(adata.n_obs, adata.n_vars))

logger.add_to_log("Execution of integrate_samples.py is complete.")

logging.shutdown()
if not sandbox_mode:
    # Uploading log file to S3.
    sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{} --exclude "*" --include {}'.format(
        data_dir, s3_url, output_prefix, version, logger_file)
    os.system(sync_cmd)
