## For speed, this script should be executed on s130 which has a GPU. However, note that integrating a large number of samples may require high memory, which may not always be available on s130.
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
import celltypist
import urllib.request
import zipfile
import gc

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
tissue_integration = configs["integration_level"] == "tissue"

sys.path.append(configs["code_path"])
from utils import *
from vdj_utils import *
from logger import SimpleLogger

output_destination = configs["output_destination"]
s3_access_file = configs["s3_access_file"]

# sort the sample ids lexicographically (and accordingly the versions) in order to avoid generating a new version if the exact same set of samples were previously used but in a different order
order = np.argsort(configs["sample_ids"].split(","))
all_sample_ids = np.array(configs["sample_ids"].split(","))[order]
processed_sample_configs_version = np.array(configs["processed_sample_configs_version"].split(","))[order]
# the followings are required because we check if the integration of the requested samples already exist on aws
configs["sample_ids"] = ",".join(all_sample_ids)
configs["processed_sample_configs_version"] = ",".join(processed_sample_configs_version)

assert(len(all_sample_ids)>1)

VARIABLE_CONFIG_KEYS = ["data_owner","s3_access_file","code_path","output_destination"] # config changes only to these fields will not initialize a new configs version
sc.settings.verbosity = 3   # verbosity: errors (0), warnings (1), info (2), hints (3)

# apply the aws credentials to allow access though aws cli; make sure the user is authorized to run in non-sandbox mode if applicable
s3_dict = set_access_keys(s3_access_file, return_dict = True)
assert sandbox_mode or hashlib.md5(bytes(s3_dict["AWS_SECRET_ACCESS_KEY"], 'utf-8')).hexdigest() in AUTHORIZED_EXECUTERS, "You are not authorized to run this script in a non sandbox mode; please set sandbox_mode to True"
set_access_keys(s3_access_file)

# create a new directory for the data and outputs
data_dir = os.path.join(output_destination, configs["output_prefix"])
os.system("mkdir -p " + data_dir)

# check for previous versions of integrated data
s3_url = "s3://immuneaging/integrated_samples/{}_level".format(configs["integration_level"])
is_new_version, version = get_configs_status(configs, s3_url + "/" + configs["output_prefix"], "integrate_samples.configs."+configs["output_prefix"],
    VARIABLE_CONFIG_KEYS, data_dir)
output_configs_file = "integrate_samples.configs.{}.{}.txt".format(configs["output_prefix"],version)

# set up logger
logger_file = "integrate_samples.{}.{}.log".format(configs["output_prefix"],version)
logger_file_path = os.path.join(data_dir, logger_file)
if os.path.isfile(logger_file_path):
    os.remove(logger_file_path)

output_h5ad_file = "{}.{}.h5ad".format(configs["output_prefix"], version)
output_h5ad_model_file = "{}.{}.model_data.h5ad".format(configs["output_prefix"], version)

output_h5ad_file_unstim = "{}.unstim.{}.h5ad".format(configs["output_prefix"], version)
output_h5ad_model_file_unstim = "{}.unstim.{}.model_data.h5ad".format(configs["output_prefix"], version)

output_h5ad_file_stim = "{}.stim.{}.h5ad".format(configs["output_prefix"], version)
output_h5ad_model_file_stim = "{}.stim.{}.model_data.h5ad".format(configs["output_prefix"], version)

logger = SimpleLogger(filename = logger_file_path)
logger.add_to_log("Running integrate_samples.py...")
logger.add_to_log("Starting time: {}".format(get_current_time()))
with open(integrate_samples_script, "r") as f:
    logger.add_to_log("integrate_samples.py md5 checksum: {}\n".format(hashlib.md5(bytes(f.read(), 'utf-8')).hexdigest()))

logger.add_to_log("using the following configurations:\n{}".format(str(configs)))
logger.add_to_log("Configs version: " + version)
logger.add_to_log("New configs version: " + str(is_new_version))

if is_new_version:
    if not sandbox_mode:
        logger.add_to_log("Uploading new configs version to S3...")
        with open(os.path.join(data_dir,output_configs_file), 'w') as f:
            json.dump(configs, f)
        sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{} --exclude "*" --include {}'.format(
            data_dir, s3_url, configs["output_prefix"], version, output_configs_file)
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
else:
    logger.add_to_log("Checking if h5ad file already exists on S3...")
    h5ad_file_exists = False
    logger_file_exists = False
    ls_cmd = "aws s3 ls {}/{}/{} --recursive".format(s3_url, configs["output_prefix"],version)
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

logger.add_to_log("Downloading h5ad files of processed samples from S3...")
all_h5ad_files = []
# for collecting the sample IDs of unstim and stim samples for which we will generate an additional, separate integration:
unstim_sample_ids = [] 
unstim_h5ad_files = []
stim_sample_ids = [] 
stim_h5ad_files = []
for j in range(len(all_sample_ids)):
    sample_id = all_sample_ids[j]
    sample_version = processed_sample_configs_version[j]
    sample_h5ad_file = "{}_GEX.processed.{}.h5ad".format(sample_id,sample_version)
    sample_h5ad_path = os.path.join(data_dir,sample_h5ad_file)
    all_h5ad_files.append(sample_h5ad_path)
    stim_status = samples["Stimulation"][samples["Sample_ID"] == sample_id].values[0]
    if stim_status == "Nonstim":
        unstim_sample_ids.append(sample_id)
        unstim_h5ad_files.append(sample_h5ad_path)
    else:
        stim_sample_ids.append(sample_id)
        stim_h5ad_files.append(sample_h5ad_path)
    if not os.path.exists(sample_h5ad_path):
        sync_cmd = 'aws s3 sync --no-progress s3://immuneaging/processed_samples/{}_GEX/{}/ {} --exclude "*" --include {}'.format(
        sample_id,sample_version,data_dir,sample_h5ad_file)
        logger.add_to_log("syncing {}...".format(sample_h5ad_file))
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        aws_response = os.popen(sync_cmd).read()
        logger.add_to_log("aws response: {}\n".format(aws_response))
    if not os.path.exists(sample_h5ad_path):
        logger.add_to_log("h5ad file does not exist on aws for sample {}. Terminating execution.".format(sample_id))
        sys.exit()

if not tissue_integration:
    compartment = configs["output_prefix"]
    logger.add_to_log(f"Downloading csv file containing cell barcodes for the {compartment} compartment from S3...")
    barcodes_csv_file = f"{compartment}-compartment-barcodes.csv"
    aws_sync("s3://immuneaging/per-compartment-barcodes/", data_dir, barcodes_csv_file, logger)
    barcodes_csv_path = os.path.join(data_dir, barcodes_csv_file)
    if not os.path.exists(barcodes_csv_path):
        logger.add_to_log(f"per-compartment-barcodes csv file does not exist on aws for compartment {compartment}. Terminating execution.")
        sys.exit()
    logger.add_to_log(f"Downloading csv file containing cell barcodes to exclude, if any...")
    exclude_barcodes_csv_file = f"Exclude-barcodes.csv"
    aws_sync("s3://immuneaging/per-compartment-barcodes/", data_dir, exclude_barcodes_csv_file, logger)
    exclude_barcodes_csv_path = os.path.join(data_dir, exclude_barcodes_csv_file)

############################################
###### SAMPLE INTEGRATION BEGINS HERE ######
############################################

# run the following integration pipeline three times if there is a combination of stim and unstim samples - once using the stim and unstim samples, once using unstim only, once using stim only
integration_modes = ["stim_unstim"]
if (len(unstim_sample_ids) > 0) and (len(all_sample_ids)-len(unstim_sample_ids)>0):
    integration_modes += ["unstim", "stim"]

valid_libs = get_all_libs("GEX")
for integration_mode in integration_modes:
    logger.add_to_log("Running {} integration pipeline...".format(integration_mode))
    if integration_mode == "stim_unstim":
        sample_ids = all_sample_ids
        h5ad_files = all_h5ad_files
        mode_suffix = ""
    elif integration_mode == "unstim":
        sample_ids = unstim_sample_ids
        h5ad_files = unstim_h5ad_files
        mode_suffix = ".unstim"
    else:
        sample_ids = stim_sample_ids
        h5ad_files = stim_h5ad_files
        mode_suffix = ".stim"
    prefix = configs["output_prefix"] + mode_suffix
    dotplot_dirname = "dotplots" + mode_suffix
    output_h5ad_file = "{}.{}.h5ad".format(prefix, version)
    logger.add_to_log("Reading h5ad files of processed samples...")
    barcodes = pd.read_csv(barcodes_csv_path, header=None)
    exclude_barcodes = pd.read_csv(exclude_barcodes_csv_path, header=None) if os.path.isfile(exclude_barcodes_csv_path) else None
    try:
        if adata_dict is not None:
            logger.add_to_log("Cleaning up adata_dict to free up some memory before the next run...")
            del adata_dict
            gc.collect()
    except NameError:
        pass
    adata_dict = {}
    for j in range(len(h5ad_files)):
        h5ad_file = h5ad_files[j]
        sample_id = sample_ids[j]
        if tissue_integration:
            adata_dict[sample_id] = sc.read_h5ad(h5ad_file)
        else:
            adata_temp = sc.read_h5ad(h5ad_file)
            # replace adata with the subset of adata that is limited to the compartment barcodes
            trimmed_adata_index = adata_temp.obs.index.map(lambda x: strip_integration_markers(x, valid_libs))
            idx = trimmed_adata_index.isin(barcodes[0])
            if exclude_barcodes is not None:
                # also, exclude any cells that are marked as Exclude
                # Note: The majority of cells to be excluded are removed during earlier stages of the pipeline
                # such as sample and library processing. However, we have this additional step here in case that
                # we find cells to exclude later on and don't want to go back and re-run earlier processing stages
                # (ideally though we should in order to keep everything consistent) or if for any other reason we
                # need to exclude some cells here.
                idx = idx & ~trimmed_adata_index.isin(exclude_barcodes[0])
            if np.sum(idx) > 0: # i.e. if there is at least one cell that passes the condition above
                adata_dict[sample_id] = adata_temp[idx, :].copy()
                del adata_temp
                gc.collect()
            else:
                del adata_temp
                gc.collect()
                continue
        # add the tissue to the adata since we may need it as part of a composite batch_key later
        adata_dict[sample_id].obs["tissue"] = samples["Organ"][samples["Sample_ID"] == sample_id].values[0]
        adata_dict[sample_id].obs["sample_id"] = sample_id
        for field in ('X_pca', 'X_scVI', 'X_totalVI', 'X_umap_pca', 'X_umap_scvi', 'X_umap_totalvi'):
            if field in adata_dict[sample_id].obsm:
                del adata_dict[sample_id].obsm[field]
    if tissue_integration and len(adata_dict) == 0:
        msg = f"No cells found for compartment {compartment} from all samples in {integration_mode} integration mode..."
        if integration_mode != "stim_unstim":
            logger.add_to_log(msg, level="warning")
            # move on to the next integration mode
            continue
        else:
            logger.add_to_log(msg + " Terminating execution.", level="error")
            logging.shutdown()
            if not sandbox_mode:
                # upload log to S3
                aws_sync(data_dir, "{}/{}/{}/".format(s3_url, configs["output_prefix"], version), logger_file, logger, do_log=False)
            sys.exit()
    sample_ids = list(adata_dict.keys())
    # get the names of all proteins and control proteins across all samples (in case the samples were collected using more than one protein panel)
    proteins = set()
    proteins_ctrl = set()
    for sample_id in sample_ids:
        if "protein_expression" in adata_dict[sample_id].obsm:
            proteins.update(adata_dict[sample_id].obsm["protein_expression"].columns)
        if "protein_expression_Ctrl" in adata_dict[sample_id].obsm:
            proteins_ctrl.update(adata_dict[sample_id].obsm["protein_expression_Ctrl"].columns)
    # each sample should include all proteins; missing values are set as NaN
    if len(proteins) > 0 or len(proteins_ctrl) > 0:
        proteins = list(proteins)
        proteins_ctrl = list(proteins_ctrl)
        for sample_id in sample_ids:
            df = pd.DataFrame(columns = proteins, index = adata_dict[sample_id].obs.index)
            if "protein_expression" in adata_dict[sample_id].obsm:
                df[adata_dict[sample_id].obsm["protein_expression"].columns] = adata_dict[sample_id].obsm["protein_expression"].copy()
            adata_dict[sample_id].obsm["protein_expression"] = df
            df_ctrl = pd.DataFrame(columns = proteins_ctrl, index = adata_dict[sample_id].obs.index)
            if "protein_expression_Ctrl" in adata_dict[sample_id].obsm:
                df_ctrl[adata_dict[sample_id].obsm["protein_expression_Ctrl"].columns] = adata_dict[sample_id].obsm["protein_expression_Ctrl"].copy()
            adata_dict[sample_id].obsm["protein_expression_Ctrl"] = df_ctrl
    logger.add_to_log("Concatenating all datasets...")
    adata = adata_dict[sample_ids[0]]
    if len(sample_ids) > 1:
        adata = adata.concatenate([adata_dict[sample_ids[j]] for j in range(1,len(sample_ids))], join="outer", index_unique=None)
    # Move the summary statistics of the genes (under .var) to a separate csv file
    cols_to_varm = [j for j in adata.var.columns if "n_cells" in j] + \
    [j for j in adata.var.columns if "mean_counts" in j] + \
    [j for j in adata.var.columns if "pct_dropout_by_counts" in j] + \
    [j for j in adata.var.columns if "total_counts" in j]
    output_gene_stats_csv_file = "{}.{}.gene_stats.csv".format(prefix, version)
    adata.var.iloc[:,adata.var.columns.isin(cols_to_varm)].to_csv(os.path.join(data_dir,output_gene_stats_csv_file))
    adata.var = adata.var.drop(labels = cols_to_varm, axis = "columns")
    # protein QC
    protein_levels_max_sds = configs["protein_levels_max_sds"] if "protein_levels_max_sds" in configs else None
    if (len(proteins) > 0) and protein_levels_max_sds is not None:
        logger.add_to_log("Running protein QC...")
        n_cells_before, n_proteins_before = adata.obsm["protein_expression"].shape
        # (1) remove cells that demonstrate extremely high (or low) protein library size; normalize library size by the total number of non NaN proteins available for the cell.
        prot_exp = adata.obsm["protein_expression"]
        normalized_cell_lib_size = prot_exp.sum(axis=1)/(prot_exp.shape[1] - prot_exp.isnull().sum(axis=1))
        keep, lower_bound, upper_bound = detect_outliers(normalized_cell_lib_size, protein_levels_max_sds)
        # save the qc thresholds in .uns so they can later be applied to other data if needed (specifically, in case of integrating new data into the model we learn here)
        adata.uns["protein_qc"] = {}
        adata.uns["protein_qc"]["normalized_cell_lib_size_lower_bound"] = lower_bound
        adata.uns["protein_qc"]["normalized_cell_lib_size_upper_bound"] = upper_bound
        logger.add_to_log("Removing {} cells with extreme protein values (normalized library size more extreme than {} standard deviations)...".format(np.sum(~keep), protein_levels_max_sds))
        adata = adata[keep,].copy()
        # (2) remove proteins that have extreme library size (across all cells; consider only cells that have non-missing values and normalize by the number of cells with no missing values for the protein)
        prot_exp = adata.obsm["protein_expression"]
        normalized_protein_lib_sizes = prot_exp.sum()/(prot_exp.shape[0] - prot_exp.isnull().sum(axis=0))
        keep, lower_bound, upper_bound = detect_outliers(normalized_protein_lib_sizes, protein_levels_max_sds)
        adata.uns["protein_qc"]["normalized_protein_lib_size_lower_bound"] = lower_bound
        adata.uns["protein_qc"]["normalized_protein_lib_size_upper_bound"] = upper_bound
        logger.add_to_log("Removing {} proteins with extreme total number of reads across cells (normalized levels - by the number of cells with no missing values for the protein - more extreme than {} standard deviations)...".format(
            np.sum(~keep), protein_levels_max_sds))
        cols_keep = adata.obsm["protein_expression"].columns[keep.values]
        cols_remove = adata.obsm["protein_expression"].columns[~keep.values]
        logger.add_to_log("Removing proteins: {} ".format(", ".join(cols_remove.values)))
        extend_removed_features_df(adata, "removed_proteins", adata.obsm["protein_expression"][cols_remove])
        adata.obsm["protein_expression"] = adata.obsm["protein_expression"][cols_keep]
    # remove proteins that were left with zero reads
    non_zero_proteins = adata.obsm["protein_expression"].sum() > 0
    num_zero_proteins = np.sum(~non_zero_proteins)
    cols_keep = adata.obsm["protein_expression"].columns[non_zero_proteins.values]
    cols_remove = adata.obsm["protein_expression"].columns[~non_zero_proteins.values]
    if num_zero_proteins > 0:
        logger.add_to_log("Removing {} proteins with zero reads across all cells...".format(num_zero_proteins))
        logger.add_to_log("Removing proteins: {} ".format(", ".join(cols_remove.values)))
    extend_removed_features_df(adata, "removed_proteins", adata.obsm["protein_expression"][cols_remove])
    adata.obsm["protein_expression"] = adata.obsm["protein_expression"][cols_keep]
    # summarize protein qc
    logger.add_to_log("Protein QC summary: a total of {} cells ({}%) and {} proteins ({}%) were filtered out owing to extreme values.".format(
        n_cells_before-adata.n_obs, round(100*(n_cells_before-adata.n_obs)/n_cells_before,2),
        n_proteins_before-adata.obsm["protein_expression"].shape[1], round(100*(n_proteins_before-adata.obsm["protein_expression"].shape[1])/n_proteins_before,2)))
    # end protein QC
    if "is_solo_singlet" in adata.obs:
        del adata.obs["is_solo_singlet"]
    if "Classification" in adata.obs:
        del adata.obs["Classification"]
    logger.add_to_log("A total of {} cells and {} genes are available after merge.".format(adata.n_obs, adata.n_vars))
    # set the train size; this train size was justified by an experiment that is described here https://yoseflab.github.io/scvi-tools-reproducibility/runtime_analysis_reproducibility/
    configs["train_size"] = 0.9 if 0.1 * adata.n_obs < 20000 else 1-(20000/adata.n_obs)
    if (not tissue_integration) and adata.n_obs > 200000:
        configs["early_stopping"] = "True"
    try:
        is_cite = "protein_expression" in adata.obsm
        if is_cite:
            logger.add_to_log("Detected Antibody Capture features.")
        # iterate over batch keys
        batch_keys = ["batch"] if "batch_key" not in configs else configs["batch_key"].split(",")
        scvi_model_files = {}
        totalvi_model_files = {}
        run_pca = True
        for batch_key in batch_keys:
            logger.add_to_log("Running for batch_key {}...".format(batch_key))
            rna = adata.copy()
            if batch_key not in rna.obs:
                if batch_key == "donor_id+tissue":
                    rna.obs[batch_key] = rna.obs["donor_id"].astype("str") + "+" + rna.obs["tissue"].astype("str")
                    # we need to remove cells, if any, that belong to a batch that has a size one
                    # or else the scanpy hvg call below fails. It is ok if there is only ever one or
                    # two such cells, but not if there is a lot of them, which is why we emit a warning
                    # log.
                    batch_vc = rna.obs["donor_id+tissue"].value_counts()
                    for b in batch_vc[batch_vc == 1].index.values:
                        barcode = rna[rna.obs[batch_key] == b].obs.index[0]
                        logger.add_to_log(f"Removing cell {barcode} where {batch_key} = {b} b/c it's the only one of its batch.", level="warning")
                        rna = rna[rna.obs[batch_key] != b,:].copy()
                        # need to also remove it from adata
                        adata = adata[adata.obs[batch_key] != b,:].copy()
                else:
                    logger.add_to_log(f"Batch key {batch_key} not found in adata columns. Terminating execution.", level="error")
                    logging.shutdown()
                    if not sandbox_mode:
                        # upload log to S3
                        aws_sync(data_dir, "{}/{}/{}/".format(s3_url, configs["output_prefix"], version), logger_file, logger, do_log=False)
                    sys.exit()
            logger.add_to_log("Filtering out vdj genes...")
            rna = filter_vdj_genes(rna, configs["vdj_genes"], data_dir, logger)
            rna.layers["rounded_decontaminated_counts"] = rna.layers["decontaminated_counts"].copy()
            rna.layers["rounded_decontaminated_counts"].data = np.round(rna.layers["rounded_decontaminated_counts"].data)
            rna.X = rna.layers["rounded_decontaminated_counts"].copy()
            sc.pp.log1p(rna)
            rna.layers["log1p_transformed"] = rna.X.copy()
            logger.add_to_log("Detecting highly variable genes...")
            if configs["highly_variable_genes_flavor"] == "seurat_v3":
                # highly_variable_genes requires counts data in this case
                rna.X = rna.layers["rounded_decontaminated_counts"]
            sc.pp.highly_variable_genes(
                rna,
                n_top_genes=configs["n_highly_variable_genes"],
                subset=True,
                flavor=configs["highly_variable_genes_flavor"],
                batch_key=batch_key,
                span = 1.0
            )
            rna.X = rna.layers["rounded_decontaminated_counts"]
            # scvi
            key = f"X_scvi_integrated_batch_key_{batch_key}"
            _, scvi_model_file = run_model(rna, configs, batch_key, None, "scvi", prefix, version, data_dir, logger, key)
            scvi_model_files[batch_key] = scvi_model_file
            logger.add_to_log("Calculate neighbors graph and UMAP based on scvi components...")
            neighbors_key = f"scvi_integrated_neighbors_batch_key_{batch_key}"
            sc.pp.neighbors(rna, n_neighbors=configs["neighborhood_graph_n_neighbors"], use_rep=key, key_added=neighbors_key) 
            rna.obsm[f"X_umap_scvi_integrated_batch_key_{batch_key}"] = sc.tl.umap(
                rna,
                min_dist=configs["umap_min_dist"],
                spread=float(configs["umap_spread"]),
                n_components=configs["umap_n_components"],
                neighbors_key=neighbors_key,
                copy=True
            ).obsm["X_umap"]
            if is_cite:
                # totalVI
                key = f"X_totalVI_integrated_batch_key_{batch_key}"
                # if there is no protein information for some of the cells set them to zero (instead of NaN)
                rna.obsm["protein_expression"] = rna.obsm["protein_expression"].fillna(0)
                # there are known spurious failures with totalVI (such as "invalid parameter loc/scale")
                # so we try a few times then carry on with the rest of the script as we can still mine the
                # rest of the data regardless of CITE info
                retry_count = 4
                try:
                    _, totalvi_model_file = run_model(rna, configs, batch_key, "protein_expression", "totalvi", prefix, version, data_dir, logger, latent_key=key, max_retry_count=retry_count)
                    totalvi_model_files[batch_key] = totalvi_model_file
                    logger.add_to_log("Calculate neighbors graph and UMAP based on totalVI components...")
                    neighbors_key = f"totalvi_integrated_neighbors_batch_key_{batch_key}"
                    sc.pp.neighbors(rna, n_neighbors=configs["neighborhood_graph_n_neighbors"],use_rep=key, key_added=neighbors_key) 
                    rna.obsm[f"X_umap_totalvi_integrated_batch_key_{batch_key}"] = sc.tl.umap(
                        rna,
                        min_dist=configs["umap_min_dist"],
                        spread=float(configs["umap_spread"]),
                        n_components=configs["umap_n_components"],
                        neighbors_key=neighbors_key,
                        copy=True
                    ).obsm["X_umap"]
                except Exception as err:
                    logger.add_to_log("Execution of totalVI failed with the following error (latest) with retry count {}: {}. Moving on...".format(retry_count, err), "warning")
            if run_pca:
                logger.add_to_log("Calculating PCA...")
                rna.X = rna.layers["log1p_transformed"]
                sc.pp.pca(rna)
                logger.add_to_log("Calculating neighborhood graph and UMAP based on PCA...")
                key = "pca_neighbors"
                sc.pp.neighbors(rna, n_neighbors=configs["neighborhood_graph_n_neighbors"], use_rep="X_pca", key_added=key) 
                rna.obsm["X_umap_pca"] = sc.tl.umap(
                    rna,
                    min_dist=configs["umap_min_dist"],
                    spread=float(configs["umap_spread"]),
                    n_components=configs["umap_n_components"],
                    neighbors_key=key,
                    copy=True
                ).obsm["X_umap"]
                run_pca = False
            # update the adata with the components of the dim reductions and umap coordinates
            adata.obsm.update(rna.obsm)
            # save the identity of the most variable genes used
            adata.var[f"is_highly_variable_gene_batch_key_{batch_key}"] = adata.var.index.isin(rna.var.index)
    except Exception as err:
        logger.add_to_log("Execution failed with the following error: {}.\n{}".format(err, traceback.format_exc()), "critical")
        logger.add_to_log("Terminating execution prematurely.")
        if not sandbox_mode:
            # upload log to S3
            sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{}/ --exclude "*" --include {}'.format( \
                data_dir, s3_url, configs["output_prefix"], version, logger_file)
            os.system(sync_cmd)
        print(err)
        sys.exit()
    logger.add_to_log("Using CellTypist for annotations...")
    # remove celltypist predictions and related metadata that were added at the sample-level processing
    celltypist_cols = [j for j in adata.obs.columns if "celltypist" in j]
    adata.obs = adata.obs.drop(labels = celltypist_cols, axis = "columns")
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
            components_key = f"X_scvi_integrated_batch_key_{batch_key}",
            neighbors_key = f"neighbors_scvi",
            n_neighbors = configs["neighborhood_graph_n_neighbors"],
            resolutions = leiden_resolutions,
            model_name = f"scvi_batch_key_{batch_key}" + mode_suffix,
            dotplot_min_frac = celltypist_dotplot_min_frac,
            logger = logger,
            save_all_outputs = True
        )
        totalvi_key = f"X_totalVI_integrated_batch_key_{batch_key}"
        if totalvi_key in adata.obsm:
            dotplot_paths += annotate(
                adata,
                model_paths = celltypist_model_paths,
                model_urls = celltypist_model_urls,
                components_key = totalvi_key,
                neighbors_key = "neighbors_totalvi",
                n_neighbors = configs["neighborhood_graph_n_neighbors"],
                resolutions = leiden_resolutions,
                model_name = f"totalvi_batch_key_{batch_key}" + mode_suffix,
                dotplot_min_frac = celltypist_dotplot_min_frac,
                logger = logger,
            )
    dotplot_dir = os.path.join(data_dir,dotplot_dirname)
    os.system("rm -r -f {}".format(dotplot_dir))
    os.system("mkdir {}".format(dotplot_dir))
    for dotplot_path in dotplot_paths:
        os.system("mv {} {}".format(dotplot_path, dotplot_dir))
    dotplots_zipfile = "{}.{}.celltypist_dotplots.zip".format(prefix, version)
    zipf = zipfile.ZipFile(os.path.join(data_dir,dotplots_zipfile), 'w', zipfile.ZIP_DEFLATED)
    zipdir(dotplot_dir, zipf)
    zipf.close()
    logger.add_to_log("Adding bcr_lib_id and tcr_lib_id to adata where applicable...")
    gex_to_bcr, gex_to_tcr = get_gex_lib_to_vdj_lib_mapping()
    add_vdj_lib_ids_to_adata(adata, gex_to_bcr, gex_to_tcr)
    logger.add_to_log("Saving h5ad files...")
    adata.obs["age"] = adata.obs["age"].astype(str)
    adata.obs["BMI"] = adata.obs["BMI"].astype(str)
    adata.obs["height"] = adata.obs["height"].astype(str)
    write_anndata_with_object_cols(adata, data_dir, output_h5ad_file)
    # OUTPUT UPLOAD TO S3 - ONLY IF NOT IN SANDBOX MODE
    if not sandbox_mode:
        logger.add_to_log("Uploading h5ad file to S3...")
        sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{} --exclude "*" --include {}'.format(
            data_dir, s3_url, configs["output_prefix"], version, output_h5ad_file)
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
        logger.add_to_log("Uploading model files (a single .zip file for each model) and CellTypist dot plots to S3...")
        inclusions = [f"--include {file}" for file in list(scvi_model_files.values()) + list(totalvi_model_files.values())]
        sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{}/ --exclude "*" {}'.format(
            data_dir, s3_url, configs["output_prefix"], version, " ".join(inclusions))
        sync_cmd += ' --include {}'.format(dotplots_zipfile)         
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
        logger.add_to_log("Uploading gene stats csv file to S3...")
        sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{} --exclude "*" --include {}'.format(
            data_dir, s3_url, configs["output_prefix"], version, output_gene_stats_csv_file)
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
    logger.add_to_log("Number of cells: {}, number of genes: {}.".format(adata.n_obs, adata.n_vars))

logger.add_to_log("Execution of integrate_samples.py is complete.")

logging.shutdown()
if not sandbox_mode:
    # Uploading log file to S3.
    sync_cmd = 'aws s3 sync --no-progress {} {}/{}/{} --exclude "*" --include {}'.format(
        data_dir, s3_url, configs["output_prefix"], version, logger_file)
    os.system(sync_cmd)
