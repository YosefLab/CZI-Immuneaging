import sys

from anndata._core.anndata import AnnData
process_sample_script = sys.argv[0]
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
import time
import celltypist
import urllib.request
import pickle

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
output_destination = configs["output_destination"]
donor = configs["donor"]
seq_run = configs["seq_run"]
sample_id = configs["sample_id"]
library_type = configs["library_type"]

sys.path.append(configs["code_path"])
from utils import *

AUTHORIZED_EXECUTERS = ["b750bd0287811e901c88dc328187e25f", "1c75133ab6a1fc3ed9233d3fe40b3d73"] # md5 checksums of the AWS_SECRET_ACCESS_KEY value of those that are authorized to run process_library and upload outputs to the server; note that individuals with upload permission to aws can bypass that by changing the code - this is just designed to alert users that they should only use sandbox mode.
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
assert sandbox_mode or hashlib.md5(bytes(s3_dict["AWS_SECRET_ACCESS_KEY"], 'utf-8')).hexdigest() in AUTHORIZED_EXECUTERS, "You are not authorized to run this script in a non sanbox mode; please set sandbox_mode to True"
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
    cp_cmd = "cp {} {}".format(configs_file, os.path.join(data_dir,output_configs_file))
    os.system(cp_cmd)
    if not sandbox_mode:
        add_to_log("Uploading new configs version to S3...")
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
processed_libraries_dir = configs["processed_libraries_dir"]
if len(processed_libraries_dir) > 0:
    add_to_log("Copying h5ad files of processed libraries from {}...".format(processed_libraries_dir))
    cp_cmd = "cp -r {}/ {}".format(processed_libraries_dir.rstrip("/"), data_dir)
    os.system(cp_cmd)
else:
    add_to_log("Downloading h5ad files of processed libraries from S3...")
    add_to_log("*** Note: This can take some time. If you already have the processed libraries, you can halt this process and provide processed_libraries_dir in the config file in order to use your existing h5ad files. ***")   
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

summary = ["\n{0}\nExecution summary\n{0}".format("="*25)]

add_to_log("Downloading the Donors sheet from the Google spreadsheets...")
donors = read_immune_aging_sheet("Donors")
add_to_log("Downloading the Samples sheet from the Google spreadsheets...")
samples = read_immune_aging_sheet("Samples")

############################################
###### SAMPLE PROCESSING BEGINS HERE #######
############################################

alerts = []

add_to_log("Reading h5ad files of processed libraries...")
adata_dict = {}
initial_n_obs = 0
solo_genes = set()
library_ids_keep = []
for j in range(len(library_ids)):
    library_id = library_ids[j]
    library_version = library_versions[j]
    lib_h5ad_file = os.path.join(data_dir, "{}_{}_{}_{}.processed.{}.h5ad".format(donor, seq_run,
        library_type, library_id, library_version))
    adata_dict[library_id] = sc.read_h5ad(lib_h5ad_file)
    adata_dict[library_id].obs["library_id"] = library_id
    if "Classification" in adata_dict[library_id].obs.columns:
        adata_dict[library_id] = adata_dict[library_id][adata_dict[library_id].obs["Classification"] == sample_id]
        if "min_cells_per_library" in configs and configs["min_cells_per_library"] > adata_dict[library_id].n_obs:
            # do not consider cells from this library
            msg = "Cells from library {} were not included - there are {} cells from the sample, however, min_cells_per_library was set to {}.".format(
                library_id,adata_dict[library_id].n_obs,configs["min_cells_per_library"])
            alerts.append(msg)
            del adata_dict[library_id]
            continue
    library_ids_keep.append(library_id)
    solo_genes_j = np.logical_and(sc.pp.filter_genes(adata_dict[library_id], min_cells=configs["solo_filter_genes_min_cells"], inplace=False)[0], (adata_dict[library_id].var["feature_types"] == "Gene Expression").values)
    solo_genes_j = adata_dict[library_id].var.index[solo_genes_j]
    initial_n_obs = initial_n_obs + adata_dict[library_id].n_obs
    if len(solo_genes)==0:
        solo_genes = solo_genes_j
    else:
        solo_genes = np.intersect1d(solo_genes,solo_genes_j)

library_ids = library_ids_keep

if len(library_ids)==0:
    for i in alerts:
        add_to_log(i)
    add_to_log("No cells passed the filtering steps. Terminating execution.")
    logging.shutdown()
    if not sandbox_mode:
        # Uploading log file to S3...
        sync_cmd = 'aws s3 sync {} s3://immuneaging/processed_samples/{}/{}/ --exclude "*" --include {}'.format(
            data_dir, prefix, version, logger_file)
        os.system(sync_cmd)
    sys.exit()

add_to_log("Concatenating all cells of sample {} from available libraries...".format(sample_id))
adata = adata_dict[library_ids[0]]
if len(library_ids) > 1:
    adata = adata.concatenate([adata_dict[library_ids[j]] for j in range(1,len(library_ids))])

add_to_log("A total of {} cells and {} genes were found.".format(adata.n_obs, adata.n_vars))
summary.append("Started with a total of {} cells and {} genes coming from {} libraries.".format(
    initial_n_obs, adata.n_vars, len(library_ids)))

add_to_log("Adding metadata...")
donor_index = donors["Donor ID"] == donor
add_to_log("Adding donor-level metadata...")
for k in DONORS_FIELDS.keys():
    adata.obs[DONORS_FIELDS[k]] = donors[k][donor_index].values[0]

add_to_log("Adding sample-level metadata...")
sample_index = samples["Sample_ID"] == sample_id
for k in SAMPLES_FIELDS.keys():
    adata.obs[SAMPLES_FIELDS[k]] = samples[k][sample_index].values[0]

if library_type == "GEX":
    adata.obs["GEX_chem"] = samples["GEX chem"][sample_index].values[0]
    adata.obs["HTO_chem"] = samples["HTO chem"][sample_index].values[0]
    adata.obs["ADT_chem"] = samples["CITE chem"][sample_index].values[0]

if adata.n_obs > 0:
    no_cells = False
    add_to_log("Current number of cells: {}.".format(adata.n_obs))
    add_to_log("Current number of genes: {}.".format(adata.n_vars))
else:
    no_cells = True
    add_to_log("Detected no cells. Skipping data processing steps.")
    summary.append("Detected no cells; skipped data processing steps.")

if not no_cells:
    try:
        # save raw counts
        is_cite = "Antibody Capture" in np.unique(adata.var.feature_types)
        if is_cite:
            add_to_log("Detected Antibody Capture features.")
            # Get the CITE data. Only keep the subset of proteins that are not HTO tags. We already store these in an obs column
            hto_tag = configs["donor"]+"-"
            protein = adata[:, (adata.var["feature_types"] == "Antibody Capture") & ~(adata.var_names.str.startswith(hto_tag))].copy()
            # copy protein data from X into adata.obsm["protein_expression"]
            protein_df = protein.to_df()
            protein_expression_obsm_key = "protein_expression"
            if not protein_df.empty:
                adata.obsm[protein_expression_obsm_key] = protein_df
            else:
                add_to_log("All detected Antibody Capture features were due to HTO, no proteins of interest to analyze.")
                is_cite = False
            rna = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
        else:
            rna = adata.copy()
        add_to_log("Running decontX for estimating contamination levels from ambient RNA...")
        decontx_data_dir = os.path.join(data_dir,"decontx")
        os.system("mkdir -p " + decontx_data_dir)
        raw_counts_file = os.path.join(decontx_data_dir, "{}.raw_counts.txt".format(prefix))
        decontaminated_counts_file = os.path.join(decontx_data_dir, "{}.decontx.decontaminated.txt".format(prefix))
        contamination_levels_file = os.path.join(decontx_data_dir, "{}.decontx.contamination.txt".format(prefix))
        decontx_model_file = os.path.join(decontx_data_dir, "{}.processed.{}.decontx_model.RData".format(prefix, version))
        r_script_file = os.path.join(decontx_data_dir, "{}.decontx.script.R".format(prefix))
        df = pd.DataFrame(rna.X.todense().T)
        df.index = rna.var.index
        df.columns = rna.obs.index
        df.to_csv(raw_counts_file)
        if len(library_ids)>1:
            batch_key = "batch"
            batch_file = os.path.join(decontx_data_dir, "{}.batch.txt".format(prefix))
            pd.DataFrame(rna.obs[batch_key].values.astype(str)).to_csv(batch_file, header=False, index=False)
        else:
            batch_key = None
            batch_file = None
        # R commands for running and outputing decontx
        l = ["library('celda')",
            "x <- as.matrix(read.csv('{}', row.names=1))".format(raw_counts_file),
            "batch <- if ('{0}' == 'None') NULL else as.character(read.table('{0}', header=FALSE)$V1)".format(batch_file),
            "res <- decontX(x=x, batch=batch)",
            "write.table(res$contamination, file ='{}',quote = FALSE,row.names = FALSE,col.names = FALSE)".format(contamination_levels_file),
            "write.table(as.matrix(res$decontXcounts), file ='{}',quote = FALSE,sep = ',')".format(decontaminated_counts_file),
            "decontx_model <- list('estimates'=res$estimates, 'z'= res$z)",
            "save(decontx_model, file='{}')".format(decontx_model_file)]
        with open(r_script_file,'w') as f: 
            f.write("\n".join(l))
        add_to_log("Running the script in {}".format(decontx_model_file))
        os.system("Rscript " + r_script_file)
        add_to_log("Adding decontaminated counts and contamination levels to data object...")
        contamination_levels = pd.read_csv(contamination_levels_file, index_col=0, header=None).index
        decontaminated_counts = pd.read_csv(decontaminated_counts_file).transpose()
        decontaminated_counts.index = rna.obs.index # note that decontx replaces "-" with "." in the cell names
        rna.obs["contamination_levels"] = contamination_levels
        rna.X = np.round(decontaminated_counts)
        # keep only the genes in solo_genes, which is required to prevent errors with solo in case some genes are expressed in only a subset of the batches.
        rna = rna[:,rna.var.index.isin(solo_genes)]
        # remove empty cells after decontaminations
        n_obs_before = rna.n_obs
        rna = rna[rna.X.sum(axis=1) >= configs["filter_decontaminated_cells_min_genes"],:]
        n_decon_cells_filtered = n_obs_before-rna.n_obs
        msg = "Removed {} cells with total decontaminated counts below filter_decontaminated_cells_min_genes={}".format(
            n_decon_cells_filtered, configs["filter_decontaminated_cells_min_genes"])
        add_to_log(msg)
        summary.append(msg)
        add_to_log("Detecting highly variable genes...")
        sc.pp.highly_variable_genes(rna, n_top_genes=configs["n_highly_variable_genes"], subset=True,
            flavor=configs["highly_variable_genes_flavor"], batch_key=batch_key, span = 1.0)
        if is_cite:
            rna, totalvi_model, totalvi_model_file = run_model(rna, configs, batch_key, protein_expression_obsm_key, "totalvi", prefix, version, data_dir)
        rna, scvi_model, scvi_model_file = run_model(rna, configs, batch_key, None, "scvi", prefix, version, data_dir)
        add_to_log("Running solo for detecting doublets...")
        if len(library_ids)>1:
            batches = pd.unique(rna.obs[batch_key])
            add_to_log("Running solo on the following batches separately: {}".format(batches))
            is_solo_singlet = np.ones((rna.n_obs,), dtype=bool)
            for batch in batches:
                add_to_log("Running solo on batch {}...".format(batch))
                solo_batch = scvi.external.SOLO.from_scvi_model(scvi_model, restrict_to_batch=batch)
                solo_batch.train(max_epochs=configs["solo_max_epochs"])
                is_solo_singlet[(rna.obs["batch"] == batch).values] = solo_batch.predict(soft=False) == "singlet"
                rna.obs["is_solo_singlet"] = is_solo_singlet
        else:
            add_to_log("Running solo...")
            solo = scvi.external.SOLO.from_scvi_model(scvi_model)
            solo.train(max_epochs=configs["solo_max_epochs"])
            is_solo_singlet = solo.predict(soft=False) == "singlet"
        add_to_log("Removing doublets...")
        n_obs_before = rna.n_obs
        rna = rna[is_solo_singlet,]
        add_to_log("Removed {} estimated doublets; {} droplets remained.".format(n_obs_before-rna.n_obs,rna.n_obs))
        summary.append("Removed {} estimated doublets.".format(n_obs_before-rna.n_obs))
        if rna.n_obs == 0:
            add_to_log("No cells left after doublet detection. Skipping the next processing steps.")
            summary.append("No cells left after doublet detection.")
        else:
            add_to_log("Normalizing RNA...")
            sc.pp.normalize_total(rna, target_sum=configs["normalize_total_target_sum"])
            sc.pp.log1p(rna)
            sc.pp.scale(rna)
            rna.raw = rna
            add_to_log("Calculating PCA...")
            sc.pp.pca(rna)
            add_to_log("Calculating neighborhood graph and UMAP based on PCA...")
            key = "pca_neighbors"
            sc.pp.neighbors(rna, n_neighbors=configs["neighborhood_graph_n_neighbors"],
                use_rep="X_pca", key_added=key) 
            rna.obsm["X_umap_pca"] = sc.tl.umap(rna, min_dist=configs["umap_min_dist"], spread=float(configs["umap_spread"]),
                n_components=configs["umap_n_components"], neighbors_key=key, copy=True).obsm["X_umap"]
            add_to_log("Calculating neighborhood graph and UMAP based on SCVI components...")
            key = "scvi_neighbors"
            sc.pp.neighbors(rna, n_neighbors=configs["neighborhood_graph_n_neighbors"],
                use_rep="X_scVI", key_added=key) 
            rna.obsm["X_umap_scvi"] = sc.tl.umap(rna, min_dist=configs["umap_min_dist"], spread=float(configs["umap_spread"]),
                n_components=configs["umap_n_components"], neighbors_key=key, copy=True).obsm["X_umap"]
            add_to_log("Calculating neighborhood graph and UMAP based on TOTALVI components...")
            if is_cite:
                key = "totalvi_neighbors"
                sc.pp.neighbors(rna, n_neighbors=configs["neighborhood_graph_n_neighbors"],
                    use_rep="X_totalVI", key_added=key) 
                rna.obsm["X_umap_totalvi"] = sc.tl.umap(rna, min_dist=configs["umap_min_dist"], spread=float(configs["umap_spread"]),
                    n_components=configs["umap_n_components"], neighbors_key=key, copy=True).obsm["X_umap"]
        add_to_log("Gathering data...")
        # keep only the cells that passed all filters
        keep = adata.obs.index.isin(rna.obs.index)
        adata = adata[keep,]
        adata.obsm.update(rna.obsm)
        adata.obs[rna.obs.columns] = rna.obs
        add_to_log("Remove protein counts from adata.X...")
        # keep only the subset of X that does not have protein data, since we already moved these to obsm
        adata = adata[:, adata.var["feature_types"] == "Gene Expression"]
        # save raw rna counts
        adata.layers["raw_counts"] = adata.X.copy()
        # save decontaminated counts (only applies to rna data; consider only cells that we keep after filters)
        adata.layers["decontaminated_counts"] = adata.layers["raw_counts"]
        adata.layers["decontaminated_counts"][:,adata.var.index.isin(decontaminated_counts.columns)] = np.array(decontaminated_counts.loc[adata.obs.index])
        if adata.n_obs > 0:
            add_to_log("Normalize rna counts in adata.X...")
            sc.pp.normalize_total(adata, target_sum=configs["normalize_total_target_sum"])
            sc.pp.log1p(adata)
            add_to_log("Predict cell type labels using celltypist...")
            model_urls = configs["celltypist_model_urls"].split(",")
            # run prediction using every specified model (url)
            for i in range(len(model_urls)):
                model_file = model_urls[i].split("/")[-1]
                model_path = os.path.join(data_dir,model_file)
                # download reference data
                urllib.request.urlretrieve(model_urls[i], model_path)
                model = celltypist.models.Model.load(model = model_path)
                # for some reason celltypist changes adata in a way that then doesn't allow to copy it (which is needed later); a fix is to use a copy of adata.
                adata_copy = adata.copy()
                if configs["normalize_total_target_sum"] != 10000:
                    # renormalize the data in the copy of adata in case configs["normalize_total_target_sum"] != 10000 (which is the scale required by celltypist)
                    add_to_log("normalizing data for celltypist...")
                    adata_copy.X = adata.layers["raw_counts"]
                    sc.pp.normalize_total(adata_copy, target_sum=10000)
                    sc.pp.log1p(adata_copy)
                predictions = celltypist.annotate(adata_copy, model = model, majority_voting = True)
                add_to_log("Saving celltypist annotations for model {}...".format(model_file))
                adata.obs["celltypist_predicted_labels."+str(i+1)] = predictions.predicted_labels["predicted_labels"]
                adata.obs["celltypist_over_clustering."+str(i+1)] = predictions.predicted_labels["over_clustering"]
                adata.obs["celltypist_majority_voting."+str(i+1)] = predictions.predicted_labels["majority_voting"]
                adata.obs["celltypist_model."+str(i+1)] = model_urls[i]
    except Exception as err:
        add_to_log("Execution failed with the following error:\n{}".format(err))
        add_to_log("Terminating execution prematurely.")
        if not sandbox_mode:
            # upload log to S3
            sync_cmd = 'aws s3 sync {} s3://immuneaging/processed_samples/{}/{}/ --exclude "*" --include {}'.format(
                data_dir, prefix, version, logger_file)
            os.system(sync_cmd)
        print(err)
        sys.exit()

add_to_log("Saving h5ad file...")
adata.write(os.path.join(data_dir,h5ad_file))

###############################################################
###### OUTPUT UPLOAD TO S3 - ONLY IF NOT IN SANDBOX MODE ######
###############################################################

if not sandbox_mode:
    add_to_log("Uploading h5ad file to S3...")
    sync_cmd = 'aws s3 sync {} s3://immuneaging/processed_samples/{}/{}/ --exclude "*" --include {}'.format(
        data_dir, prefix, version, h5ad_file)
    add_to_log("sync_cmd: {}".format(sync_cmd))
    add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
    if not no_cells:
        add_to_log("Uploading model files (a single .zip file for each model) to S3...")
        sync_cmd = 'aws s3 sync {} s3://immuneaging/processed_samples/{}/{}/ --exclude "*" --include {} --include {}'.format(
            data_dir, prefix, version, scvi_model_file, totalvi_model_file)
        add_to_log("sync_cmd: {}".format(sync_cmd))
        add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))        
        add_to_log("Uploading decontx model file to S3...")
        sync_cmd = 'aws s3 sync {} s3://immuneaging/processed_samples/{}/{}/ --exclude "*" --include {}'.format(
            decontx_data_dir, prefix, version, decontx_model_file.split("/")[-1])
        add_to_log("sync_cmd: {}".format(sync_cmd))
        add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

add_to_log("Execution of process_sample.py is complete.")

summary.append("Final number of cells: {}, final number of genes: {}.".format(
    adata.n_obs, adata.n_vars))
for i in summary:
    add_to_log(i)

if len(alerts) == 0:
    add_to_log("No alerts.")
else:
    add_to_log("Alerts:")
    for i in alerts:
        add_to_log(i)

logging.shutdown()
if not sandbox_mode:
    # Uploading log file to S3...
    sync_cmd = 'aws s3 sync {} s3://immuneaging/processed_samples/{}/{}/ --exclude "*" --include {}'.format(
        data_dir, prefix, version, logger_file)
    os.system(sync_cmd)