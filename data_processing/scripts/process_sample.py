#########################################################
###### INITIALIZATIONS AND PREPARATIONS BEGIN HERE ######
#########################################################

import sys
import logging
import os
import json
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse.csr import csr_matrix
import scvi
import hashlib
import celltypist
import urllib.request
import traceback
import scirpy as ir

logging.getLogger('numba').setLevel(logging.WARNING)

process_sample_script = sys.argv[0]
configs_file = sys.argv[1]

sc.settings.verbosity = 3   # verbosity: errors (0), warnings (1), info (2), hints (3)

with open(configs_file) as f: 
    data = f.read()	

configs = json.loads(data)
sandbox_mode = configs["sandbox_mode"] == "True"
output_destination = configs["output_destination"]
donor = configs["donor"]
seq_run = configs["seq_run"]
sample_id = configs["sample_id"]

sys.path.append(configs["code_path"])

from utils import *
from logger import SimpleLogger
init_scvi_settings()

# config changes only to these fields will not initialize a new configs version
VARIABLE_CONFIG_KEYS = ["data_owner","s3_access_file","code_path","output_destination"]

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

# apply the aws credentials to allow access though aws cli; make sure the user is authorized to run in non-sandbox mode if applicable
s3_dict = set_access_keys(configs["s3_access_file"], return_dict = True)
assert sandbox_mode or hashlib.md5(bytes(s3_dict["AWS_SECRET_ACCESS_KEY"], 'utf-8')).hexdigest() in AUTHORIZED_EXECUTERS, "You are not authorized to run this script in a non sanbox mode; please set sandbox_mode to True"
set_access_keys(configs["s3_access_file"])

# create a new directory for the data and outputs
data_dir = os.path.join(output_destination, "_".join([donor, seq_run]))
os.system("mkdir -p " + data_dir)

# check for previous versions of the processed sample
prefix = sample_id
s3_path = "s3://immuneaging/processed_samples/" + prefix
is_new_version, version = get_configs_status(configs, s3_path, "process_sample.configs." + prefix,
    VARIABLE_CONFIG_KEYS, data_dir)
output_configs_file = "process_sample.configs.{}.{}.txt".format(prefix,version)

# set up logger
logger_file = "process_sample.{}.{}.log".format(prefix,version)
logger_file_path = os.path.join(data_dir, logger_file)
if os.path.isfile(logger_file_path):
    os.remove(logger_file_path)

logger = SimpleLogger(filename = logger_file_path)
logger.add_to_log("Running process_sample.py...")
logger.add_to_log(QC_STRING_START_TIME.format(get_current_time()))
with open(process_sample_script, "r") as f:
    logger.add_to_log("process_sample.py md5 checksum: {}\n".format(hashlib.md5(bytes(f.read(), 'utf-8')).hexdigest()))

logger.add_to_log("using the following configurations:\n{}".format(str(configs)))
logger.add_to_log("Configs version: " + version)
logger.add_to_log("New configs version: " + str(is_new_version))

h5ad_file = "{}.processed.{}.h5ad".format(prefix, version)
if is_new_version:
    cp_cmd = "cp {} {}".format(configs_file, os.path.join(data_dir,output_configs_file))
    os.system(cp_cmd)
    if not sandbox_mode:
        logger.add_to_log("Uploading new configs version to S3...")
        sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/processed_samples/{}/{}/ --exclude "*" --include {}'.format(
            data_dir, prefix, version, output_configs_file)
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
else:
    logger.add_to_log("Checking if h5ad file already exists on S3...")
    h5ad_file_exists = False
    logger_file_exists = False
    ls_cmd = "aws s3 ls s3://immuneaging/processed_samples/{}/{} --recursive".format(prefix,version)
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

library_ids = configs["library_ids"].split(',')
library_types = configs["library_types"].split(',')
library_versions = configs["processed_library_configs_version"].split(',')
if ("processed_libraries_dir" in configs) and (len(configs.get("processed_libraries_dir")) > 0):
    processed_libraries_dir = configs["processed_libraries_dir"]
    logger.add_to_log("Copying h5ad files of processed libraries from {}...".format(processed_libraries_dir))
    cp_cmd = "cp -r {}/ {}".format(processed_libraries_dir.rstrip("/"), data_dir)
    os.system(cp_cmd)
else:
    logger.add_to_log("Downloading h5ad files of processed libraries from S3...")
    logger.add_to_log("*** Note: This can take some time. If you already have the processed libraries, you can halt this process and provide processed_libraries_dir in the config file in order to use your existing h5ad files. ***")   
    for j in range(len(library_ids)):
        library_id = library_ids[j]
        library_type = library_types[j]
        library_version = library_versions[j]
        lib_h5ad_file = "{}_{}_{}_{}.processed.{}.h5ad".format(donor, seq_run,
            library_type, library_id, library_version)
        sync_cmd = 'aws s3 sync --no-progress s3://immuneaging/processed_libraries/{}_{}_{}_{}/{}/ {} --exclude "*" --include {}'.format(
            donor, seq_run, library_type, library_id, library_version, data_dir, lib_h5ad_file)
        logger.add_to_log("syncing {}...".format(lib_h5ad_file))
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

summary = ["\n{0}\nExecution summary\n{0}".format("="*25)]

logger.add_to_log("Downloading the Donors sheet from the Google spreadsheets...")
donors = read_immune_aging_sheet("Donors")
logger.add_to_log("Downloading the Samples sheet from the Google spreadsheets...")
samples = read_immune_aging_sheet("Samples")

############################################
###### SAMPLE PROCESSING BEGINS HERE #######
############################################

logger.add_to_log("Reading h5ad files of processed libraries for GEX libs...")
adata_dict = {}
initial_n_obs = 0
solo_genes = set()
library_ids_gex = []
for j in range(len(library_ids)):
    if library_types[j] != "GEX":
        continue
    library_id = library_ids[j]
    library_type = library_types[j]
    library_version = library_versions[j]
    lib_h5ad_file = os.path.join(data_dir, "{}_{}_{}_{}.processed.{}.h5ad".format(donor, seq_run,
        library_type, library_id, library_version))
    adata_dict[library_id] = sc.read_h5ad(lib_h5ad_file)
    adata_dict[library_id].obs["library_id"] = library_id
    if "Classification" in adata_dict[library_id].obs.columns:
        adata_dict[library_id] = adata_dict[library_id][adata_dict[library_id].obs["Classification"] == sample_id].copy()
        if "min_cells_per_library" in configs and configs["min_cells_per_library"] > adata_dict[library_id].n_obs:
            # do not consider cells from this library
            msg = "Cells from library {} were not included - there are {} cells from the sample, however, min_cells_per_library was set to {}.".format(
                library_id,adata_dict[library_id].n_obs,configs["min_cells_per_library"])
            logger.add_to_log(msg, "warning")
            del adata_dict[library_id]
            continue
    library_ids_gex.append(library_id)
    solo_genes_j = np.logical_and(sc.pp.filter_genes(adata_dict[library_id], min_cells=configs["solo_filter_genes_min_cells"], inplace=False)[0], (adata_dict[library_id].var["feature_types"] == "Gene Expression").values)
    solo_genes_j = adata_dict[library_id].var.index[solo_genes_j]
    initial_n_obs += adata_dict[library_id].n_obs
    if len(solo_genes)==0:
        solo_genes = solo_genes_j
    else:
        solo_genes = np.intersect1d(solo_genes,solo_genes_j)

if len(library_ids_gex)==0:
    logger.add_to_log("No cells passed the filtering steps. Terminating execution.", "error")
    logging.shutdown()
    if not sandbox_mode:
        # Uploading log file to S3...
        sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/processed_samples/{}/{}/ --exclude "*" --include {}'.format(
            data_dir, prefix, version, logger_file)
        os.system(sync_cmd)
    sys.exit()

# TODO make sure concatenation doesn't produce a "batch" column that would cause issues downstream

logger.add_to_log("Concatenating all cells of sample {} from available GEX libraries...".format(sample_id))
adata = adata_dict[library_ids_gex[0]]
if len(library_ids_gex) > 1:
    adata = adata.concatenate([adata_dict[library_ids_gex[j]] for j in range(1,len(library_ids_gex))])

def build_adata_from_ir_libs(lib_type: str) -> AnnData:
    assert lib_type in ["BCR", "TCR"]
    logger.add_to_log("Reading h5ad files of processed libraries for {} libs...".format(lib_type))
    adata_dict = {}
    library_ids_ir = []
    for j in range(len(library_ids)):
        if library_types[j] != lib_type:
            continue
        library_id = library_ids[j]
        library_type = library_types[j]
        library_version = library_versions[j]
        lib_h5ad_file = os.path.join(data_dir, "{}_{}_{}_{}.processed.{}.h5ad".format(donor, seq_run,
            library_type, library_id, library_version))
        adata_dict[library_id] = sc.read_h5ad(lib_h5ad_file)
        adata_dict[library_id].obs["library_id"] = library_id
        library_ids_ir.append(library_id)

    logger.add_to_log("Concatenating all cells of sample {} from available {} libraries...".format(sample_id, lib_type))
    adata = adata_dict[library_ids_ir[0]]
    if len(library_ids_ir) > 1:
        adata = adata.concatenate([adata_dict[library_ids_ir[j]] for j in range(1,len(library_ids_ir))])
    return adata

adata_bcr = build_adata_from_ir_libs("BCR")
adata_tcr = build_adata_from_ir_libs("TCR")

# TODO use merge_airr_chains?
# TODO use scirpy.tl.chain_qc() here too?
adata_bcr.obs = adata_bcr.obs.join(adata_tcr.obs, how="outer").copy()
adata_ir = adata_bcr.copy()

ir.pp.merge_with_ir(adata, adata_ir)

logger.add_to_log("A total of {} cells and {} genes were found.".format(adata.n_obs, adata.n_vars))
summary.append("Started with a total of {} cells and {} genes coming from {} libraries.".format(
    initial_n_obs, adata.n_vars, len(library_ids)))

logger.add_to_log("Adding metadata...")
donor_index = donors["Donor ID"] == donor
logger.add_to_log("Adding donor-level metadata...")
for k in DONORS_FIELDS.keys():
    adata.obs[DONORS_FIELDS[k]] = donors[k][donor_index].values[0]

logger.add_to_log("Adding sample-level metadata...")
sample_index = samples["Sample_ID"] == sample_id
for k in SAMPLES_FIELDS.keys():
    adata.obs[SAMPLES_FIELDS[k]] = samples[k][sample_index].values[0]

if library_type == "GEX":
    adata.obs["GEX_chem"] = samples["GEX chem"][sample_index].values[0]
    adata.obs["HTO_chem"] = samples["HTO chem"][sample_index].values[0]
    adata.obs["ADT_chem"] = samples["CITE chem"][sample_index].values[0]

if adata.n_obs > 0:
    no_cells = False
    logger.add_to_log("Current number of cells: {}.".format(adata.n_obs))
    logger.add_to_log("Current number of genes: {}.".format(adata.n_vars))
else:
    no_cells = True
    logger.add_to_log("Detected no cells. Skipping data processing steps.", "error")
    summary.append("Detected no cells; skipped data processing steps.")

if not no_cells:
    try:
        # save raw counts
        is_cite = "Antibody Capture" in np.unique(adata.var.feature_types)
        if is_cite:
            logger.add_to_log("Detected Antibody Capture features.")
            # Get the CITE data. Only keep the subset of proteins that are not HTO tags. We already store these in an obs column
            hto_tag = configs["donor"]+"-"
            protein = adata[:, (adata.var["feature_types"] == "Antibody Capture") & ~(adata.var_names.str.startswith(hto_tag))].copy()
            # copy protein data from X into adata.obsm["protein_expression"]
            protein_df = protein.to_df()
            protein_expression_obsm_key = "protein_expression"
            protein_expression_ctrl_obsm_key = "protein_expression_Ctrl"
            if not protein_df.empty:
                if np.median(protein_df.fillna(0).sum(axis=1)) == 0:
                    logger.add_to_log("median coverage (total number of protein reads per cell) across cells is 0. Removing protein information from data.", level = "warning")
                    is_cite = False
                else:
                    # switch the protein names to their internal names defined in the protein panels (in the Google Spreadsheet)
                    protein_df.columns = get_internal_protein_names(protein_df)
                    # save control and non-control proteins in different obsm structures
                    is_ctrl_protein = np.array([i.endswith("Ctrl") for i in protein_df.columns])
                    adata.obsm[protein_expression_ctrl_obsm_key] = protein_df[protein_df.columns[is_ctrl_protein]].copy()
                    adata.obsm[protein_expression_obsm_key] = protein_df[protein_df.columns[np.logical_not(is_ctrl_protein)]].copy()
            else:
                logger.add_to_log("All detected Antibody Capture features were due to HTO, no proteins of interest to analyze.")
                is_cite = False
            rna = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
        else:
            rna = adata.copy()
        logger.add_to_log("Running decontX for estimating contamination levels from ambient RNA...")
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
        logger.add_to_log("Running the script in {}".format(decontx_model_file))
        os.system("Rscript " + r_script_file)
        logger.add_to_log("Adding decontaminated counts and contamination levels to data object...")
        contamination_levels = pd.read_csv(contamination_levels_file, index_col=0, header=None).index
        decontaminated_counts = pd.read_csv(decontaminated_counts_file).transpose()
        decontaminated_counts.index = rna.obs.index # note that decontx replaces "-" with "." in the cell names
        rna.obs["contamination_levels"] = contamination_levels
        rna.X = csr_matrix(np.round(decontaminated_counts))
        # keep only the genes in solo_genes, which is required to prevent errors with solo in case some genes are expressed in only a subset of the batches.
        rna = rna[:,rna.var.index.isin(solo_genes)].copy()
        # remove empty cells after decontaminations
        n_obs_before = rna.n_obs
        rna = rna[rna.X.sum(axis=1) >= configs["filter_decontaminated_cells_min_genes"],:].copy()
        n_decon_cells_filtered = n_obs_before-rna.n_obs
        percent_removed = 100*n_decon_cells_filtered/n_obs_before
        level = "warning" if percent_removed > 10 else "info"
        msg = QC_STRING_AMBIENT_RNA.format(n_decon_cells_filtered, percent_removed, configs["filter_decontaminated_cells_min_genes"])
        logger.add_to_log(msg, level=level)
        summary.append(msg)
        logger.add_to_log("Filtering out vdj genes...")
        rna = filter_vdj_genes(rna, configs["vdj_genes"], data_dir, logger)
        logger.add_to_log("Detecting highly variable genes...")
        rna.layers["rounded_decontaminated_counts_copy"] = rna.X.copy()
        if configs["highly_variable_genes_flavor"] != "seurat_v3":
            # highly_variable_genes requires log-transformed data in this case
            sc.pp.log1p(rna)
        sc.pp.highly_variable_genes(rna, n_top_genes=configs["n_highly_variable_genes"], subset=True,
            flavor=configs["highly_variable_genes_flavor"], batch_key=batch_key, span = 1.0)
        rna.X = rna.layers["rounded_decontaminated_counts_copy"]
        logger.add_to_log("Predict cell type labels using celltypist...")
        model_urls = configs["celltypist_model_urls"].split(",")
        if configs["rbc_model_url"] != "":
            model_urls.append(configs["rbc_model_url"])
        # run prediction using every specified model (url)
        rbc_model_index = -1
        for i in range(len(model_urls)):
            model_file = model_urls[i].split("/")[-1]
            model_path = os.path.join(data_dir,model_file)
            # download reference data
            if model_urls[i].startswith("s3://"):
                model_folder = model_urls[i][:-len(model_file)] # remove the model_file suffix
                aws_sync(model_folder, data_dir, model_file, logger)
            else:
                urllib.request.urlretrieve(model_urls[i], model_path)
            model = celltypist.models.Model.load(model = model_path)
            # for some reason celltypist changes the anndata object in a way that then doesn't allow to copy it (which is needed later); a fix is to use a copy of the anndata object.
            rna_copy = rna.copy()
            # normalize the copied data with a scale of 10000 (which is the scale required by celltypist)
            logger.add_to_log("normalizing data for celltypist...")
            sc.pp.normalize_total(rna_copy, target_sum=10000)
            sc.pp.log1p(rna_copy)
            predictions = celltypist.annotate(rna_copy, model = model, majority_voting = True)
            # save the index for the RBC model if one exists, since we will need it further below
            if model_file.startswith("RBC_model"):
                rbc_model_index = i
            logger.add_to_log("Saving celltypist annotations for model {}, model description:\n{}".format(model_file, json.dumps(model.description, indent=2)))
            rna.obs["celltypist_predicted_labels."+str(i+1)] = predictions.predicted_labels["predicted_labels"]
            rna.obs["celltypist_over_clustering."+str(i+1)] = predictions.predicted_labels["over_clustering"]
            rna.obs["celltypist_majority_voting."+str(i+1)] = predictions.predicted_labels["majority_voting"]
            rna.obs["celltypist_model."+str(i+1)] = model_urls[i]
        # filter out RBC's
        if rbc_model_index != -1:
            n_obs_before = rna.n_obs
            rna = rna[rna.obs["celltypist_predicted_labels."+str(rbc_model_index+1)] != "RBC", :].copy()
            percent_removed = 100*(n_obs_before-rna.n_obs)/n_obs_before
            level = "warning" if percent_removed > 20 else "info"
            logger.add_to_log(QC_STRING_RBC.format(n_obs_before-rna.n_obs, percent_removed, rna.n_obs), level=level)
        if is_cite:
            # there are known spurious failures with totalVI (such as "invalid parameter loc/scale")
            # so we try a few times then carry on with the rest of the script as we can still mine the
            # rest of the data regardless of CITE info
            retry_count = 4
            try:
                _, totalvi_model_file = run_model(rna, configs, batch_key, protein_expression_obsm_key, "totalvi", prefix, version, data_dir, logger, max_retry_count=retry_count)
            except Exception as err:
                logger.add_to_log("Execution of totalVI failed with the following error (latest) with retry count {}: {}. Moving on...".format(retry_count, err), "warning")
                is_cite = False
        scvi_model, scvi_model_file = run_model(rna, configs, batch_key, None, "scvi", prefix, version, data_dir, logger)
        logger.add_to_log("Running solo for detecting doublets...")
        if len(library_ids)>1:
            batches = pd.unique(rna.obs[batch_key])
            logger.add_to_log("Running solo on the following batches separately: {}".format(batches))
            is_solo_singlet = np.ones((rna.n_obs,), dtype=bool)
            for batch in batches:
                logger.add_to_log("Running solo on batch {}...".format(batch))
                solo_batch = scvi.external.SOLO.from_scvi_model(scvi_model, restrict_to_batch=batch)
                solo_batch.train(max_epochs=configs["solo_max_epochs"])
                is_solo_singlet[(rna.obs["batch"] == batch).values] = solo_batch.predict(soft=False) == "singlet"
                rna.obs["is_solo_singlet"] = is_solo_singlet
        else:
            logger.add_to_log("Running solo...")
            solo = scvi.external.SOLO.from_scvi_model(scvi_model)
            solo.train(max_epochs=configs["solo_max_epochs"])
            is_solo_singlet = solo.predict(soft=False) == "singlet"
        logger.add_to_log("Removing doublets...")
        n_obs_before = rna.n_obs
        rna = rna[is_solo_singlet,].copy()
        percent_removed = 100*(n_obs_before-rna.n_obs)/n_obs_before
        level = "warning" if percent_removed > 40 else "info"
        logger.add_to_log(QC_STRING_DOUBLETS.format(n_obs_before-rna.n_obs, percent_removed, rna.n_obs), level=level)
        summary.append("Removed {} estimated doublets.".format(n_obs_before-rna.n_obs))
        if rna.n_obs == 0:
            logger.add_to_log("No cells left after doublet detection. Skipping the next processing steps.", "error")
            summary.append("No cells left after doublet detection.")
        else:
            logger.add_to_log("Normalizing RNA...")
            sc.pp.normalize_total(rna, target_sum=configs["normalize_total_target_sum"])
            sc.pp.log1p(rna)
            sc.pp.scale(rna)
            rna.raw = rna
            logger.add_to_log("Calculating PCA...")
            sc.pp.pca(rna)
            logger.add_to_log("Calculating neighborhood graph and UMAP based on PCA...")
            key = "pca_neighbors"
            sc.pp.neighbors(rna, n_neighbors=configs["neighborhood_graph_n_neighbors"],
                use_rep="X_pca", key_added=key) 
            rna.obsm["X_umap_pca"] = sc.tl.umap(rna, min_dist=configs["umap_min_dist"], spread=float(configs["umap_spread"]),
                n_components=configs["umap_n_components"], neighbors_key=key, copy=True).obsm["X_umap"]
            logger.add_to_log("Calculating neighborhood graph and UMAP based on SCVI components...")
            key = "scvi_neighbors"
            sc.pp.neighbors(rna, n_neighbors=configs["neighborhood_graph_n_neighbors"],
                use_rep="X_scVI", key_added=key) 
            rna.obsm["X_umap_scvi"] = sc.tl.umap(rna, min_dist=configs["umap_min_dist"], spread=float(configs["umap_spread"]),
                n_components=configs["umap_n_components"], neighbors_key=key, copy=True).obsm["X_umap"]
            logger.add_to_log("Calculating neighborhood graph and UMAP based on TOTALVI components...")
            if is_cite:
                key = "totalvi_neighbors"
                sc.pp.neighbors(rna, n_neighbors=configs["neighborhood_graph_n_neighbors"],
                    use_rep="X_totalVI", key_added=key) 
                rna.obsm["X_umap_totalvi"] = sc.tl.umap(rna, min_dist=configs["umap_min_dist"], spread=float(configs["umap_spread"]),
                    n_components=configs["umap_n_components"], neighbors_key=key, copy=True).obsm["X_umap"]
        logger.add_to_log("Gathering data...")
        # keep only the cells that passed all filters
        keep = adata.obs.index.isin(rna.obs.index)
        adata = adata[keep,].copy()
        adata.obsm.update(rna.obsm)
        adata.obs[rna.obs.columns] = rna.obs
        logger.add_to_log("Remove protein counts from adata.X...")
        # keep only the subset of X that does not have protein data, since we already moved these to obsm
        adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
        # save raw rna counts
        adata.layers["raw_counts"] = adata.X.copy()
        # save decontaminated counts (only applies to rna data; consider only cells that we keep after filters)
        adata.layers["decontaminated_counts"] = csr_matrix(decontaminated_counts.loc[adata.obs.index,adata.var.index])
        if adata.n_obs > 0:
            logger.add_to_log("Normalize rna counts in adata.X...")
            sc.pp.normalize_total(adata, target_sum=configs["normalize_total_target_sum"])
            sc.pp.log1p(adata)
    except Exception as err:
        logger.add_to_log("Execution failed with the following error: {}.\n{}".format(err, traceback.format_exc()), "critical")
        logger.add_to_log("Terminating execution prematurely.")
        if not sandbox_mode:
            # upload log to S3
            sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/processed_samples/{}/{}/ --exclude "*" --include {}'.format(
                data_dir, prefix, version, logger_file)
            os.system(sync_cmd)
        print(err)
        sys.exit()

logger.add_to_log("Saving h5ad file...")
adata.write(os.path.join(data_dir,h5ad_file), compression="lzf")

###############################################################
###### OUTPUT UPLOAD TO S3 - ONLY IF NOT IN SANDBOX MODE ######
###############################################################

if not sandbox_mode:
    logger.add_to_log("Uploading h5ad file to S3...")
    sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/processed_samples/{}/{}/ --exclude "*" --include {}'.format(data_dir, prefix, version, h5ad_file)
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
    if not no_cells:
        logger.add_to_log("Uploading model files (a single .zip file for each model) to S3...")
        sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/processed_samples/{}/{}/ --exclude "*" --include {}'.format(data_dir, prefix, version, scvi_model_file)
        if is_cite:
            # also upload the totalvi file if we had CITE data
            sync_cmd += ' --include {}'.format(totalvi_model_file)
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))        
        logger.add_to_log("Uploading decontx model file to S3...")
        sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/processed_samples/{}/{}/ --exclude "*" --include {}'.format(
            decontx_data_dir, prefix, version, decontx_model_file.split("/")[-1])
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

logger.add_to_log("Execution of process_sample.py is complete.")

summary.append(QC_STRING_COUNTS.format(adata.n_obs, adata.n_vars))
for i in summary:
    logger.add_to_log(i)

logging.shutdown()
if not sandbox_mode:
    # Uploading log file to S3...
    sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/processed_samples/{}/{}/ --exclude "*" --include {}'.format(
        data_dir, prefix, version, logger_file)
    os.system(sync_cmd)
