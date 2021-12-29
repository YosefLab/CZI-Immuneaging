import sys
import logging
import os
import json
import scanpy as sc
import numpy as np
import hashlib
import scirpy as ir

process_lib_script = sys.argv[0]
configs_file = sys.argv[1]

#########################################################
###### INITIALIZATIONS AND PREPARATIONS BEGIN HERE ######
#########################################################

with open(configs_file) as f: 
    data = f.read()	

configs = json.loads(data)

sys.path.append(configs["code_path"])

from logger import SimpleLogger
from utils import *

VARIABLE_CONFIG_KEYS = ["data_owner","s3_access_file","code_path","output_destination"] # config changes only to these fields will not initialize a new configs version
sc.settings.verbosity = 3   # verbosity: errors (0), warnings (1), info (2), hints (3)

timestamp = get_current_time()
sandbox_mode = configs["sandbox_mode"] == "True"

# apply the aws credentials to allow access though aws cli; make sure the user is authorized to run in non-sandbox mode if applicable
s3_dict = set_access_keys(configs["s3_access_file"], return_dict = True)
assert sandbox_mode or hashlib.md5(bytes(s3_dict["AWS_SECRET_ACCESS_KEY"], 'utf-8')).hexdigest() in AUTHORIZED_EXECUTERS, "You are not authorized to run this script in a non sanbox mode; please set sandbox_mode to True"
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
    sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/processed_libraries/{}/{}/ --exclude "*" --include {}'.format(
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

if configs["library_type"] == "GEX":
    logger.add_to_log("Downloading h5ad file of aligned library from S3...")
    aligned_h5ad_file = "{}_{}.{}.{}.h5ad".format(configs["donor"], configs["seq_run"],
        configs["library_id"], configs["aligned_library_configs_version"])
    sync_cmd = 'aws s3 sync --no-progress s3://immuneaging/aligned_libraries/{}/{}_{}_{}_{}/ {} --exclude "*" --include {}'.format(
        configs["aligned_library_configs_version"], configs["donor"], configs["seq_run"], configs["library_type"],
        configs["library_id"], data_dir, aligned_h5ad_file)
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

    logger.add_to_log("Reading aligned h5ad file...")
    aligned_h5ad_file = os.path.join(data_dir, aligned_h5ad_file)
    adata = sc.read_h5ad(aligned_h5ad_file)
    summary.append("Started with a total of {} cells and {} genes.".format(adata.n_obs, adata.n_vars))

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
    adata = adata[adata.obs['pct_counts_mt'] <= configs["filter_cells_max_pct_counts_mt"], :].copy()
    logger.add_to_log("Filtered out {} cells with more than {}\% counts coming from mitochondrial genes.".format(n_cells_before-adata.n_obs, configs["filter_cells_max_pct_counts_mt"]))
    n_cells_before = adata.n_obs
    adata = adata[adata.obs['pct_counts_ribo'] >= configs["filter_cells_min_pct_counts_ribo"], :].copy()
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
    adata = adata[:,genes_to_keep].copy()
    logger.add_to_log("Filtered out the following {} genes: {}".format(n_genes_before-adata.n_vars, ", ".join(genes_to_exclude_names)))

    cell_hashing = [i for i in adata.var_names[np.where(adata.var_names.str.startswith(configs["donor"]+"-"))]]
    if len(cell_hashing)>1:
        logger.add_to_log("Demultiplexing is needed; using hashsolo...")
        # add the HTO counts to .obs
        adata.obs[cell_hashing] = adata[:,cell_hashing].X.toarray()
        hashsolo_priors = [float(i) for i in configs["hashsolo_priors"].split(',')]
        sc.external.pp.hashsolo(adata, cell_hashing_columns = cell_hashing, priors = hashsolo_priors, inplace = True,
            number_of_noise_barcodes = len(cell_hashing)-1)
        num_doublets = sum(adata.obs["Classification"] == "Doublet")
        percent_doublets = 100*num_doublets/adata.n_obs
        level = "error" if percent_doublets > 40 else "info"
        logger.add_to_log("Removing {:.2f}% of the droplets ({} droplets out of {}) called by hashsolo as doublets...".format(percent_doublets, num_doublets, adata.n_obs), level=level)
        adata = adata[adata.obs["Classification"] != "Doublet"].copy()
        logger.add_to_log("Adding the library ID to the cell barcode name (will allow to distinguish between differnet cells with the same barcode when integrating differnet libraries that were used to collect the same samples)...")
        adata.obs_names = adata.obs_names + "_" + configs["library_id"]

    summary.append("Final number of cells: {}, final number of genes: {}.".format(adata.n_obs, adata.n_vars))
elif configs["library_type"] == "BCR" or configs["library_type"] == "TCR":
    logger.add_to_log("Downloading aligned library from S3...")
    aligned_csv_file = "{}_{}_{}_{}.cellranger.filtered_contig_annotations.csv".format(
        configs["donor"],
        configs["seq_run"],
        configs["library_type"],
        configs["library_id"]
    )
    sync_cmd = 'aws s3 sync --no-progress s3://immuneaging/aligned_libraries/{}/{}_{}_{}_{}/ {} --exclude "*" --include {}'.format(
        configs["aligned_library_configs_version"], configs["donor"], configs["seq_run"], configs["library_type"],
        configs["library_id"], data_dir, aligned_csv_file)
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

    logger.add_to_log("Reading aligned csv file...")
    aligned_csv_file = os.path.join(data_dir, aligned_csv_file)
    adata = ir.io.read_10x_vdj(aligned_csv_file)
    summary.append("Started with a total of {} cells.".format(adata.n_obs))

    logger.add_to_log("Filtering out cell calls that are marked low confidence by cellranger.")
    # for more info about what this and other cellranger vdj output fields mean, see https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/output/annotation
    n_low_confidence_cells = sum(adata.obs["high_confidence"] == "False")
    n_low_confidence_cells_pct = (n_low_confidence_cells/adata.n_obs) * 100
    adata = adata[adata.obs["high_confidence"] == "True"].copy()
    level = "warning" if n_low_confidence_cells_pct > 5 else "info"
    logger.add_to_log("Removed {} ({:.2f}%) cells called with low confidence by cellranger.".format(n_low_confidence_cells, n_low_confidence_cells_pct), level=level)

    logger.add_to_log("Applying basic filters and chain QC using scirpy...")
    ir.tl.chain_qc(adata)
    # filter out cells that are multichain or ambiguous as these likely represent doublets
    # plus, filter our "orphan chain" cells
    # (Note: Orphan chain cells can also be matched to clonotypes on a single chain only, by
    # using receptor_arms=”any” when running scirpy.tl.define_clonotypes(). TODO should we do this?)
    n_cells_before = adata.n_obs
    n_multichain_cells = sum(adata.obs["multi_chain"] == "True")
    n_multichain_cells_pct = (n_multichain_cells/adata.n_obs) * 100
    n_ambiguous_cells = sum(adata.obs["chain_pairing"] == "ambiguous")
    n_ambiguous_cells_pct = (n_ambiguous_cells/adata.n_obs) * 100
    n_orphan_cells = sum(adata.obs["chain_pairing"] == "orphan VJ" | adata.obs["chain_pairing"] == "orphan VDJ")
    n_orphan_cells_pct = (n_orphan_cells/adata.n_obs) * 100
    # TODO should we also filter out "single pair" and "extra V(D)J" cells?
    # TODO should we select "cells with a single pair of productive αβ TCR chains" (same as celltypist studies)
    indices = (adata.obs["multi_chain"] == "False") & (adata.obs["chain_pairing"] != "ambiguous") & (adata.obs["chain_pairing"] != "orphan VJ") & (adata.obs["chain_pairing"] != "orphan VDJ")
    adata = adata[indices].copy()
    level = "warning" if n_multichain_cells_pct > 5 else "info"
    logger.add_to_log("Removed multichains (individual count and percentage: {}, {:.2f}%).".format(n_multichain_cells, n_multichain_cells_pct), level=level)
    level = "warning" if n_ambiguous_cells_pct > 5 else "info"
    logger.add_to_log("Removed ambiguous cells (individual count and percentage: {}, {:.2f}%).".format(n_ambiguous_cells, n_ambiguous_cells_pct), level=level)
    level = "warning" if n_orphan_cells_pct > 50 else "info"
    logger.add_to_log("Removed orphan V(D)J cells (individual count and percentage: {}, {:.2f}%).".format(n_orphan_cells, n_orphan_cells_pct), level=level)
    logger.add_to_log("Original cell count: {}, cell count after all the filtering: {}.".format(n_cells_before, adata.n_obs))

    logger.add_to_log("Validating that there are no non-cells, no non-IR cells, and that cells' receptor_type match the lib type.")
    n_no_cells = sum(adata.obs["is_cell"] == "False")
    n_no_ir_cells = sum(adata.obs["chain_pairing"] == "no IR" | adata.obs["has_ir"] == "False")
    n_unexpected_ir_type_cells = sum(adata.obs["receptor_type"] != configs["library_type"]) # BCR or TCR
    if n_no_cells > 0:
        logger.add_to_log("Detected {} barcodes not called as a cell.".format(n_no_cells), level = "error")
        sys.exit()
    elif n_no_ir_cells > 0:
        logger.add_to_log("Detected {} cells with no IR detected.".format(n_no_ir_cells), level = "error")
        sys.exit()
    elif n_unexpected_ir_type_cells > 0:
        logger.add_to_log("Detected {} cells with a different receptor_type than expected. unique receptor types found: {}.".format(n_unexpected_ir_type_cells, np.unique(adata.obs["receptor_type"])), level = "error")
        sys.exit()

    logger.add_to_log("Pre-pending obs column names with the library type.")
    adata.obs = adata.obs.add_prefix("{}-".format(configs["library_type"]))

    summary.append("Final number of cells: {}.".format(adata.n_obs))
else:
    raise ValueError("Unrecognized lib type: {}".format(configs["library_type"]))

logger.add_to_log("Saving h5ad file...")
adata.write(os.path.join(data_dir,h5ad_file), compression="lzf")

if not sandbox_mode:
    logger.add_to_log("Uploading h5ad file to S3...")
    sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/processed_libraries/{}/{}/ --exclude "*" --include {}'.format(
        data_dir, prefix, version, h5ad_file)
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

logger.add_to_log("Execution of process_library.py is complete.")

for i in summary:
    logger.add_to_log(i)

logging.shutdown()
if not sandbox_mode:
    sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/processed_libraries/{}/{}/ --exclude "*" --include {}'.format(
        data_dir, prefix, version, logger_file)
    os.system(sync_cmd)
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
