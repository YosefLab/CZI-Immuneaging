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
sc.settings.verbosity = 1   # verbosity: errors (0), warnings (1), info (2), hints (3)

timestamp = get_current_time()
sandbox_mode = configs["sandbox_mode"] == "True"

# apply the aws credentials to allow access though aws cli; make sure the user is authorized to run in non-sandbox mode if applicable
s3_dict = set_access_keys(configs["s3_access_file"], return_dict = True)
assert sandbox_mode or hashlib.md5(bytes(s3_dict["AWS_SECRET_ACCESS_KEY"], 'utf-8')).hexdigest() in AUTHORIZED_EXECUTERS, "You are not authorized to run this script in a non sandbox mode; please set sandbox_mode to True"
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

def download_aligned_lib_artifact(file_name: str, data_dir: str) -> str:
    sync_cmd = 'aws s3 sync --no-progress s3://immuneaging/aligned_libraries/{}/{}_{}_{}_{}/ {} --exclude "*" --include {}'.format(
        configs["aligned_library_configs_version"], configs["donor"], configs["seq_run"], configs["library_type"],
        configs["library_id"], data_dir, file_name
    )
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
    file_path = os.path.join(data_dir, file_name)
    if not os.path.isfile(file_path):
        msg = "Failed to download file {} from S3.".format(file_name)
        logger.add_to_log(msg, level="error")
        raise ValueError(msg)
    return file_path

def store_lib_alignment_metrics(adata, data_dir):
    logger.add_to_log("Downloading cellranger metrics csv file of aligned library from S3...")
    metrics_csv_file_name = "{}_{}_{}_{}.cellranger.metrics_summary.csv".format(configs["donor"], configs["seq_run"], configs["library_type"], configs["library_id"])
    metrics_csv_file = download_aligned_lib_artifact(metrics_csv_file_name, data_dir)

    logger.add_to_log("Storing alignment metrics from cellranger in adata.uns...")
    lib_metrics = pd.read_csv(metrics_csv_file)
    adata.uns["lib_metrics"] = {}
    for metric in lib_metrics.columns:
        adata.uns["lib_metrics"][metric] = lib_metrics[metric][0] # lib_metrics[metric] is a pandas series with a single row 

def flush_logs_and_upload():
    for i in summary:
        logger.add_to_log(i)
    logging.shutdown()
    if not sandbox_mode:
        sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/processed_libraries/{}/{}/ --exclude "*" --include {}'.format(
            data_dir, prefix, version, logger_file)
        os.system(sync_cmd)
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

if configs["library_type"] == "GEX":
    logger.add_to_log("Downloading h5ad file of aligned library from S3...")
    aligned_h5ad_file_name = "{}_{}.{}.{}.h5ad".format(configs["donor"], configs["seq_run"], configs["library_id"], configs["aligned_library_configs_version"])
    aligned_h5ad_file = download_aligned_lib_artifact(aligned_h5ad_file_name, data_dir)

    logger.add_to_log("Reading aligned h5ad file...")
    adata = sc.read_h5ad(aligned_h5ad_file)
    summary.append("Started with a total of {} cells and {} genes.".format(adata.n_obs, adata.n_vars))

    store_lib_alignment_metrics(adata, data_dir)
    
    logger.add_to_log("Dropping out if this is a known poor quality library...")
    logger.add_to_log("Downloading list of poor_quality libraries from AWS...")
    cell_filtering_aws_dir = "s3://immuneaging/cell_filtering/"
    poor_quality_libs_df = read_csv_from_aws(data_dir, cell_filtering_aws_dir, "poor_quality_libs.csv", logger)
    if configs["donor"] + "_" + configs["library_id"] in poor_quality_libs_df.columns:
        logger.add_to_log("This is a known poor quality library: {}. Exiting...".format(configs["library_id"]), level="warning")
        flush_logs_and_upload()
        sys.exit()

    # Add the library ID to the cell barcode name. This will enable two things:
    # a) distinguish between different cells with the same barcode when integrating different libraries (of the same type) that were used to collect the same samples
    # b) successfully merge ir data with gex data during sample processing, the merging is contingent on the obs_names being exactly the same
    # We need to do this before the next step where we filter out non-immune cells since the black-listed cell barcodes include the library-id.
    logger.add_to_log("Adding the library ID to the cell barcode name...")
    adata.obs_names = adata.obs_names + "_" + configs["library_id"]

    logger.add_to_log("Filtering out non-immune cells...")
    tissues = ["BAL", "BLO", "ILN", "JEJEPI", "JEJLP", "LIV", "MLN", "SKN", "TLN"]
    for tissue in tissues:
        logger.add_to_log("Downloading list of non-immune cells for tissue {}...".format(tissue))
        non_immune_cells_df = read_csv_from_aws(data_dir, cell_filtering_aws_dir, "{}_blacklist.csv".format(tissue), logger)
        # the blacklist for some tissues includes whether the cell barcode should be discarded entirely
        # or whether it should merely be discarded from analysis but still kept around (mast cell are an
        # example). In such cases, only remove the cells that are marked "exclude from the dataset".
        if "Exclude from dataset" in non_immune_cells_df.columns:
            exclude_cells_barcodes = non_immune_cells_df[non_immune_cells_df["Exclude from dataset"] == "Yes"]["cell_barcode"]
            # also add the "Exclude from Aging analysis" as an obs column to adata
            exclude_key = "Exclude from Aging analysis"
            inter = np.intersect1d(adata.obs_names, non_immune_cells_df["cell_barcode"])
            adata.obs[exclude_key] = non_immune_cells_df.set_index("cell_barcode").loc[inter][exclude_key]
        else:
            exclude_cells_barcodes = non_immune_cells_df["cell_barcode"]
        n_obs_before = adata.n_obs
        n_obs_after = np.sum(~adata.obs_names.isin(exclude_cells_barcodes.values))
        n_cells_filtered = n_obs_before-adata.n_obs
        percent_removed = 100*n_cells_filtered/n_obs_before
        level = "warning" if percent_removed > 20 else "info"
        logger.add_to_log("Filtered out {} non-immune cells for tissue {}, percent_removed: {}".format(n_cells_filtered, tissue, percent_removed), level=level)

    # move protein/hto data out of adata.X into adata.obsm/obs
    hto_tag = configs["donor"]+"-"
    cell_hashing = [i for i in adata.var_names[np.where(adata.var_names.str.startswith(hto_tag))]]
    if len(cell_hashing) >= 1:
        logger.add_to_log("Moving hto data out of adata.X into adata.obs...")
        adata.obs[cell_hashing] = adata[:,cell_hashing].X.toarray()
        adata = adata[:, ~adata.var_names.str.startswith(hto_tag)].copy()
        # rename TLN to LLN
        new_cell_hashing = {c: c.replace("TLN", "LLN") for c in cell_hashing}
        adata.obs.rename(columns=new_cell_hashing, errors='raise', inplace=True)
        cell_hashing = list(new_cell_hashing.values())
        # The 694B donor had some of its libraries aligned with HTO tags that end in -1
        # and some with HTO tags that end in -X where X corresponds to each unique sample
        # ID (in short, this was caused by a failure that was fixed mid-way through lib
        # alignment and we didn't want to re-align all the libraries after the fix)
        # We patch that here by renaming the HTO tags to the new ones that match the IA
        # sample spreadsheet.
        # Update: Owner 759B needs the same treatment...
        def old_name_to_new_name(old_name: str) -> str:
            # replace BLD with BLO since that is something else that we renamed after having aligned
            # the libraries
            old_name = old_name.replace("BLD", "BLO")
            if old_name == "694B-BLO-1":
                if configs["library_id"] in ["CZI-IA11512685", "CZI-IA11512686", "CZI-IA11512687"]:
                    return "694B-BLO-203"
                else:
                    return "694B-BLO-210"
            elif old_name == "694B-BMA-1":
                if configs["library_id"] in ["CZI-IA11512688", "CZI-IA11512689", "CZI-IA11485873"]:
                    return "694B-BMA-204"
                else:
                    return "694B-BMA-211"
            elif old_name == "694B-SPL-1":
                if configs["library_id"] in ["CZI-IA11512688", "CZI-IA11512689", "CZI-IA11485873"]:
                    return "694B-SPL-205"
                else:
                    return "694B-SPL-212"
            elif old_name == "694B-JEJEPI-1":
                if configs["library_id"] in ["CZI-IA11512685", "CZI-IA11512686", "CZI-IA11512687"]:
                    return "694B-JEJEPI-207"
                else:
                    return "694B-JEJEPI-213"
            elif old_name == "694B-JEJLP-1":
                if configs["library_id"] in ["CZI-IA11512685", "CZI-IA11512686", "CZI-IA11512687"]:
                    return "694B-JEJLP-208"
                else:
                    return "694B-JEJLP-214"
            elif old_name == "694B-SKN-1":
                raise ValueError("This is unexpected!")
            # 759B
            elif (old_name == "759B-MLN-1") and (configs["library_id"] in ["CZI-IA12953908", "CZI-IA12953909", "CZI-IA12953910", "CZI-IA12953911"]):
                return "759B-MLN-263"
            elif (old_name == "759B-LLN-1") and (configs["library_id"] in ["CZI-IA12953908", "CZI-IA12953909", "CZI-IA12953910", "CZI-IA12953911"]):
                return "759B-LLN-264"
            elif (old_name == "759B-SPL-1") and (configs["library_id"] in ["CZI-IA12953908", "CZI-IA12953909", "CZI-IA12953910", "CZI-IA12953911"]):
                return "759B-SPL-265"
            elif (old_name == "759B-BLO-1") and (configs["library_id"] in ["CZI-IA12953908", "CZI-IA12953909", "CZI-IA12953910", "CZI-IA12953911"]):
                return "759B-BLO-266"
            elif (old_name == "759B-BMA-1") and (configs["library_id"] in ["CZI-IA12953908", "CZI-IA12953909", "CZI-IA12953910", "CZI-IA12953911"]):
                return "759B-BMA-267"
            elif (old_name == "759B-JEJLP-1") and (configs["library_id"] in ["CZI-IA12953912"]):
                return "759B-JEJLP-268"
            elif (old_name == "759B-JEJEPI-1") and (configs["library_id"] in ["CZI-IA12953912"]):
                return "759B-JEJEPI-269"
            elif (old_name == "759B-SKN-1") and (configs["library_id"] in ["CZI-IA12953914"]):
                return "759B-SKN-270"
            else:
                return old_name
        new_cell_hashing = {c: old_name_to_new_name(c) for c in cell_hashing}
        adata.obs.rename(columns=new_cell_hashing, errors='raise', inplace=True)
        cell_hashing = list(new_cell_hashing.values())

    if "Antibody Capture" in np.unique(adata.var.feature_types):
        logger.add_to_log("Moving protein data out of adata.X into adata.obsm...")
        protein = adata[:, adata.var["feature_types"] == "Antibody Capture"].copy()
        protein_df = protein.to_df()
        # switch the protein names to their internal names defined in the protein panels (in the Google Spreadsheet)
        protein_df.columns = get_internal_protein_names(protein_df)
        # save control and non-control proteins in different obsm structures
        is_ctrl_protein = np.array([i.endswith("Ctrl") for i in protein_df.columns])
        protein_expression_obsm_key = "protein_expression"
        protein_expression_ctrl_obsm_key = "protein_expression_Ctrl"
        adata.obsm[protein_expression_ctrl_obsm_key] = protein_df[protein_df.columns[is_ctrl_protein]].copy()
        adata.obsm[protein_expression_obsm_key] = protein_df[protein_df.columns[np.logical_not(is_ctrl_protein)]].copy()
        adata = adata[:, adata.var["feature_types"] != "Antibody Capture"]

    logger.add_to_log("Applying basic filters...")
    n_cells_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=configs["filter_cells_min_genes"])
    logger.add_to_log("Filtered out {} cells that have less than {} genes expressed.".format(n_cells_before-adata.n_obs, configs["filter_cells_min_genes"]))
    if "filter_cells_min_umi" in configs:
        n_cells_before = adata.n_obs
        adata = adata[adata.X.sum(axis=-1) > configs["filter_cells_min_umi"]].copy()
        logger.add_to_log("Filtered out {} cells that have less than {} total umi's.".format(n_cells_before-adata.n_obs, configs["filter_cells_min_umi"]))
    n_genes_before = adata.n_vars
    gene_subset, _ = sc.pp.filter_genes(adata, min_cells=configs["filter_genes_min_cells"], inplace=False)
    extend_removed_features_df(adata, "removed_genes", adata[:,~gene_subset].copy().to_df())
    adata = adata[:,gene_subset].copy()
    logger.add_to_log("Filtered out {} genes that are detected in less than {} cells.".format(n_genes_before-adata.n_vars, configs["filter_genes_min_cells"]))

    if adata.n_obs == 0:
        logger.add_to_log("No cells left after basic filtering steps. Exiting...", level="error")
        flush_logs_and_upload()
        sys.exit()

    adata.var['mt'] = adata.var_names.str.startswith('MT-') # mitochondrial genes
    adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL")) # ribosomal genes
    hb_genes = list(adata.var_names[adata.var_names.str.contains(("^HB[^(P)]"))]) + ['ALAS2', 'EPOR']
    adata.var['hb'] = [i in hb_genes for i in adata.var_names]
    hsp_genes = ['HSPA6', 'HSPA1A', 'HSPA1B', 'HSPB1', 'HSPH1']
    adata.var['hsp'] = [i in hsp_genes for i in adata.var_names]
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo', 'hb', 'hsp'], percent_top=None, log1p=False, inplace=True)
    n_cells_before = adata.n_obs
    adata = adata[adata.obs['pct_counts_mt'] <= configs["filter_cells_max_pct_counts_mt"], :].copy()
    logger.add_to_log("Filtered out {} cells with more than {}\% counts coming from mitochondrial genes.".format(n_cells_before-adata.n_obs, configs["filter_cells_max_pct_counts_mt"]))
    n_cells_before = adata.n_obs
    adata = adata[adata.obs['pct_counts_ribo'] >= configs["filter_cells_min_pct_counts_ribo"], :].copy()
    logger.add_to_log("Filtered out {} cells with less than {}\% counts coming from ribosomal genes.".format(n_cells_before-adata.n_obs, configs["filter_cells_min_pct_counts_ribo"]))

    if adata.n_obs == 0:
        logger.add_to_log("No cells left after filtering. Exiting...", level="error")
        flush_logs_and_upload()
        sys.exit()

    genes_to_exclude = set()
    if configs["genes_to_exclude"] != "None":
        for gene in configs["genes_to_exclude"].split(','):
            gene_names = set(adata.var_names[adata.var_names.str.startswith(gene)])
            genes_to_exclude.update(gene_names)
    if configs["exclude_mito_genes"] == "True":
            gene_names = set(adata.var_names[adata.var_names.str.startswith('MT-')])
            genes_to_exclude.update(gene_names)
    n_genes_before = adata.n_vars
    genes_to_exclude_idx = adata.var_names.isin(genes_to_exclude)
    # add the genes to exclude to the removed_genes obsm df
    exclude_df = adata[:, genes_to_exclude_idx].copy().to_df()
    extend_removed_features_df(adata, "removed_genes", exclude_df)
    adata = adata[:, ~genes_to_exclude_idx].copy()
    logger.add_to_log("Filtered out the following {} genes: {}".format(n_genes_before-adata.n_vars, ", ".join(genes_to_exclude)))

    if len(cell_hashing) > 1:
        logger.add_to_log("Demultiplexing is needed; using hashsolo...")
        hashsolo_priors = [float(i) for i in configs["hashsolo_priors"].split(',')]
        sc.external.pp.hashsolo(adata, cell_hashing_columns = cell_hashing, priors = hashsolo_priors, inplace = True,
            number_of_noise_barcodes = len(cell_hashing)-1)
        num_doublets = sum(adata.obs["Classification"] == "Doublet")
        percent_doublets = 100*num_doublets/adata.n_obs
        level = "error" if percent_doublets > 40 else "info"
        logger.add_to_log("Removing {:.2f}% of the droplets ({} droplets out of {}) called by hashsolo as doublets...".format(percent_doublets, num_doublets, adata.n_obs), level=level)
        adata = adata[adata.obs["Classification"] != "Doublet"].copy()

    summary.append("Final number of cells: {}, final number of genes: {}.".format(adata.n_obs, adata.n_vars))
elif configs["library_type"] == "BCR" or configs["library_type"] == "TCR":
    logger.add_to_log("Downloading aligned library from S3...")
    aligned_csv_file_name = "{}_{}_{}_{}.cellranger.filtered_contig_annotations.{}.csv".format(
        configs["donor"],
        configs["seq_run"],
        configs["library_type"],
        configs["library_id"],
        configs["aligned_library_configs_version"]
    )
    aligned_csv_file = download_aligned_lib_artifact(aligned_csv_file_name, data_dir)

    logger.add_to_log("Reading aligned csv file...")
    adata = ir.io.read_10x_vdj(aligned_csv_file)
    summary.append("Started with a total of {} cells.".format(adata.n_obs))

    store_lib_alignment_metrics(adata, data_dir)

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
    # (Note: We're not filtering our orphan-chain cells. They can be matched to clonotypes
    # on a single chain only, by using receptor_arms=”any” when running scirpy.tl.define_clonotypes()
    # - see scirpy docs)
    n_cells_before = adata.n_obs
    n_multichain_cells = sum(adata.obs["multi_chain"] == "True")
    n_multichain_cells_pct = (n_multichain_cells/adata.n_obs) * 100
    n_ambiguous_cells = sum(adata.obs["chain_pairing"] == "ambiguous")
    n_ambiguous_cells_pct = (n_ambiguous_cells/adata.n_obs) * 100
    indices = (adata.obs["multi_chain"] == "False") & (adata.obs["chain_pairing"] != "ambiguous")
    level = "warning" if n_multichain_cells_pct > 5 else "info"
    logger.add_to_log("Removed multichains (individual count and percentage: {}, {:.2f}%).".format(n_multichain_cells, n_multichain_cells_pct), level=level)
    level = "warning" if n_ambiguous_cells_pct > 5 else "info"
    logger.add_to_log("Removed ambiguous cells (individual count and percentage: {}, {:.2f}%).".format(n_ambiguous_cells, n_ambiguous_cells_pct), level=level)
    logger.add_to_log("Original cell count: {}, cell count after the filtering: {}.".format(n_cells_before, adata.n_obs))

    logger.add_to_log("Validating that there are no non-cells, no non-IR cells, and that cells' receptor_type match the lib type.")
    n_no_cells = sum(adata.obs["is_cell"] == "False")
    n_no_ir_cells = sum((adata.obs["chain_pairing"] == "no IR") | (adata.obs["has_ir"] == "False"))    
    n_unexpected_ir_type_cells = sum(adata.obs["receptor_type"] != configs["library_type"]) # BCR or TCR
    if n_no_cells > 0:
        logger.add_to_log("Detected {} barcodes not called as a cell.".format(n_no_cells), level = "error")
        sys.exit()
    elif n_no_ir_cells > 0:
        logger.add_to_log("Detected {} cells with no IR detected.".format(n_no_ir_cells), level = "error")
        sys.exit()
    elif n_unexpected_ir_type_cells > 0:
        logger.add_to_log("Detected {} cells with a different receptor_type than expected. unique receptor types found: {}.".format(
            n_unexpected_ir_type_cells, np.unique(adata.obs["receptor_type"])), level = "error")
        sys.exit()

    # this is for clarity and also to distinguish similarly-named BCR and TCR obs columns
    logger.add_to_log("Pre-pending obs column names with the library type.")
    adata.obs = adata.obs.add_prefix("{}-".format(configs["library_type"]))

    # Add the corresponding GEX library ID to the cell barcode name for the same reasons as we do for GEX (see above)
    logger.add_to_log("Adding the corresponding GEX library ID to the cell barcode name...")
    adata.obs_names = adata.obs_names + "_" + configs["corresponding_gex_lib"]

    summary.append("Final number of cells: {}.".format(adata.n_obs))
else:
    raise ValueError("Unrecognized lib type: {}".format(configs["library_type"]))

logger.add_to_log("Saving h5ad file...")
adata.obs[f'library_pipeline_version_{configs["library_type"]}'] = f"{configs['library_type']}__{configs['library_id']}__{configs['pipeline_version']}"
adata.obs[f'library_code_version__{configs["library_type"]}'] =  f"{configs['library_type']}__{configs['library_id']}__{configs['code_version']}"
write_anndata_with_object_cols(adata, data_dir, h5ad_file)

if not sandbox_mode:
    logger.add_to_log("Uploading h5ad file to S3...")
    sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/processed_libraries/{}/{}/ --exclude "*" --include {}'.format(
        data_dir, prefix, version, h5ad_file)
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

logger.add_to_log("Execution of process_library.py is complete.")

flush_logs_and_upload()