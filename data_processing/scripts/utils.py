import json
import os
import re
import time
import glob
import warnings
import pandas as pd
import numpy as np
import scvi
import zipfile
import anndata
from anndata import AnnData
from math import floor
import csv
from typing import Type, List, NamedTuple, Optional
import traceback
from datetime import datetime
from logger import BaseLogger
import scanpy as sc
# import celltypist
import logging

logging.getLogger('numba').setLevel(logging.WARNING)

AUTHORIZED_EXECUTERS = ["b750bd0287811e901c88dc328187e25f", "1c75133ab6a1fc3ed9233d3fe40b3d73"] # md5 checksums of the AWS_SECRET_ACCESS_KEY value of those that are authorized to upload outputs of processing scripts to the server; note that individuals with upload permission to aws can bypass that by changing the code - this is just designed to alert users that they should only use sandbox mode.

class CELLRANGER_METRICS_NT(NamedTuple):
    MEDIAN_GENES_PER_CELL: str = "Median Genes per Cell"
    MEDIAN_UMI_COUNTS_PER_CELL: str = "Median UMI Counts per Cell"
    SEQUENCING_SATURATION: str = "Sequencing Saturation"

CELLRANGER_METRICS = CELLRANGER_METRICS_NT()

# these are formatted strings that we use to log QC stats and use them to parse those lines back
# if you edit these make sure to update all call sites that use them
QC_STRING_DOUBLETS = "Removed {} estimated doublets (percent removed: {:.2f}%); {} droplets remained."
QC_STRING_AMBIENT_RNA = "Removed {} cells (percent removed: {:.2f}%) with total decontaminated counts below filter_decontaminated_cells_min_genes={}"
QC_STRING_VDJ = "Removed {} vdj genes (percent removed: {:.2f}%); {} genes remained."
QC_STRING_RBC = "Removed {} red blood cells (percent removed: {:.2f}%); {} droplets remained."
QC_STRING_COUNTS = "Final number of cells: {}, final number of genes: {}."
QC_STRING_START_TIME = "Starting time: {}"

def init_scvi_settings():
    # This does two things:
    # 1. Makes the logger look good in a log file
    # 2. Changes a bit how torch pins memory when copying to GPU, which allows you to more easily run models in parallel with an estimated 1-5% time hit
    scvi.settings.reset_logging_handler()
    scvi.settings.dl_pin_memory_gpu_training = False

init_scvi_settings()

TIME_FORMAT = "%H:%M, %m-%d-%Y"
DATE_FORMAT = "%m/%d/%Y"

def get_current_time():
    return time.strftime(TIME_FORMAT)

def get_date_from_time(t: str):
    dt = datetime.strptime(t, TIME_FORMAT)
    return dt.strftime(DATE_FORMAT)

def set_access_keys(filepath, return_dict = False):
	"""
	Sets the user's access keys to the AWS S3 bucket.

	Assumes this file includes only two uncommented lines with keys:
	export AWS_ACCESS_KEY_ID=<key>
	export AWS_SECRET_ACCESS_KEY=<key>

	If return_dict == True then only returns the dictionary.
	"""
	keys = {}
	with open(filepath) as fp:
		for i in fp:
			if len(i.rstrip()) and i[0]!="#":
				cmd = i.rstrip().split(" ")
				assert(cmd[0] == "export")
				cmd.pop(0)
				for i in range(len(cmd)-1,-1,-1):
					if len(cmd[i]) == 0:
						cmd.pop(i)
				pair = "".join(cmd).split("=")
				keys[pair[0]] = pair[1]
	if return_dict:
		return(keys)
	for k in keys:
		os.environ[k] = keys[k]

def load_configs(filename):
    with open(filename) as f: 
        data = f.read()	
    configs = json.loads(data)
    return configs

def get_configs_status(configs, s3_path, configs_file_prefix, variable_config_keys, data_dir):
	ls_cmd = 'aws s3 ls {} --recursive'.format(s3_path)
	files = os.popen(ls_cmd).read()
	latest_configs_file = None
	latest_version_num = 0
	is_new_version = False
	for f in files.rstrip().split('\n'):
		j = f.split('/')[-1].split('.')[0:-1]
		if ".".join(j[0:-1]) == configs_file_prefix:
			v = int(j[-1][1:])
			if v > latest_version_num:
				latest_version_num = v
				latest_configs_file = "/".join(f.split('/')[-2:])
	if latest_configs_file is None:
		version = 1
	else:
		cp_cmd = 'aws s3 cp {}/{} {}'.format(s3_path, latest_configs_file, data_dir)
		os.system(cp_cmd)
		configs_invariant = {i:configs[i] for i in configs if i not in variable_config_keys}
		with open(os.path.join(data_dir, latest_configs_file.split('/')[-1])) as f:
			latest_configs = json.loads(f.read().rstrip())
			latest_configs_invariant = {i:latest_configs[i] for i in latest_configs if i not in variable_config_keys}
			if configs_invariant == latest_configs_invariant:
				version = latest_version_num
			else:
				version = latest_version_num + 1
	if version > latest_version_num:
		is_new_version = True
	return [is_new_version,"v"+str(version)]

def get_latest_lib_version(s3_access_file: str, s3_lib_path: str, lib_folder_name: Optional[str] = None):
    set_access_keys(s3_access_file)
    latest_version = -1
    ls_cmd = 'aws s3 ls {} --recursive'.format(s3_lib_path)
    ls  = os.popen(ls_cmd).read()
    if len(ls) != 0:
        # search for patterns of /vX/. If there is a match, the regex group
        # one is "X" (X can be any integer >0)
        pattern = "/v(\d+)/"
        filenames = ls.split("\n")
        for filename in filenames:
            if lib_folder_name is not None and lib_folder_name not in filename:
                continue
            m = re.search(pattern, filename)
            if bool(m):
                version = int(m[1])
                if latest_version < version:
                    latest_version = version
    if latest_version == -1:
        print(f"Could not find latest version for lib. s3_lib_path: {s3_lib_path}")
    return "v" + str(latest_version)

def get_configs_version_alignment(configs, data_dir, configs_dir_remote, configs_file_remote_prefix, variable_config_keys):
	# This function checks if the configs file is using configs that were already used (while disregarding fields variable_config_keys) and are therefore documented on S3.
	# If this is the first time these configs are used then it creates a new configs version and uploads the new configs to S3.
	# Finally, the version of the configs is returned; to be used for stamping results files with the configs version.
	configs_dir_local = os.path.join(data_dir,"configs")
	os.system("mkdir -p " + configs_dir_local)
	cmd = 'aws s3 sync {0} {1} --exact-timestamps'.format(configs_dir_remote, configs_dir_local)
	os.system(cmd)
	version = None
	max_version = 0
	configs_invariant = {i:configs[i] for i in configs if i not in variable_config_keys}
	for configs_file in glob.glob(os.path.join(configs_dir_local,"*.txt")):
		with open(configs_file) as f:
			f_configs = json.loads(f.read().rstrip())
		if configs_invariant == f_configs:
			version = configs_file.split('/')[-1][0:-4].split('.')[1]
			break
		version_num = int(configs_file.split('/')[-1][0:-4].split('.')[1][1:])
		if (version_num > max_version):
			max_version = version_num
	if version is None:
		version = "v{0}".format(max_version+1)
		configs_file = os.path.join(configs_dir_local,"{0}{1}.txt".format(configs_file_remote_prefix, version))
		with open(configs_file, 'w') as f:
			json.dump(configs_invariant, f)
		cmd = 'aws s3 sync {0} {1} --exclude "*" --include {2}'.format(configs_dir_local,configs_dir_remote,configs_file)
		#print("Uploading configs_file to S3 : {0}".format(os.popen(cmd).read()))
		os.system(cmd)
	return version

def zipdir(path, ziph):
    for root, dirs, files in os.walk(path):
        for file in files:
            ziph.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), 
                os.path.join(path, '..')))

def read_immune_aging_sheet(sheet, output_fn=None, sheet_name=None, quiet=False):
    url = "https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/gviz/tq?tqx=out:csv&sheet={}".format(
        sheet
    )

    # testing url of bad sample sheet
    # url = "https://docs.google.com/spreadsheets/d/1YO1HLGLnO3PPUiK1vKZd52yoCInwpLl60zoi4zxkOrE/gviz/tq?tqx=out:csv&sheet={}".format(
    #     sheet
    # )

    try:
        import gdown
    except ImportError as e:
        raise ImportError(
            "gdown is not installed. Please install gdown via: pip install gdown"
        )

    with warnings.catch_warnings(record=True) as w:
        warnings.filterwarnings("error")
        # output_fn = gdown.download(url, output_fn, quiet=quiet, fuzzy=True)
        output_fn = gdown.download(url, output_fn, quiet=quiet)

        if len(w) == 1:
            msg = w[0]
            warnings.showwarning(
                msg.message, msg.category, msg.filename, msg.lineno, msg.line
            )

    data = pd.read_csv(output_fn)
    os.remove(output_fn)
    return data

def draw_separator_line():
    try:
        width = os.get_terminal_size().columns / 5
        print(" " * floor(width)  + "\u2014" * 3 * floor(width) + " " * floor(width) + "\n")
    except:
        # we might end up here if we can't get the terminal size. In this case, just draw a line
        # with a hard-coded length
        width = 50
        print("\u2014" * width + "\n")

def write_anndata_with_object_cols(adata: AnnData, data_dir: str, h5ad_file: str) -> None:
    try:
        adata.write(os.path.join(data_dir,h5ad_file), compression="lzf")
    except:
        # There can be some BCR-/TCR- columns that have dtype "object" due to being all NaN, thus causing
        # the write to fail. We replace them with 'nan'. Note this isn't ideal, however, since some of those
        # columns can be non-string types (e.g. they can be integer counts), but is something we can handle
        # in future processing layers.
        obj_cols = adata.obs.select_dtypes(include='object').columns
        # Also call .astype("str") since there can be other values than NaN in the column that contribute to
        # the "object" type
        adata.obs.loc[:, obj_cols] = adata.obs.loc[:, obj_cols].fillna('nan').astype("str")
        adata.write(os.path.join(data_dir,h5ad_file), compression="lzf")
        
def run_model(
        adata: AnnData,
        configs: dict,
        batch_key: str,
        protein_expression_obsm_key: str,
        model_name: str,
        prefix: str,
        version: str,
        data_dir: str,
        logger: Type[BaseLogger],
        latent_key: str = None,
        max_retry_count: int = 0,
    ):
    """
    Wrapper for `_run_model_impl` that retries the call up to `max_retry_count` times
    """
    n_tries = 0
    while n_tries < max_retry_count + 1:
        try:
            return _run_model_impl(
                adata,
                configs,
                batch_key,
                protein_expression_obsm_key,
                model_name,
                prefix,
                version,
                data_dir,
                logger,
                latent_key
            )
        except Exception as err:
            n_tries += 1
            logger.add_to_log("Execution failed with the following error: {}. Tried {} time(s)...\n{}".format(err, n_tries, traceback.format_exc()))
            if n_tries == max_retry_count + 1:
                raise
    
def _run_model_impl(
        adata: AnnData,
        configs: dict,
        batch_key: str,
        protein_expression_obsm_key: str,
        model_name: str,
        prefix: str,
        version: str,
        data_dir: str,
        logger: Type[BaseLogger],
        latent_key: str = None,
    ):
    """
    Runs scvi or totalvi model depending on the given model_name.

    Parameters
    ----------
    adata
        The anndata object containing the data to train on.
    configs
        The configs dictionary which includes model parameters.
    batch_key
        Name of the column in adata.obs containing batch information.
    protein_expression_obsm_key
        Name of the column in adata.obs containing protein expression information.
    model_name
        One of "scvi" or "totalvi". Indicates which model to run.
    prefix
        String containing sample id and other information, used in file/dir naming.
    version
        String indicating a version number for the execution.
    data_dir
        Path to the local data directory where processing output is saved.
    logger
        Logger object to use when adding logs.
	latent_key
		key to be used for saving the latent representation in adata.obsm.

    Returns
    -------
    A tuple containing the trained model and the name of the zip file where the model is saved.
    """
    assert model_name in ["scvi", "totalvi"]
    logger.add_to_log("Setting up {}...".format(model_name))
    scvi.data.setup_anndata(adata, batch_key=batch_key, protein_expression_obsm_key=protein_expression_obsm_key)
    empirical_protein_background_prior = None if "empirical_protein_background_prior" not in configs else configs["empirical_protein_background_prior"] == "True"
    model_params_keys = ["use_layer_norm", "use_batch_norm", "n_latent"]
    model_params = dict()
    for i in model_params_keys:
        if i in configs:
            model_params[i] = configs[i]
    logger.add_to_log("Training {} model with batch key {}...".format(model_name, batch_key))
    model = scvi.model.SCVI(adata, **model_params) if model_name=="scvi" else scvi.model.TOTALVI(adata, \
        empirical_protein_background_prior = empirical_protein_background_prior, **model_params)
    train_params = dict()
    train_params_keys = ["lr","early_stopping","train_size","early_stopping_patience","batch_size","limit_train_batches"]
    for i in train_params_keys:
        if i in configs:
            train_params[i] = configs[i]
    max_epochs_config_key = "scvi_max_epochs" if model_name=="scvi" else "totalvi_max_epochs"
    if max_epochs_config_key in configs:
        train_params["max_epochs"] = configs[max_epochs_config_key]
    elif model_name == "totalvi":
        # totalvi has a bug where it fails to heuristically determine max_epochs
        # so for now we replicate that code below (https://github.com/scverse/scvi-tools/issues/1638)
        train_params["max_epochs"] = np.min([round((20000 / adata.n_obs) * 400), 400])
    model.train(**train_params)
    if "elbo_train" in model.history_:
        logger.add_to_log("Number of {} training epochs: {}...".format(model_name, len(model.history_["elbo_train"])))
    logger.add_to_log("Saving {} latent representation...".format(model_name))
    latent = model.get_latent_representation()
    if latent_key is None:
        latent_key = "X_scVI" if model_name=="scvi" else "X_totalVI"
    adata.obsm[latent_key] = latent
    model_file = "{}.{}.{}_model_batch_key_{}.zip".format(prefix, version, model_name, batch_key)
    logger.add_to_log("Saving the model into {}...".format(model_file))
    model_file_path = os.path.join(data_dir, model_file)
    model_dir_path = os.path.join(data_dir,"{}.{}_model_batch_key_{}/".format(prefix, model_name, batch_key))
    if os.path.isdir(model_dir_path):
        os.system("rm -r " + model_dir_path)
    model.save(model_dir_path)
    # save the data used for fitting the model; this is useful for applying reference-based integration on query data later on (based on the current model and data).
    logger.add_to_log("Saving the data used for fitting the model...")
    os.path.join(data_dir, model_file)
    data_file = "{}.{}.{}_model_batch_key_{}.data.h5ad".format(prefix, version, model_name, batch_key)
    adata_copy = adata.copy()
    write_anndata_with_object_cols(adata_copy, model_dir_path, data_file)
    # zip the dir with all the model outputs
    zipf = zipfile.ZipFile(model_file_path, 'w', zipfile.ZIP_DEFLATED)
    zipdir(model_dir_path, zipf)
    zipf.close()
    return model, model_file

def aws_sync(source: str, target: str, include: str, logger: Type[BaseLogger], do_log: bool = True) -> None:
    sync_cmd = 'aws s3 sync --no-progress {} {} --exclude "*" --include {}'.format(source, target, include)
    if do_log:
        logger.add_to_log("syncing {}...".format(include))
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    aws_response = os.popen(sync_cmd).read()
    if do_log:
        logger.add_to_log("aws response: {}\n".format(aws_response))

def filter_vdj_genes(rna: AnnData, aws_file_path: str, data_dir: str, logger: Type[BaseLogger]) -> AnnData:
    file_path_components = aws_file_path.split("/")
    file_name = file_path_components[-1]
    aws_dir_path = "/".join(file_path_components[:-1])
    aws_sync(aws_dir_path, data_dir, file_name, logger)
    local_file_path = os.path.join(data_dir, file_name)
    with open(local_file_path) as csvfile:
        reader = csv.reader(csvfile)
        genes = [row[0] for row in reader]
    n_var_before = rna.n_vars
    rna = rna[:, ~rna.var.index.isin(genes)].copy()
    percent_removed = 100*(n_var_before-rna.n_vars)/n_var_before
    level = "warning" if percent_removed > 50 else "info"
    logger.add_to_log(QC_STRING_VDJ.format(n_var_before-rna.n_vars, percent_removed, rna.n_vars), level=level)
    return rna

# gets dataframe of CITE data extracted from the anndata and converts the protein names to internal names (based on information from the protein panels in the google spreadsheet)
def get_internal_protein_names(df):
    panel_sizes = []
    protein_panels = {}
    protein_panel_columns = np.array(["id", "name", "read", "pattern", "sequence", "feature_type", "internal_name"])
    # find the protein panel that describes the proteins in df
    i = 1
    while True:
        protein_panel_i_name = "Protein panel {}".format(str(i))
        protein_panel_i = read_immune_aging_sheet(protein_panel_i_name)
        if np.sum(np.in1d(protein_panel_columns, protein_panel_i.columns)) != len(protein_panel_columns):
            break
        protein_panels[protein_panel_i_name] = protein_panel_i
        panel_sizes.append(protein_panel_i.shape[0])
        i += 1
    assert len(np.unique(panel_sizes)) == len(panel_sizes) # make sure we can identify the panel solely based on its size
    assert np.sum(df.shape[1] == np.array(panel_sizes)) == 1 # make sure that we can identify the protein panel used
    for protein_panel in protein_panels.values():
        if df.shape[1] == protein_panel.shape[0]:
            return protein_panel["internal_name"]

def dir_and_files_exist(dir_name: str, file_names: List[str]) -> bool:
    exists = os.path.isdir(dir_name)
    if exists:
        for fn in file_names:
            if not os.path.isfile(fn):
                exists = False
                break
    return exists

def is_immune_type(df: pd.DataFrame) -> pd.DataFrame:
    known_immune_types = ["ILC", "T cells", "B cells", "monocytes", "Monocytes", "Macrophages", "macrophages", "NK cells", "T\(", "Mast cells", "Treg\(diff\)", "DC"]
    return df.isin(known_immune_types)

def strip_integration_markers(barcode: str, valid_libs: List[str] = None) -> str:
    # for example from "AACTGTCAAGTCGT-1_CZI-IA11512684-1-2-10" returns "AACTGTCAAGTCGT-1_CZI-IA11512684"
    # Note this function assumes that the library id encoded in the barcode is of the form xxx-yyy (for example "CZI-IA11512684"),
    # otherwise, the behavior is undefined
    parts = barcode.split("_")
    cell_barcode = parts[0]
    lib_id_plus_integration_markers = parts[1].split("-")
    lib_id = "-".join([lib_id_plus_integration_markers[0], lib_id_plus_integration_markers[1]])
    if valid_libs is not None and lib_id not in valid_libs:
        raise ValueError(f"lib_id {lib_id} is not a valid library. Are you sure your barcode {barcode} contains a lib_id of the form xxx_yyy?")
    return "_".join([cell_barcode, lib_id])
   
def get_all_libs(lib_type: str, donor_id: Optional[str] = None) -> set:
    if lib_type not in ["GEX", "BCR", "TCR"]:
        raise ValueError("Unsupported lib_type: {}. Must be one of: GEX, BCR, TCR".format(lib_type))
    all_libs = set()
    samples = read_immune_aging_sheet("Samples")
    column_name = "{} lib".format(lib_type)
    libs_all = samples[column_name] if donor_id is None else samples[samples["Donor ID"] == donor_id][column_name]
    for i in range(len(libs_all)):
        if libs_all.iloc[i] is np.nan:
            continue
        libs = libs_all.iloc[i].split(",")
        all_libs.update(libs)
    return all_libs

def extend_removed_features_df(adata, obsm_key, exclude_df):
    if obsm_key not in adata.obsm:
        adata.obsm[obsm_key] = exclude_df.copy()
    else:
        assert(np.all(adata.obsm[obsm_key].index == exclude_df.index))
        adata.obsm[obsm_key] = pd.merge(
            left=adata.obsm[obsm_key],
            right=exclude_df,
            left_index=True,
            right_index=True,
            how="left", # this should not matter since we are asserting that the left and right df's have equal index
            validate="one_to_one",
            suffixes=("_left_merged", "_right_merged")
        )

# Utility function to rename TLN folders and files to LLN on S3
# as a consequence of renaming these samples in our IA Sample spreadsheet
def handle_sample_tln_to_lln_renaming(donor_id: str, delete_lln: bool = False):
    path_prefix = "s3://immuneaging/processed_samples/"
    def get_tln_lln_folders():
        ls_cmd = 'aws s3 ls {}'.format(path_prefix)
        files = os.popen(ls_cmd).read()
        lln_folders, tln_folders = [], []
        for f in files.rstrip().split('\n'):
            folder_name = f.split('PRE ')[-1]
            if folder_name.startswith(donor_id + "-LLN"):
                lln_folders.append(folder_name)
            elif folder_name.startswith(donor_id + "-TLN"):
                tln_folders.append(folder_name)
        return lln_folders, tln_folders
    lln_folders, tln_folders = get_tln_lln_folders()
    if len(lln_folders) > 0:
        if delete_lln:
            for lf in lln_folders:
                ls_cmd = 'aws s3 rm --recursive {}'.format(path_prefix+lf)
                os.popen(ls_cmd).read()
        else:
            print("LLN folder(s) present: {}".format(", ".join(lln_folders)))
            return
    if len(tln_folders) > 0:
        for tf in tln_folders:
            dest = tf.replace("TLN", "LLN")
            ls_cmd = 'aws s3 mv --recursive {} {}'.format(path_prefix+tf, path_prefix+dest)
            os.popen(ls_cmd).read()
    # rename TLN to LLN in file names
    lln_folders, _ = get_tln_lln_folders()
    for lf in lln_folders:
        ls_cmd = 'aws s3 ls --recursive {}'.format(path_prefix+lf)
        files = os.popen(ls_cmd).read()
        for f in files.rstrip().split('\n'):
            dir_name = f.split(' ')[-1]
            file_name = dir_name.split('/')[-1]
            dir_name_no_file_name = dir_name.removesuffix(file_name)
            new_file_name = file_name.replace("TLN", "LLN")
            ia_path = "s3://immuneaging/"
            ls_cmd = 'aws s3 mv {} {}'.format(ia_path+dir_name, ia_path+dir_name_no_file_name+new_file_name)
            print(ls_cmd)
            os.popen(ls_cmd).read()

def read_csv_from_aws(
        data_dir: str,
        aws_dir: str,
        aws_filename: str,
        logger: Type[BaseLogger],
    ):
    sync_cmd = 'aws s3 sync --no-progress {} {} --exclude "*" --include {}'.format(aws_dir, data_dir, aws_filename)
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
    file_path = os.path.join(data_dir, aws_filename)
    if not os.path.isfile(file_path):
        msg = "Failed to download file {} from S3.".format(aws_filename)
        logger.add_to_log(msg, level="error")
        raise ValueError(msg)
    return pd.read_csv(file_path)

# detects which datapoints in x have extreme values (does not count missing data as outliers)
def detect_outliers(x, num_sds):
    # a value is considered as an outlier if it is more extreme that the mean plus (or minus) num_sds times the standard deviation
    assert np.sum(x<0) == 0
    m, s = np.nanmean(x), np.nanstd(x)
    lower_bound = np.maximum(0, m - s*num_sds)
    upper_bound = m + s*num_sds
    is_not_outlier = np.logical_and(x >= lower_bound, x <= upper_bound)
    # do not consider missing values as outliers
    return np.logical_or(x.isna(), is_not_outlier), lower_bound, upper_bound

# find and return, among the given labels, those that constitute at least the given fraction (frac) of all the labels
def find_abundant_labels(labels, frac):
    labels.unique()
    labels_include = []
    th = frac*len(labels)
    unique_labels = np.unique(labels)
    for l in unique_labels:
        if np.sum(labels == l) > th:
            labels_include.append(l)
    return labels_include

def annotate(
    adata,
    model_paths,
    model_urls,
    components_key,
    neighbors_key,
    n_neighbors,
    resolutions,
    model_name,
    dotplot_min_frac,
    logger,
    save_all_outputs: bool = False
):
    adata_new = adata.copy()
    dotplot_paths = []
    for r in range(len(resolutions)):
        resolution = resolutions[r]
        logger.add_to_log("Running Leiden clustering using resolution={0}...".format(resolution))
        sc.pp.neighbors(adata_new, n_neighbors = n_neighbors, use_rep = components_key, key_added = neighbors_key)
        sc.tl.leiden(adata_new, resolution = resolution, key_added = 'leiden', neighbors_key = neighbors_key)
        # save the leiden clusters and majority voting results in the original anndata
        leiden_key_added = f"{model_name}.leiden_resolution_{str(resolution)}"
        adata.obs[leiden_key_added] = adata_new.obs["leiden"]
        for m in range(len(model_paths)):
            celltypist_model_name = model_urls[m].split("/")[-1].split(".")[0]
            model_path = model_paths[m]
            logger.add_to_log("Running CellTypist annotation using model {0}...".format(model_path))
            predictions = celltypist.annotate(adata_new, model = model_path, majority_voting = True, over_clustering = 'leiden')
            logger.add_to_log("Saving CellTypist outputs for model {0}...".format(model_path))
            adata.obs["celltypist_majority_voting.{0}.{1}.leiden_resolution_{2}".format(celltypist_model_name, model_name, str(resolution))] = predictions.predicted_labels["majority_voting"]
            if r == 0 and save_all_outputs:
                # save the rest of the outputs; these outputs do not change with different leiden resolution so should save only once
                adata.obs["celltypist_model_url.{0}".format(celltypist_model_name)] = model_urls[m]
                adata.obs["celltypist_predicted_labels.{0}".format(celltypist_model_name)] = predictions.predicted_labels["predicted_labels"]
                adata.obsm["celltypist_probability_matrix.{0}".format(celltypist_model_name)] = predictions.probability_matrix
            # generate and save a dotplot only for the abundant cell types
            logger.add_to_log("Generating a dotplot based on the CellTypist outputs...")
            abundant_cell_types = find_abundant_labels(labels = predictions.predicted_labels["predicted_labels"], frac = dotplot_min_frac)
            dotplot_predictions = celltypist.annotate(
                adata_new[predictions.predicted_labels["predicted_labels"].isin(abundant_cell_types),:].copy(),
                model = model_path,
                majority_voting = True,
                over_clustering = 'leiden'
            )
            celltypist.dotplot(
                dotplot_predictions,
                use_as_reference = 'leiden',
                use_as_prediction = 'predicted_labels',
                show = False,
                return_fig = False
            )
            dotplot_filename = "celltypist_dotplot.{}.{}.leiden_resolution_{}.min_frac_{}".format(
                celltypist_model_name,
                model_name,
                str(resolution),
                str(dotplot_min_frac)
            )
            sc.pl._utils.savefig_or_show(dotplot_filename, show = False, save = True)
            # note that celltypist (which uses scanpy for plotting) will only output the figures into a "figures" directory under the current working directory.
            dotplot_paths.append(os.path.join(os.getcwd(), "figures", dotplot_filename + ".pdf"))
    return dotplot_paths

def get_donor_id_for_lib(library_type, library_id, samples=None):
    if library_type not in ["GEX", "BCR", "TCR"]:
        raise ValueError("Unsupported lib_type: {}. Must be one of: GEX, BCR, TCR".format(library_type))
    if samples is None:
        samples = read_immune_aging_sheet("Samples")
    column_name = "{} lib".format(library_type)
    donor_id_col_name = "Donor ID"
    libs_all = samples[[column_name, donor_id_col_name]]
    for i in range(len(libs_all)):
        if libs_all.loc[i, column_name] is np.nan:
            continue
        libs = libs_all.loc[i, column_name].split(",")
        if library_id in libs:
            # found it
            return libs_all.loc[i, donor_id_col_name]
    return ""

def read_library(library_type, library_id, s3_access_file, working_dir, stage, logger, remove_adata=True, samples=None, donor_id=None):
    if donor_id is None:
        donor_id = get_donor_id_for_lib(library_type, library_id, samples)        
    # really hacky way to work around the fact that we don't have seq runs
    # at this layer. Almost all donors have 001 as the seq run, but a few
    # have 002 or 003 so try for those too
    for seq_run in ["001", "002", "003"]:
        if stage == "processed":
            file_name_partial = "{}_{}_{}_{}".format(donor_id, seq_run, library_type, library_id)
            s3_processed_lib_path = "s3://immuneaging/processed_libraries/{}".format(file_name_partial)
            version = get_latest_lib_version(s3_access_file, s3_processed_lib_path)
            file_name = "{}.processed.{}.h5ad".format(file_name_partial, version)
            sync_cmd = 'aws s3 sync --no-progress {}/{}/ {} --exclude "*" --include {}'.format(
                s3_processed_lib_path, version, working_dir, file_name
            )
        elif stage == "aligned":
            if library_type != "GEX":
                raise NotImplementedError()
            file_name_partial = "{}_{}.{}".format(donor_id, seq_run, library_id)
            s3_aligned_lib_path = "s3://immuneaging/aligned_libraries"
            folder_name = "{}_{}_{}_{}".format(donor_id, seq_run, library_type, library_id)
            version = get_latest_lib_version(s3_access_file, s3_aligned_lib_path, lib_folder_name=folder_name)
            file_name = "{}.{}.h5ad".format(file_name_partial, version)
            sync_cmd = 'aws s3 sync --no-progress {}/{}/{} {} --exclude "*" --include {}'.format(
                s3_aligned_lib_path, version, folder_name, working_dir, file_name
            )

        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        adata_file = os.path.join(working_dir, file_name)
        if os.path.isfile(adata_file):
            print(f"file {adata_file} already downloaded, skipping download")
            break

        logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
        if not os.path.isfile(adata_file):
            logger.add_to_log("Failed to download file {} from S3 using seq_run: {}".format(file_name, seq_run))
        else:
            break # don't need to try the other seq runs
    
    if not os.path.isfile(adata_file):
        logger.add_to_log("Failed to download file {} from S3.".format(file_name, seq_run))
        return None
        
    adata = anndata.read_h5ad(adata_file)
    if remove_adata:
        os.remove(adata_file)
    return adata
