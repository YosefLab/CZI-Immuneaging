import io
import json
from multiprocessing.sharedctypes import Value
import os
import time
import glob
import warnings
from jax import value_and_grad
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
from data_processing.scripts.logger import SimpleLogger
import gc

from logger import BaseLogger

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
    model_params_keys = ["use_layer_norm", "use_batch_norm"]
    model_params = dict()
    for i in model_params_keys:
        if i in configs:
            model_params[i] = configs[i]
    logger.add_to_log("Training {} model...".format(model_name))
    model = scvi.model.SCVI(adata, **model_params) if model_name=="scvi" else scvi.model.TOTALVI(adata, \
        empirical_protein_background_prior = empirical_protein_background_prior, **model_params)
    max_epochs_config_key = "scvi_max_epochs" if model_name=="scvi" else "totalvi_max_epochs"
    train_params_keys = ["lr","early_stopping","train_size","early_stopping_patience","batch_size","limit_train_batches"]
    train_params = dict()
    train_params["max_epochs"] = configs[max_epochs_config_key]
    for i in train_params_keys:
        if i in configs:
            train_params[i] = configs[i]
    model.train(**train_params)
    logger.add_to_log("Saving {} latent representation...".format(model_name))
    latent = model.get_latent_representation()
    if latent_key is None:
        latent_key = "X_scVI" if model_name=="scvi" else "X_totalVI"
    adata.obsm[latent_key] = latent
    model_file = "{}.{}.{}_model.zip".format(prefix, version, model_name)
    logger.add_to_log("Saving the model into {}...".format(model_file))
    model_file_path = os.path.join(data_dir, model_file)
    model_dir_path = os.path.join(data_dir,"{}.{}_model/".format(prefix, model_name))
    if os.path.isdir(model_dir_path):
        os.system("rm -r " + model_dir_path)
    model.save(model_dir_path)
    # save the data used for fitting the model; this is useful for applying reference-based integration on query data later on (based on the current model and data).
    logger.add_to_log("Saving the data used for fitting the model...")
    os.path.join(data_dir, model_file)
    data_file = "{}.{}.{}_model.data.h5ad".format(prefix, version, model_name)
    adata_copy = adata.copy()
    write_anndata_with_object_cols(adata_copy, model_dir_path, data_file)
    # zip the dir with all the model outputs
    zipf = zipfile.ZipFile(model_file_path, 'w', zipfile.ZIP_DEFLATED)
    zipdir(model_dir_path, zipf)
    zipf.close()
    return model, model_file

def aws_sync(source: str, target: str, include: str, logger: Type[BaseLogger]) -> None:
    sync_cmd = 'aws s3 sync --no-progress {} {} --exclude "*" --include {}'.format(source, target, include)
    logger.add_to_log("syncing {}...".format(include))
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    aws_response = os.popen(sync_cmd).read()
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
    assert np.sum(df.shape[1] == np.array(panel_sizes)) == 1 # make sure that twe can identify the protein panel used
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

def add_vdj_lib_ids_to_integrated_data(tissues_dir: str):
    bcr_to_gex, tcr_to_gex = get_vdj_lib_to_gex_lib_mapping()

    # make sure there is a 1:1 mapping betwen gex and bcr, and same for tcr
    gex = list(bcr_to_gex.values())
    if len(gex) != len(np.unique(gex)):
        raise ValueError("There is not a 1:1 mapping between gex libs and bcr libs")
    gex = list(tcr_to_gex.values())
    if len(gex) != len(np.unique(gex)):
        raise ValueError("There is not a 1:1 mapping between gex libs and tcr libs")

    gex_to_bcr = {v:k for k,v in bcr_to_gex.items()}
    gex_to_tcr = {v:k for k,v in tcr_to_gex.items()}

    files = os.listdir(tissues_dir)
    files = [f for f in files if f.endswith(".h5ad")]
    for file in files:
        file_path = os.path.join(tissues_dir, file)
        adata = anndata.read_h5ad(file_path)
        adata.obs["bcr_library_id"] = np.array([gex_to_bcr[l] if l in gex_to_bcr else "NA" for l in adata.obs["library_id"].values])
        adata.obs["tcr_library_id"] = np.array([gex_to_tcr[l] if l in gex_to_tcr else "NA" for l in adata.obs["library_id"].values])
        adata.write(file_path, compression="lzf")
        del adata
        gc.collect()

def get_vdj_lib_to_gex_lib_mapping():
    # Returns a mapping of all vdj libraries to their corresponding gex libraries
    # in the form of two dictionaries, one for bcr libs and one for tcr libs
    samples = read_immune_aging_sheet("Samples")

    def add_lib(lib_type: str, all_libs: dict) -> None:
        if lib_type not in ["BCR", "TCR"]:
            raise ValueError("Unsupported lib_type: {}. Must be one of: BCR, TCR".format(lib_type))
        column_name = "{} lib".format(lib_type)
        libs_all = samples[column_name]
        gex_libs_all = samples["GEX lib"]
        for i in range(len(libs_all)):
            if libs_all.iloc[i] is np.nan:
                continue
            libs = libs_all.iloc[i].split(",")
            gex_libs = gex_libs_all.iloc[i].split(",")
            for j in range(len(libs)):
                lib = libs[j]
                corresponding_gex_lib = gex_libs[j]
                all_libs[lib] = corresponding_gex_lib

    all_bcr_libs = {}
    add_lib("BCR", all_bcr_libs)

    all_tcr_libs = {}
    add_lib("TCR", all_tcr_libs)

    return all_bcr_libs, all_tcr_libs

def report_vdj_vs_cell_label_metrics_all_libs(ir_lib_type: str, tissues_dir: str):
    if ir_lib_type not in ["BCR", "TCR"]:
        raise ValueError("Unsupported lib_type: {}. Must be one of: BCR, TCR".format(ir_lib_type))
    is_b = ir_lib_type == "BCR"
    ir_libs = list(get_vdj_lib_to_gex_lib_mapping()[0 if is_b else 1].keys())
    samples = read_immune_aging_sheet("Samples")
    first = True
    for lib in ir_libs:
        # get some other metadata associated with this lib
        col_name = "{} lib".format(ir_lib_type)
        samples["Pick?"] = samples[col_name].fillna('').apply(lambda x: "Yes" if lib in x.split(",") else "No")
        idx = samples["Pick?"] == "Yes"
        sites = samples[idx]["Site"]
        donors = samples[idx]["Donor ID"]
        if len(set(sites)) > 1:
            raise ValueError("More than one site was found for lib id {} of type {}".format(lib, ir_lib_type))
        if len(set(donors)) > 1:
            raise ValueError("More than one donor was found for lib id {} of type {}".format(lib, ir_lib_type))
        site, donor = set(sites).pop(), set(donors).pop()
        # get adatas
        adatas = []
        files = os.listdir(tissues_dir)
        files = [f for f in files if f.endswith(".h5ad")]
        for file in files:
            file_path = os.path.join(tissues_dir, file)
            adata = anndata.read_h5ad(file_path)
            lib_id_col = "{}cr_library_id".format("b" if is_b else "t")
            if lib_id_col not in adata.obs.columns.values:
                raise ValueError("Make sure all tissue adata objects have {}".format(lib_id_col))
            adata_s = adata[adata.obs[lib_id_col] == lib, :].copy()
            if len(adata_s) > 0:
                adatas.append(adata_s)
            del adata
            gc.collect()
        if len(adatas) == 0:
            print("No cells found for library {}, moving on to the next library".format(lib))
            continue
        elif len(adatas) == 1:
            adata = adatas[0]
        else:
            adata = adatas[0].concatenate(adatas[1:], join="outer")
        report_vdj_vs_cell_label_metrics(adata, lib, ir_lib_type, site, donor, get_csv=True, skip_header=not first)
        if first:
            first = False
        del adata
        del adatas
        gc.collect()

def report_vdj_vs_cell_label_metrics_all_libs_v2(ir_lib_type: str, tissue_adatas: List[AnnData]):
    if ir_lib_type not in ["BCR", "TCR"]:
        raise ValueError("Unsupported lib_type: {}. Must be one of: BCR, TCR".format(ir_lib_type))
    is_b = ir_lib_type == "BCR"
    ir_libs = list(get_vdj_lib_to_gex_lib_mapping()[0 if is_b else 1].keys())
    samples = read_immune_aging_sheet("Samples")
    first = True
    for lib in ir_libs:
        # get some other metadata associated with this lib
        col_name = "{} lib".format(ir_lib_type)
        samples["Pick?"] = samples[col_name].fillna('').apply(lambda x: "Yes" if lib in x.split(",") else "No")
        idx = samples["Pick?"] == "Yes"
        sites = samples[idx]["Site"]
        donors = samples[idx]["Donor ID"]
        if len(set(sites)) > 1:
            raise ValueError("More than one site was found for lib id {} of type {}".format(lib, ir_lib_type))
        if len(set(donors)) > 1:
            raise ValueError("More than one donor was found for lib id {} of type {}".format(lib, ir_lib_type))
        site, donor = set(sites).pop(), set(donors).pop()
        # get adatas
        adatas = []
        for tdata in tissue_adatas:
            lib_id_col = "{}cr_library_id".format("b" if is_b else "t")
            if lib_id_col not in tdata.obs.columns.values:
                raise ValueError("Make sure all tissue adata objects have {}".format(lib_id_col))
            adata_s = tdata[tdata.obs[lib_id_col] == lib, :].copy()
            if len(adata_s) > 0:
                adatas.append(adata_s)
        if len(adatas) == 0:
            print("No cells found for library {}".format(lib))
            continue
        elif len(adatas) == 1:
            adata = adatas[0]
        else:
            adata = adatas[0].concatenate(adatas[1:], join="outer")
        report_vdj_vs_cell_label_metrics(adata, lib, ir_lib_type, site, donor, get_csv=True, skip_header=not first, skip_debug_print=True)
        if first:
            first = False
        del adata
        del adatas
        gc.collect()

def report_vdj_vs_cell_label_metrics(
    adata: AnnData,
    ir_lib_id: str,
    ir_lib_type: str,
    site: str,
    donor: str,
    get_csv: bool = False,
    skip_header: bool = False,
    skip_debug_print: bool = False
):
    def debug_print(msg: str):
        if not skip_debug_print:
            print(msg)

    ct_key_high = "celltypist_majority_voting.Immune_All_High.totalvi.leiden_resolution_2.0"
    ct_key_low = "celltypist_majority_voting.Immune_All_Low.totalvi.leiden_resolution_2.0"

    # Define known cell label categories
    # high hierarchy aka low resolution
    known_b_cells_high = ["B cells", "B-cell lineage", "Plasma cells"]
    known_t_cells_high = ["Double-negative thymocytes", "Double-positive thymocytes", "ETP", "T cells"]
    known_non_b_non_t_cells_high = [
        "DC",
        "DC precursor",
        "Endothelial cells",
        "Epithelial cells",
        "Erythrocytes",
        "Erythroid",
        "Fibroblasts",
        "Granulocytes",
        "ILC",
        "Monocyte precursor",
        "Monocytes",
        "HSC/MPP",
        "Macrophages",
        "Mast cells",
        "Megakaryocyte precursor",
        "Megakaryocytes/platelets",
        "Early MK",
        "MNP",
        "Mono-mac",
        "Myelocytes",
        "pDC",
        "pDC precursor",
        "Promyelocyte"
    ]
    known_non_b_cells_high = known_t_cells_high + known_non_b_non_t_cells_high
    known_non_t_cells_high = known_b_cells_high + known_non_b_non_t_cells_high
    # low hierarchy aka high resolution
    known_b_cells_low = ["Cycling B cells"]
    known_t_cells_low = ["Cycling gamma-delta T cells", "Cycling T cells", "Early lymphoid/T lymphoid"]
    known_non_b_non_t_cells_low = [
        "Cycling DCs",
        "Cycling monocytes",
        "Cycling NK cells"
    ]
    known_non_b_cells_low = known_t_cells_low + known_non_b_non_t_cells_low
    known_non_t_cells_low = known_b_cells_low + known_non_b_non_t_cells_low

    def run_analysis(lib_id: str, receptor_type: str):
        if receptor_type not in ["BCR", "TCR"]:
            raise ValueError("Unknown receptor_type passed: {}".format(receptor_type))
        cell_type = "b" if receptor_type == "BCR" else "t"
        is_b = receptor_type == "BCR"
        debug_print("Running for {} of type {}...".format(lib_id, receptor_type))

        # Look at labels of cells that have receptors of the given type 
        debug_print("Looking at labels of cells that have {} cell receptors...".format(cell_type))
        ir_cells_by_receptor = adata[adata.obs["{}-has_ir".format(receptor_type)] == "True", :].copy()
        num_ir_cells_by_receptor_1 = len(ir_cells_by_receptor)
        debug_print("Examining high hierarchy cell type labels...")
        ir_cells_by_type = ir_cells_by_receptor.obs[ct_key_high]
        # known b/t cell types high
        indices = ir_cells_by_type.isin(known_b_cells_high if is_b else known_t_cells_high)
        num_known_ir_cells_high = indices.sum()
        ir_cells_by_type = ir_cells_by_type[~indices]
        # known non b/t cell types high
        indices = ir_cells_by_type.isin(known_non_b_cells_high if is_b else known_non_t_cells_high)
        num_known_non_ir_cells_high = indices.sum()
        ir_cells_by_type = ir_cells_by_type[~indices]
        if len(ir_cells_by_type) == 0:
            debug_print("No cells left to examine in low hierarchy")
            # continue anyway - all the steps below will just return 0 in this case
        else:
            debug_print("Examining low hierarchy cell type labels...")
        # low hierarchy
        ir_cells_by_type_low = ir_cells_by_receptor[ir_cells_by_type.index, :].obs[ct_key_low]
        # known b/t cell types low
        indices = ir_cells_by_type_low.isin(known_b_cells_low if is_b else known_t_cells_low)
        num_known_ir_cells_low = indices.sum()
        ir_cells_by_type_low = ir_cells_by_type_low[~indices]
        # known non b/t cell types low
        indices = ir_cells_by_type_low.isin(known_non_b_cells_low if is_b else known_non_t_cells_low)
        num_known_non_ir_cells_low = indices.sum()
        if num_ir_cells_by_receptor_1 == 0:
            pct_known_ir_cells_high_n_low, pct_known_non_ir_cells_high_n_low = -1, -1
        else:
            pct_known_ir_cells_high_n_low = ((num_known_ir_cells_high + num_known_ir_cells_low) * 100) / num_ir_cells_by_receptor_1
            pct_known_non_ir_cells_high_n_low = ((num_known_non_ir_cells_high + num_known_non_ir_cells_low) * 100) / num_ir_cells_by_receptor_1

        draw_separator_line()

        # Look at presence of b/t cell receptors for cells that have b/t cell labels (will give us a hint about V(D)J library seq. saturation)
        debug_print("Looking at presence of {}'s for cells that have {} cell labels ...".format(receptor_type, cell_type))
        # known b/t cell types high and low
        indices_high = adata.obs[ct_key_high].isin(known_b_cells_high if is_b else known_t_cells_high)
        indices_low = adata.obs[ct_key_low].isin(known_b_cells_low if is_b else known_t_cells_low)
        indices = indices_high | indices_low
        ir_cells_by_type = adata[indices, :].copy()
        num_ir_cells_by_type = len(ir_cells_by_type)
        # look at b/t cell receptors
        ir_cells_by_receptor = ir_cells_by_type[ir_cells_by_type.obs["{}-has_ir".format(receptor_type)] == "True", :].copy()
        num_ir_cells_by_receptor_2 = len(ir_cells_by_receptor)
        # num_ir_cells_by_type can be 0 e.g. in the case of SKIN
        pct_ir_cells_by_receptor_2 = -1 if num_ir_cells_by_type == 0 else (num_ir_cells_by_receptor_2 * 100) / num_ir_cells_by_type
        # look for any cell receptors of the "opposite" cell type
        key = "TCR-has_ir" if is_b else "BCR-has_ir"
        opposite_cells_by_receptor = ir_cells_by_type[ir_cells_by_type.obs[key] == "True", :].copy()
        num_opposite_cells_by_receptor = len(opposite_cells_by_receptor)
        pct_opposite_cells_by_receptor = -1 if num_ir_cells_by_type == 0 else (num_opposite_cells_by_receptor * 100) / num_ir_cells_by_type

        # Look at presence of b/t cell receptors for cells that DO NOT have b/t cell labels  (will give us a hint about false positive receptors)
        debug_print("Looking at presence of {}'s for cells that DO NOT have {} cell labels ...".format(receptor_type, cell_type))
        # known non b/t cell types high and low
        indices_high = adata.obs[ct_key_high].isin(known_b_cells_high if is_b else known_t_cells_high)
        indices_low = adata.obs[ct_key_low].isin(known_b_cells_low if is_b else known_t_cells_low)
        indices = indices_high | indices_low
        non_ir_cells_by_type = adata[~indices, :].copy()
        num_non_ir_cells_by_type = len(non_ir_cells_by_type)
        # look at b/t cell receptors
        ir_cells_by_receptor = non_ir_cells_by_type[non_ir_cells_by_type.obs["{}-has_ir".format(receptor_type)] == "True", :].copy()
        num_ir_cells_by_receptor_3 = len(ir_cells_by_receptor)
        pct_ir_cells_by_receptor_3 = -1 if num_non_ir_cells_by_type == 0 else (num_ir_cells_by_receptor_3 * 100) / num_non_ir_cells_by_type

        draw_separator_line()

        debug_print("Results:")
        if not get_csv:
            print("Of {} {} cells by receptor:".format(num_ir_cells_by_receptor_1, cell_type))
            print("\t{} cells are of {} type in high hierarchy".format(num_known_ir_cells_high, cell_type))
            print("\t{} cells are not of {} cell type in high hierarchy".format(num_known_non_ir_cells_high, cell_type))
            print("\t{} cells are of {} type in low hierarchy".format(num_known_ir_cells_low, cell_type))
            print("\t{} cells are not of {} cell type in low hierarchy".format(num_known_non_ir_cells_low, cell_type))

            print("Of {} {} cells by type (high or low):".format(num_ir_cells_by_type, cell_type))
            print("\t{} cells have {} cell receptors".format(num_ir_cells_by_receptor_2, cell_type))
            print("\t{} cells have {} cell receptors".format(num_opposite_cells_by_receptor, "t" if is_b else "b"))

            print("Of {} non {} cells by type (high or low):".format(num_non_ir_cells_by_type, cell_type))
            print("\t{} cells have {} cell receptors".format(num_ir_cells_by_receptor_3, cell_type))

            print("cell label keys used: {}, {}".format(ct_key_high, ct_key_low))
        elif get_csv:
            csv_rows = [
                {
                    "Lib Id": lib_id,
                    "Site": site,
                    "Donor ID": donor,

                    "# {} cells by receptor".format(cell_type): num_ir_cells_by_receptor_1,
                    "Out of {} cells by receptor: # {} cells (high)".format(cell_type, cell_type): num_known_ir_cells_high,
                    "Out of {} cells by receptor: # {} cells (low)".format(cell_type, cell_type): num_known_ir_cells_low,
                    "Out of {} cells by receptor: pct {} cells (high&low)".format(cell_type, cell_type): pct_known_ir_cells_high_n_low,

                    "Out of {} cells by receptor: # non {} cells (high)".format(cell_type, cell_type): num_known_non_ir_cells_high,
                    "Out of {} cells by receptor: # non {} cells (low)".format(cell_type, cell_type): num_known_non_ir_cells_low,
                    "Out of {} cells by receptor: pct non {} cells (high&low)".format(cell_type, cell_type): pct_known_non_ir_cells_high_n_low,

                    "# {} cells by type".format(cell_type): num_ir_cells_by_type,
                    "Out of {} cells by type: # cells w/ {} cr".format(cell_type, cell_type): num_ir_cells_by_receptor_2,
                    "Out of {} cells by type: pct cells w/ {} cr (aka seq. saturation)".format(cell_type, cell_type): pct_ir_cells_by_receptor_2,
                    "Out of {} cells by type: # cells w/ {} cr".format(cell_type, "t" if is_b else "b"): num_opposite_cells_by_receptor,
                    "Out of {} cells by type: pct cells w/ {} cr".format(cell_type, "t" if is_b else "b"): pct_opposite_cells_by_receptor,

                    "# non {} cells by type".format(cell_type): num_non_ir_cells_by_type,
                    "Out of non {} cells by type: # cells w/ {} cr".format(cell_type, cell_type): num_ir_cells_by_receptor_3,
                    "Out of non {} cells by type: pct cells w/ {} cr (aka false positive)".format(cell_type, cell_type): pct_ir_cells_by_receptor_3,
                }
            ]
            csv_file = io.StringIO()
            field_names = list(csv_rows[0].keys())
            writer = csv.DictWriter(csv_file, fieldnames=field_names)
            if not skip_header:
                writer.writeheader()
            writer.writerows(csv_rows)
            print(csv_file.getvalue())
            csv_file.close()

    run_analysis(ir_lib_id, ir_lib_type)
