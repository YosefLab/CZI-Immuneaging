import json
import os
import time
import glob
import warnings
import pandas as pd
import scvi
import zipfile
from anndata._core.anndata import AnnData
from math import floor
import csv
from typing import Type
from logger import BaseLogger

AUTHORIZED_EXECUTERS = ["b750bd0287811e901c88dc328187e25f", "1c75133ab6a1fc3ed9233d3fe40b3d73"] # md5 checksums of the AWS_SECRET_ACCESS_KEY value of those that are authorized to upload outputs of processing scripts to the server; note that individuals with upload permission to aws can bypass that by changing the code - this is just designed to alert users that they should only use sandbox mode.

# these are formatted strings that we use to log QC stats and use them to parse those lines back
# if you edit these make sure to update all call sites that use them
QC_STRING_DOUBLETS = "Removed {} estimated doublets (percent removed: {:.2f}%); {} droplets remained."
QC_STRING_AMBIENT_RNA = "Removed {} cells (percent removed: {:.2f}%) with total decontaminated counts below filter_decontaminated_cells_min_genes={}"
QC_STRING_VDJ = "Removed {} vdj genes (percent removed: {:.2f}%); {} genes remained."
QC_STRING_RBC = "Removed {} red blood cells (percent removed: {:.2f}%); {} droplets remained."

def init_scvi_settings():
    # This does two things:
    # 1. Makes the logger look good in a log file
    # 2. Changes a bit how torch pins memory when copying to GPU, which allows you to more easily run models in parallel with an estimated 1-5% time hit
    scvi.settings.reset_logging_handler()
    scvi.settings.dl_pin_memory_gpu_training = False

init_scvi_settings()

def get_current_time():
	return time.strftime("%H:%M, %m-%d-%Y")

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

def run_model(
        adata: AnnData,
        configs: dict,
        batch_key: str,
        protein_expression_obsm_key: str,
        model_name: str,
        prefix: str,
        version: str,
        data_dir: str,
        logger,
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
    model_params_keys = ["use_layer_norm", "use_batch_norm"]
    model_params = dict()
    for i in model_params_keys:
        if i in configs:
            model_params[i] = configs[i]
    logger.add_to_log("Training {} model...".format(model_name))
    model = scvi.model.SCVI(adata, **model_params) if model_name=="scvi" else scvi.model.TOTALVI(adata, **model_params)
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

def filter_vdj_genes(rna: AnnData, aws_file_path: str, data_dir: str, logger: Type[BaseLogger]) -> None:
    file_path_components = aws_file_path.split("/")
    file_name = file_path_components[-1]
    aws_dir_path = "/".join(file_path_components[:-1])
    aws_sync(aws_dir_path, data_dir, file_name, logger)
    local_file_path = os.path.join(data_dir, file_name)
    with open(local_file_path) as csvfile:
        reader = csv.reader(csvfile)
        genes = [row[0] for row in reader]
    n_var_before = rna.n_vars
    rna = rna[:, ~rna.var.index.isin(genes)]
    percent_removed = 100*(n_var_before-rna.n_vars)/n_var_before
    level = "warning" if percent_removed > 50 else "info" # TODO adjust threshold if needed
    logger.add_to_log(QC_STRING_VDJ.format(n_var_before-rna.n_vars, percent_removed, rna.n_vars), level=level)
