import json
import os
import time
import logging
import glob
import warnings
import pandas as pd
import scvi
import zipfile
from anndata._core.anndata import AnnData

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
	return

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
            warnings.showwarning(
                msg.message, msg.category, msg.filename, msg.lineno, msg.line
            )

    data = pd.read_csv(output_fn)
    os.remove(output_fn)
    return data

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
        **kwargs
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
	latent_key
		key to be used for saving the latent representation in adata.obsm.
    Returns
    -------
    A tuple containing the updated adata, trained model, and the name of the zip file where the model is saved.
    """
    assert model_name in ["scvi", "totalvi"]
    logger.add_to_log("Setting up {}...".format(model_name))
    scvi.data.setup_anndata(adata, batch_key=batch_key, protein_expression_obsm_key=protein_expression_obsm_key)
    params = dict()
    if "use_layer_norm" in configs:
        params["use_layer_norm"] = configs["use_layer_norm"]
    if "use_batch_norm" in configs:
        params["use_batch_norm"] = configs["use_batch_norm"]
    logger.add_to_log("Training {} model...".format(model_name))
    model = scvi.model.SCVI(adata, **params) if model_name=="scvi" else scvi.model.TOTALVI(adata, **params)
    max_epochs_config_key = "scvi_max_epochs" if model_name=="scvi" else "totalvi_max_epochs"
    params = dict()
    if "lr" in configs:
        model.train(max_epochs=configs[max_epochs_config_key], lr=float(configs["learning_rate"])) 
    else:
        model.train(max_epochs=configs[max_epochs_config_key])
    logger.add_to_log("Saving {} latent representation...".format(model_name))
    latent = model.get_latent_representation()
    latent_key = kwargs.get("latent_key", None)
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
    return adata, model, model_file

