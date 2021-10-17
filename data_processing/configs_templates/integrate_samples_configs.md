## Configurations for integrate_samples.py

The template for the configuration file can be found <a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/tree/main/data_processing/configs_templates/integrate_samples.configs_file.example.txt">here </a>.

The configuration file is formatted as json with the following fields:

* `"sandbox_mode"` - `"True"` or `"False` to indicate whether running in sandbox mode or not. Always set to `"True"` if you are not an admin. Note that execution will fail and you will be alerted by the script if setting `"sandbox_mode": "False"` without AWS permissions of an admin.
* `"data_owner"` - name or username of the data owner (the person executing integrate_samples.py)
* `"code_path"` - Absolute path to the data processing scripts (if cloning the repository then should end with `"data_processing/scripts"`)
* `"output_destination"` - Absolute path that will be used for saving outputs
* `"s3_access_file"` - absolute path to the aws credentials file (provided by the admin)
* `"output_prefix"` - prefix for the output files
* `"sample_ids"` - A comma-separated (no spaces) list of the sample IDs (`Sample_ID` field in the Google spreadsheet) to integrate
* `"processed_sample_configs_version"` - A comma-separated (no spaces) list of the versions of the processed sample files for the samples specified by `"sample_ids"` - these versions are determined by the configs version that was used to process the samples; the latest version of each processed sample can be found on the S3 bucket under `s3://immuneaging/processed_samples/`
* `"n_highly_variable_genes"` - The number of highly variable genes to be used prior to applying dimensionality reduction using PCA and SCVI
* `"highly_variable_genes_flavor"` - The flavor for identifying highly variable genes using `scanpy.pp.highly_variable_genes`
* `"scvi_max_epochs"` - The maximum number of epochs to be used when applying SCVI
* `"totalvi_max_epochs"` - The maximum number of epochs to be used when applying totalVI
* `"use_layer_norm"` - Whether to use layer norm in layers when fitting the SCVI and totalVI models; for default behavior set to `"none"`.
* `"use_batch_norm"` - whether to use batch norm when fitting the SCVI and totalVI models; for default behavior set to `"both"`.
* `"early_stopping"` - `"True"` or `"False"` to indicate whether to use early stopping in the model training of scvi and totalVI; for large data it is advised to use `"True"`.
* `"early_stopping_patience"` - Number of validation epochs with no improvement after which training of SCVI and totalVI will be stopped.
* `"batch_size"` - Minibatch size to use during training of SCVI and totalVI.
* `"limit_train_batches"` - Limit on the number of training batches during one epoch of training of SCVI and totalVI.
* `"neighborhood_graph_n_neighbors"` - The number of neighbors to use for computing the neighborhood graph (using `scanpy.pp.neighbors`)
* `"umap_min_dist"` - The `min_dist` argument for computing UMAP (using `scanpy.tl.umap`)
* `"umap_spread"` - The `spread` argument for computing UMAP (using `scanpy.tl.umap`)
* `"umap_n_components"` - The number of UMAP components to compute (using `scanpy.tl.umap`)
* `"python_env_version"` - The environment name to be used when running process_sample.py
* `"r_setup_version"` - Version of the setup file for additional R setups on top of those defined in `python_env_version`