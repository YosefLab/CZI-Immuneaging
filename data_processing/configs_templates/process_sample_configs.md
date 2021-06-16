## Configurations for process_library.py

The template for the configuration file can be found <a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/tree/main/data_processing/configs_templates/process_library.configs_file.example.txt">here </a>.

The configuration file is formatted as json with the following fields:

* `"sandbox_mode"` - `"True"` or `"False` to indicate whether running in sandbox mode or not. Always set to `"True"` if you are not an admin. Note that execution will fail and you will be alerted by the script if setting `"sandbox_mode": "True"` without AWS permissions of an admin.
* `"data_owner"` - name or username of the data owner (the person executing the process_library.py)
* `"code_path"` - Absolute path to the data processing scripts (i.e. if cloning the repository then should end with `"data_processing/scripts"`)
* `"output_destination"` - Absolute path that will be used for saving outputs
* `"s3_access_file"` - absolute path to the aws credentials file (provided by the admin)
* `"donor"` - Donor ID, as indicated in the Google Spreadsheet
* `"seq_run"` - Seq run, as indicated in the Google Spreadsheet
* `"library_type"` - Currently can only be set to `"GEX"`; this option will automatically take into account any ADT and/or HTO libraries associated with the library provided under `library_id`
* `"sample_id"` - The Sample_ID, as indicated in the Google Spreadsheet
* `"library_ids"` - A comma-separated (no spaces) list of the library IDs (as indicated in the Google Spreadsheet) that include the sample specified by `"sample_id"`
* `"processed_library_configs_version"` - A comma-separated (no spaces) list of the versions of the processed library files for the libraries specified by `"library_ids"` - these version numbers are determined by the configs version that was used to process the libraries; the latest alignment version of each processed library can be found on the S3 bucket under `s3://immuneaging/processed_libraries/`
* `"normalize_total_target_sum"` - A normalization factor; to be used for setting the total number of expression in each cell (if CITE seq data is available for the sample then will be applied to RNA and proteins separately)
* `"n_highly_variable_genes"` - The number of highly variable genes to be used prior to applying dimensionality reduction using PCA and SCVI
* `"highly_variable_genes_flavor"` - The flavor for identifying highly variable genes using `scanpy.pp.highly_variable_genes`
* `"scvi_max_epochs"` - The maximum number of epochs to be used when applying SCVI
* `"solo_max_epochs"` - The maximum number of epochs to be used when applying solo for doublet detection
* `"neighborhood_graph_n_neighbors"` - The number of neighbors to use for computing the neighborhood graph (using `scanpy.pp.neighbors`)
* `"umap_min_dist"` - The `min_dist` argument for computing UMAP (using `scanpy.tl.umap`)
* `"umap_spread"` - The `spread` argument for computing UMAP (using `scanpy.tl.umap`)
* `"umap_n_components"` - The number of UMAP components to compute (using `scanpy.tl.umap`)
 * `"python_env_version"` - The environment name to be used when running python commands for process_sample.py
