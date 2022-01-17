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
* `"library_type"` - Can be set to `"GEX"`, `"BCR"` or `"TCR"`; `"GEX"` will automatically take into account any ADT and/or HTO libraries associated with the library provided under `library_id`
* `"library_id"` - The library ID, as indicated in the Google Spreadsheet
* `"corresponding_gex_lib"` - The library ID of the GEX lib that corresponds to this library (i.e. sequenced GEX data for the same cells). Relevant only for BCR and TCR libraries - for other lib types this is optional (would be the same as library_id).
* `"filter_cells_min_genes"` - Cells with number of detectable genes below this threshold will be removed
* `"filter_genes_min_cells"` - Genes that appear in less cells than this threshold will be removed
* `"filter_cells_max_pct_counts_mt"` - Cells with mitochondrial content above this threshold will be removed (value should be in the range [0,100])
* `"filter_cells_min_pct_counts_ribo"` - Cells with ribosomal content above this threshold will be removed (value should be in the range [0,100])
* `"genes_to_exclude"` - A comma-separated (no spaces) list of genes to exclude regardless of other quality procedures; set to `"None"` if no such genes
* `"exclude_mito_genes"` - Set `"True"` or `"False"` to indicate whether mitochondrial genes should be excluded regardless of other quality procedures
* `"hashsolo_priors"` - A comma-separated (no spaces) list of priors for hashsolo; the values are the expected fractions of multiplets, singlets, and doublets, respectively
* `"aligned_library_configs_version"` - The alignment version of the library to process - this version number is determined by the configs version that was used to align the library; the latest alignment version of each aligned library can be found on the S3 bucket under `s3://immuneaging/aligned_libraries`
* `"python_env_version"` - The environment name to be used when running process_library.py
* `"r_setup_version"` - Version of the setup file for additional R setups on top of those defined in `python_env_version`

