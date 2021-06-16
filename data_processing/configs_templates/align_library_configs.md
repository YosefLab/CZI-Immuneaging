## Configurations for align_library.py

The template for the configuration file can be found <a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/tree/main/data_processing/configs_templates/align_library.configs_file.example.txt">here </a>.

The configuration file is formatted as json with the following fields:

* `"donor"` - Donor ID, as indicated in the Google Spreadsheet
* `"seq_run"` - Seq run, as indicated in the Google Spreadsheet
* `"output"`_destination - Absolute path that will be used for saving outputs
* `"aligner"` - The file name of the aligner executable file; currently only `"cellranger"` is allowed
* `"aligner_version"` - Version of the aligner used
* `"aligner_software_path"` - absolute path to the directory containing the executable of the aligner
* `"alignment_ref_genome_file"` - File name of the reference genome file to be used for alignment
* `"alignment_ref_genome_path"` - Absolute path to the directory containing the provided reference genome file
* `"berkeley_user"` - The Berkeley username of the person executing the align_library.py script
* `"s3_access_file"` - absolute path to the aws credentials file (provided by the admin)
* `"python_env_version"` - The environment name to be used when running python commands for align_library.py;
* `"r_env_version"` - The environment name to be used when running R commands for align_library.py; this environment is currently unused in align_library.py
* `"r_setup_version"` - Version of the setup file for additional R setups on top of those defined in `r_env_version`; this is currently unused in align_library.py
* `"code_path"` - Absolute path to the data processing scripts (i.e. if cloning the repository then should end with data_processing/scripts)
