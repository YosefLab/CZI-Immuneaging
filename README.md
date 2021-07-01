# The Immune Aging Data Hub

This repository serves as the data hub of the Immune Aging project.
Here you can find all the resources and information needed for accessing the currently available processed data of the project, as well as a complete documentation of the data collection process, including instructions for data upload (for designated data uploaders), data processing (for designated data owners) and other data management procedures (for designated system admins).

This page did not answer your question? Please <a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/issues">open an issue</a> and label it as a 'question'.

## Table of contents
1. [Preliminaries](#preliminaries)
    1. [Setting up data access](#preliminaries_access)
    2. [List of donors and samples](#preliminaries_spreadsheet)
2. [Data Download](#download)
    1. [Directory structure on S3](#download_structure)
    2. [Downloading data via the AWS console](#download_console)
    3. [Downloading data via terminal](#download_terminal)
3. [Data upload](#upload)
    1. [Prerequisites](#prerequisites)
    2. [File naming conventions](#upload_naming)
    3. [Upload to S3](#upload_upload)
4. [Data visualization](#visualization)
    1. [Live VISION Session](#visualization_live)
    2. [Local VISION Session](#visualization_local)
5. [Data Processing](#processing)
    1. [Prerequisites](#processing_prerequisites)
    2. [Processing libraries](#Processing_libraries)
    3. [Processing samples](#processing_samples)
    4. [Sandbox envorinment](#sandbox_envorinment)
    5. [Job queue](#job_queue)
6. [Data Hub maintenance](#maintenance)
    1. [Setting up an AWS IAM](#maintenance_iams)
    2. [Generating an AWS credentials file](#maintenance_credentials)

---

## <a name="preliminaries"></a> Preliminaries

### <a name="preliminaries_access"></a> Setting up data access

The project's data is stored on an Amazon Web Services (AWS) Simple Storage Service (S3; "S3 bucket"). Whether you will be uploading data to the Immune Aging S3 bucket or just using it for downloading data, you will need to set up your access as follows:

1. Get an IAM user and credentials file:
You will need an Identity and Access Management (IAM) user account, which will be associated with a specific user in the Immune Aging S3 bucket. Your user account will allow you to access the S3 bucket through the <a href="https://911998420209.signin.aws.amazon.com/console">AWS console</a>.
<br /><br />
If you will be uploading data to the bucket or will be using terminal commands for downloading data (i.e. rather than using the AWS console for downloading; terminal download, which will be described later, allows to download large amount of files more conveniently), you will also need a user-specific credentials file that will allow you to access the S3 bucket via terminal.
<br /><br />
The Yosef group will manage the generation of both IAM users and credentials file for collaborators on the Immune Aging project. 
In order to set up an IAM user and receive credentials file please email Elior Rahmani (erahmani@berkeley.edu) and cc your PI for approval of your request (the project's PIs are Menna Clatworthy, Donna Farber, Muzz Haniffa, Joanne Jones, Peter Sims, Sarah Teichmann, Roser Vento, and Nir Yosef). Note that IAM accounts and credentials files will be issued with read access only for all users, except for designated data uploaders and data processors who will also get write (i.e. upload) access. In any case, **DO NOT SHARE THE CREDENTIALS FILE WITH ANYONE**.

1. Install AWS CLI (Optional for downloading data):

* Uploading raw data (.fastq files) is done via the AWS Command Line Interface (CLI), which can be downloaded from <a href="https://aws.amazon.com/cli/">here</a>. While you will not need to explicitly use the AWS CLI for data upload, the upload script (to be described later) requires it as a dependency.
* Downloading multiple data files becomes more streamlined by using the AWS CLI rather than the <a href="https://911998420209.signin.aws.amazon.com/console">AWS console</a>, which is a web-based graphical user interface. For more details see [Data Download](#download).

**Note:** Non-US based users that wish to download data via AWS CLI or upload data using the upload tool described below may need to manually configure their location prior to any download/upload. Specifically, if your download/upload fails and you receive an error "Could not connect to the endpoint URL" then log in to the <a href="https://911998420209.signin.aws.amazon.com/console">AWS console</a> with your account and properly set the location in the top right corner of the console to match your physical location (e.g., users located in London should set `eu-west-2`).
If the problem persists, run in terminal the command `aws configure`. You will then be prompted to enter your `AWS Access Key ID` and `AWS Secret Access Key` - enter those keys as they appear in your credentials file. When prompted to enter your `Default region name`, enter the region name as you specified in the AWS console (e.g., London users should enter `eu-west-2`).

<!--
**Granting temporary/one-off access:**
To share data with temporary members (like research assistants temporarily helping out) or for one-off access, email Galen (gx2113@columbia.edu) with the S3 folder you want to share (or upload data to) and he will generate a script that will automatically upload/download data when run. This is super easy for him to do (says Galen writing this), so don't hesitate to reach out. DO NOT SHARE YOUR USER SPECIFIC CREDENTIALS!
-->

### <a name="preliminaries_spreadsheet"></a> List of donors and samples

The project's donors and samples, as well as the raw metadata, can be found in the <a href="https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2" target="_blank">AI samples Google Spreadsheet</a>.
You should never share this link with anyone outside the immune aging project. Also, **DO NOT MAKE EDITS IN THIS FILE**, unless you are a designated data uploader and you read all the instructions under [Data upload](#upload) - the data processing scripts rely on the data in this spreadsheet and some changes can compromise data integrity.

Note that:

* This spreadsheet does not reflect at any given moment the existing data on the S3 bucket or which samples have already been processed. Specifically, it may include information about donors and samples that either were not sequenced yet, are being currently sequenced, being uploaded, or pending data processing.
* The metadata of the donors and samples in the spreadsheet is also available in the processed data files (see below).


## <a name="download"></a> Data Download

### <a name="download_structure"></a> Directory structure on S3

Whether you will be downloading data from the S3 bucket [directly from the AWS consol](#download_console) or [via terminal](#download_terminal), you will need to know the bucket's directory structure.

The root directory of the Immune Aging S3 bucket includes the following sub-directories:

* aligned_libraries/
* job_queue/
* processed_libraries/
* processed_samples/
* raw_columbia/
* raw_sanger/
* test_folder/

Most of the users should mostly care about the `processed_samples`, which includes a processed data file for each sample in the project.
This directory, as well as the directories `aligned_libraries` and `processed_libraries` -  intermediate products of the data processing pipline - are discussed in detail below.
Briefly, the directories `raw_columbia` and `raw_sanger` are designated for raw data uploads (more details under [Data upload](#upload)) ,`job_queue` is a job queue to be used by data processors (more details under [Data Processing](#processing)), and `test_folder` is being used by the system admins for testing.

The `aligned_libraries` directory stores data of aligned libraries and supplementary files related to the alignment. Particularly, it includes a `h5ad` file for each sequenced library (a file format which can be used for downstream analysis in single-cell analysis tools such as <a href="https://scvi-tools.org/">scvi-tools</a> and <a href="https://satijalab.org/seurat/">Seurat</a>).
Throughout the life cycle of the Immune Aging project we expect to improve our data processing pipelines. We therefore version each intermediate step of processed data by structuring the directories accordingly and by including a `.log` file with each processed data file that describes in details the processing pipeline used. The log file is important to guarantee reproducibility at the long-term, as we expect our pipeline to update from time to time, as new best practices in the field emerge.

The structure of the `aligned_libraries` directory is as follows:

* aligned_libraries/
    * configs/
        * align_library.v1.txt
        * align_library.v2.txt
        * ...
    * v1/
        * library1
            * library1.cellranger.cloupe.cloupe
            * library1.cellranger.feature_ref.csv
            * library1.cellranger.libraries.csv
            * library1.cellranger.metrics_summary.csv
            * library1.cellranger.web_summary.html
            * library1.v1.h5ad
            * align_library.library1.v1.log
        * library2
            * library2.cellranger.cloupe.cloupe
            * library2.cellranger.feature_ref.csv
            * library2.cellranger.libraries.csv
            * library2.cellranger.metrics_summary.csv
            * library2.cellranger.web_summary.html
            * library2.v1.h5ad
            * align_library.library2.v1.log
        * ...
    * v2/
        * ...
    * ...

The `configs` directory includes versioned configuration files for the alingment pipeline. For each version of the aligned data, one designated directory (e.g., directory `v1` for version 1) includes the alingment outputs for each library. Particularly, it includes a .h5ad file, output files from cellranger, and a .log filw documenting the execution of the pipeline on the library. Note: libraries are not named arbitrarily as library1, library 2 etc. but rather take the following naming convention: `<donor_id>_<seq_run>_<library_type>_<library_id>`.

The `processed_libraries` directory stores data of processed libraries (i.e. beyond alingment), an intermediate product before the sample-level processing; its structure is as follows:

* processed_libraries/
    * library1/
        * v1/
            * library1.processed.v1.h5ad
            * process_library.library1.v1.log
            * process_library.configs.library1.v1.txt
        * v2/
            * library1.processed.v2.h5ad
            * process_library.library1.v2.log
            * process_library.configs.library1.v2.txt
        * ...
    * library2/
        * ...
    * ...

Here, files with the prefix `process_library.configs.` incude the configurations that were used in the execution of the library processing pipeline, and files with a `.log` suffix are documentation of the execution of the pipeline. As in the `aligned_libraries` directory, the actual library names follow the naming convention `<donor_id>_<seq_run>_<library_type>_<library_id>`.

The `processed_samples` directory stores data of processed samples (i.e. sample-level integration across different libraries); its structure is as follows:

* processed_samples/
    * sample1/
        * v1/
            * sample1.processed.v1.h5ad
            * sample1.processed.v1.scvi_model.zip
            * process_sample.sample1.v1.log
            * process_sample.configs.sample1.v1.txt
        * v2/
            * sample1.processed.v2.h5ad
            * sample1.processed.v2.scvi_model.zip
            * process_sample.sample1.v2.log
            * process_sample.configs.sample11.v2.txt
        * ...
    * sample2/
        * ...
    * ...

Here, files with the prefix `process_sample.configs.` incude the configurations that were used in the execution of the sample processing pipeline, and files with a `.log` suffix are documentation of the execution of the pipeline. 
Actual sample names follow the naming convention `<sample_id>_<data_type>`.

### <a name="download_console"></a> Downloading data via the AWS console

Data can be downloaded from the Immune Aging S3 bucket by logging in to AWS through <a href="https://911998420209.signin.aws.amazon.com/console.">this link</a> and navigating through the project's directory structure.

### <a name="download_terminal"></a> Downloading data via terminal

To read from the S3 bucket via terminal:
1. Run your user-specific credentials file to set the aws keys (note that the credential file is a shell file, which can be executed in a linux/mac environment).
2. Sync directories to your local machine via `aws s3 sync <source> <target> [--options]` (where source is the aws folder). You can also use `aws s3 ls <target> [--options]` to list contents of a directory. Check out more commands <a href= "https://docs.aws.amazon.com/cli/latest/userguide/cli-services-s3-commands.html">here</a>.



## <a name="upload"></a> Data upload

### <a name="upload_naming"></a> File naming conventions

New data can only be uploaded as pre-aligned, gzip-compressed `fastq` files (i.e. `.fastq.gz` files). Each file must follow the following naming convention:
```
<donor_id>_<seq_run>_<library_type>_<library_id>_<S_number>_<lane>_<read>_001.fastq.gz
```
<!--
where
* `donor_id` - the donor ID, should be consistent with the <a href="https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2">samples spreadsheet</a>; do not include a dot
* `organ` - the organ, which should be consistent with the Dictionary tab in the <a href="https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2">samples spreadsheet</a>; for example, using `SPL` for spleen and `LIV` for liver.
* `cell_type` the cell type, which should be consistent with the Dictionary tab in the <a href="https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2">samples spreadsheet</a>; for example, using `T` for T cells and `B` for B cells.
* `sequencer_output` - the standard output given by the raw data from the sequencer; for example, the output for a specific lane.
-->
Note that the suffix `<library_id>_<S_number>_<lane>_<read>_001` is the expected file name format of the raw sequencing output; for `<read>` only the values `R1` or `R2` are allowed.
The entries `<donor_id>`, `<library_type>`, and `<library_id>` must be consistent with the donors sheet and the samples sheet of the <a href="https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/">IA Google spreadsheet</a>, and that library type can take a value out of the five following options: GEX (for gene expression), ADT (for CITE-seq), TCR (T-Cell Receptor), BCR (T-Cell Receptor), HTO (hashtag; in case a hashtag library was processed separately). The field `<seq_run>` is a site-specific value that can be set to any 3-digit value in the range 001-999 (e.g., it can always be set to 001, or it can be increased by 1 for every new donor), however, in case of resequenced libraries that are designated to be processed together with an initial sequencing version of the same libraries, this number should be the same one used for the libraries from the initial sequencing.

A few examples for valid file names include:
```
582C_001_HTO_CZI-IA9924321_S1_L001_R1_001.fastq.gz
582C_002_ADT_CZI-IA9924369_S1_L001_R1_001.fastq.gz
591C_010_GEX_CZI-IA10034921_S2_L002_R2_001.fastq.gz
```

<!--- 
Tcells-rest vs Tcells-stim, cell isolation (FAC-sorting, magnetic beads, etc)  information can go in obs field of h5ad
Use established naming convention for donors and tissues
 --->

### <a name="upload_upload"></a> Upload to S3

Once the Google spreadsheet is updated and the fastq files are ready to be uploaded, you can move forward to the actual data upload, which can be done using a designated script.

First, download the <a href="https://raw.githubusercontent.com/YosefLab/Immune-Aging-Data-Hub/main/scripts/upload.py?token=ABMIR5MFHUTH6CWXMNS4YGDASVMJW"> upload.py</a> script (open the link, right click on the screen, and save as "upload.py"). Please make sure you always use the most updated version of the script. Second, make sure you have a working installation of python with the `pandas` library installed. If you do not have python with `pandas` already installed, you can simply download and install <a href="https://www.anaconda.com/products/individual">Anaconda</a>, a distribution of python that includes many popular libraries for data science, including `pandas`.

Simple usage of the upload.py script - type in terminal:
```
python upload.py --aws_keys <aws_keys> --destination <destination> --fastq <fastq>
```

Arguments:

* `<aws_keys>` - path to the user-specific credentials file (see details under [Setting up data access](#preliminaries_access)).
* `<destination>` - can either be set to `columbia` or `sanger`; `test` is also allowed - for testing purposes (will upload data to a directory named `test_folder` on the S3 bucket).
* `<fastq>` - path to folder containing fastq.gz files to be uploaded (**note:** all `fastq.gz` files in the provided folder will be uploaded).

For example:
```
python upload.py --aws_keys my_credentials.sh --destination sanger --fastq 390C_files/
```

At any time, you can see a description of the arguments accepted by the upload script by running:
```
python scripts/upload.py -h
```

The script runs several validations before performing the actual upload, including checking the metadata in the Google spreadsheet for consistency and verifying that all file names properly follow the naming convention. The script will issue detailed warnings and errors in case any of those validations fails (in which case no files will be uploaded).

Finally, after uploading new data, please notify Elior (erahmani@berkeley.edu). The Yosef team will then run the new data through the processing pipeline.

Notes:
* The way the upload.py script works is by syncing the data in the user-specified folder (provided via the `--fastq` argument) with the S3 bucket, rather than performing a naive upload. In case of a lost connection while the script is running, this mechanism allows the script to automatically resume the upload without re-uploading files that are already in the bucket.
* Deleting files that have already been uploaded can only be done by an admin. If you believe that some raw files were uploaded by mistake and should be deleted please email Elior Rahmani (erahmani@berkeley.edu).

## <a name="visualization"></a> Data visualization

### <a name="visualization_live"></a> Live VISION Session

Coming soon...
<!--
We provide a live VISION session of the latest version of the harmonized data.
VISION allows to... [ref]...
Every time a new sample or a batch of samples has been uploaded, the data will be processed and harmonized with all the existing samples in the project, followed by generating a new VISION session.

The most updated VISION sessions can be found here, separated by tissue:
* <a href="...">Lung</a> 
* <a href="...">Liver</a> 
-->

### <a name="visualization_local"></a> Local VISION Session

Coming soon...
<!--
VISION can also be launched locally given a vision object.
We keep the previous VISION objects (i.e. that were previously displayed in the live session) as well as the latest vision object on the S3 bucket.

In more detail the S3 bucket includes a directory named `vision`, which follows the following structure:
* vision/
    * ...

A VISION object can be downloaded by using the `sync` AWS CLI command as follows:
```
...
```

For example...
```
...
```

Then, in order to start a VISION session locally, we run:
```
...
```

-->

## <a name="processing"></a> Data Processing

This section documents the pipeline for data processing (post-alignment processing of libraries and sample-level integration).
Data processing can be performed by data owners in a sandbox environment (on a local machine), which allows to fine tune configurations. One such configurations are curated by a data owner, they can be uploaded to the S3 bucket into a queue, from which a system admin from the Yosef group will execute the processing (based on the configurations provided by the data owners). 

### <a name="processing_prerequisites"></a> Prerequisites

First, Download anaconda with python >= 3.7.
Second, get the .yml file of the latest version of the python environment from <a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/tree/main/envs">here</a> (i.e. immune_aging.py_env.*.yml where * is the latest version available), and use it to set up a new conda environment via terminal.

For example:
```
conda env create -f immune_aging.py_env.v1.yml
```
This will create a new environment under `$CONDA_PREFIX"/envs/immune_aging.py_env.v1.yml"`

### <a name="processing_libraries"></a> Processing libraries

The script `run_processing_pipeline.py` will...

For example, run:
```python
python run_processing_pipeline.py test_config_file.txt
```

This script initiates processing_pipeline.py which.. runs a job..

### <a name="processing_samples"></a> Processing samples

The script `harmonization_pipeline.py`... would make sense to run on the node with the GPU (s130)

### <a name="sandbox_envorinment"></a> Sandbox envorinment

### <a name="job_queue"></a> Job queue


## <a name="maintenance"></a> Data Hub maintenance


### <a name="maintenance_iams"></a> Setting up an AWS IAM

Explanations...

### <a name="maintenance_credentials"></a> Generating an AWS credentials file

Explanations...

### align libraries

install cellranger and download ref genome..

adding keys to the config files - (1) if adding keys that should not affect on versioning then update VARIABLE_CONFIG_KEYS in the processing files... (2) otherwise it will automatically generate a new version etc.. but need to make sure the script is backwards compatible...


---

## Authors & Maintainers of the Immune Aging Data Hub

Elior Rahmani <erahmani@berkeley.edu>

Galen Xing <gx2113@columbia.edu>

Allon Wagner <allonwag@berkeley.edu>

Nir Yosef <niryosef@berkeley.edu>
