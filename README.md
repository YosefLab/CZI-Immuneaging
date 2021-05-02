# The Immune Aging Data Hub

**(Work in progress...)**

This repository serves as the data hub of the Immune Aging project.
Here you can find all the resources and information needed for accessing and visualizing the currently available processed data of the project, as well as a complete documentation of the data collection process, including data upload and data processing.

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
    2. [Processing a single sample](#processing_sample)
    3. [Data Harmonization](#processing_harmonization)
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
If you will be uploading data to the bucket or will be using terminal commands for downloading data (i.e. rather using the AWS console for downloading; terminal download, which will be described later, allows to download large amount of files more conveniently), you will also need a user-specific credentials file that will allow you to access the S3 bucket via terminal.
<br /><br />
The Yosef group will manage the generation of both IAM users and credentials file for collaborators on the Immune Aging project. 
In order to set up an IAM user and receive credentials file please email Elior (erahmani@berkeley.edu) and Galen (gx2113@columbia.edu) and cc your Immune Aging PI, who will need to approve your request. Note that IAM accounts and credentials files will be issued with read access only for all users, except for designated data uploaders who will also get write access. In any case, **DO NOT SHARE THE CREDENTIALS FILE WITH ANYONE**.

1. Install AWS CLI (Optional for downloading data):

* Uploading data is done via the AWS Command Line Interface (CLI), which can be downloaded from <a href="https://aws.amazon.com/cli/">here</a>. While you will not need to explicitly use the AWS CLI for data upload, the upload script (to be described later) requires it as a dependency.
* Downloading multiple data files becomes more streamlined by using the AWS CLI rather than the <a href="https://911998420209.signin.aws.amazon.com/console">AWS console</a>, which is a web-based graphical user interface. For more details see [Data Download](#download).

<!--
**Granting temporary/one-off access:**
To share data with temporary members (like research assistants temporarily helping out) or for one-off access, email Galen (gx2113@columbia.edu) with the S3 folder you want to share (or upload data to) and he will generate a script that will automatically upload/download data when run. This is super easy for him to do (says Galen writing this), so don't hesitate to reach out. DO NOT SHARE YOUR USER SPECIFIC CREDENTIALS!
-->

### <a name="preliminaries_spreadsheet"></a> List of donors and samples

The project's donors and samples, as well as the raw metadata, can be found in the <a href="https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2" target="_blank">AI samples Google Spreadsheet</a>.
**DO NOT SHARE THIS LINK WITH ANYONE OUTSIDE THE IMMUNE AGING PROJECT**.

Note that:

* This spreadsheet does not reflect at any given moment the existing data on the S3 bucket or which samples have already been processed. Specifically, it may include information about donors and samples that either were not sequenced yet, are being currently sequenced, being uploaded, or pending data processing.
* The metadata of the donors and samples in the spreadsheet is also available in the processed data files (see below).


## <a name="download"></a> Data Download

### <a name="download_structure"></a> Directory structure on S3

Whether you will be downloading data from the S3 bucket [directly from the AWS consol](#download_console) or [via terminal](#download_terminal), you will need to know the bucket's directory structure.

The root directory if the Immune Aging S3 bucket includes the following sub-directories:

* raw_columbia/
* raw_sanger/
* processed/
* integrated/
* vision/

The first two directories are designated for raw data uploads (fastq.gz files).
The third and fourth directories, which are discussed in detail below, are processed versions of the raw data. The `vision` directory will be discussed later under Data Visualization.

The `processed` directory includes a `h5ad` file for each sample. This file can be used for downstream analysis in single-cell analysis tools such as scvi-tools [ref] and Seurat[ref]. 
Throughout the life cycle of the Immune Aging project we expect to improve our data processing pipelines. We therefore version the processed data by structuring the directories accordingly and by including a `.log` file with each processed sample that describes in details the processing pipeline used, such as the aligner used, software versions etc. The log file is important to guarantee reproducibility at the long-term, as we expect our pipeline to update from time to time, as new best practices in the field emerge.

TODO: add the `configs` dirs in processed and harmonized...
TODO: describe the cell ranger outputs and other outputs that are specific for each sample..

The structure of the `processed` directory is designed as follows:

* processed/
    * v1/
        * sample_name_1.processed.v1.h5ad
        * sample_name_1.processed.v1.h5ad.log
        * sample_name_2.processed.v1.h5ad
        * sample_name_2.processed.v1.h5ad.log
        * ...
    * v2/
        * sample_name_1.processed.v2.h5ad
        * sample_name_1.processed.v2.h5ad.log
        * sample_name_2.processed.v2.h5ad
        * sample_name_2.processed.v2.h5ad.log
        * ...
    * ...

For example, given two files `sample_name_1.processed.v1.h5ad, sample_name_1.processed.v2.h5ad`, we can tell by the file names that they were generated using a different data processing pipeline.
The exact differences and the complete information about each `h5ad` file can be found in the file's log. In addition, we further maintain a lookup table that maps each version to a concise summary of the pipeline used to generate the data file. You can find this table <a href="...">here</a>[todo].

The `harmonized` directory includes a `h5ad` with harmonized data containing multiple samples. Since the samples in the Immune Aging project are not generated and processed at the same time, for every new incoming sample or a batch of samples we will pool together and harmonize all existing data at the specific point in time.

The structure of the `harmonized` directory is designed as follows:
* harmonized/
    * v1/
        * data_name1.harmonized.v1.time_stamp1.h5ad
        * data_name1.harmonized.v1.time_stamp1.h5ad.log
        * data_name1.harmonized.v1.time_stamp2.h5ad
        * data_name1.harmonized.v1.time_stamp2.h5ad.log
        * data_name2.harmonized.v1.time_stamp1.h5ad
        * data_name2.harmonized.v1.time_stamp1.h5ad.log
        * data_name2.harmonized.v1.time_stamp2.h5ad
        * data_name2.harmonized.v1.time_stamp2.h5ad.log
        * ...
    * v2/
        * ...
    * ...

For example, given two files `data_name.harmonized.v1.03-01-2021.h5ad, data_name.harmonized.v1.03-15-2021.h5ad`, we can tell by the file names that both used version `v1` of the harmonization pipeline, however, given the two different time stamps, we know that each file includes a different subset of samples.
The exact information will appear in the matching log files of the data files, however, we also maintain a lookup table that maps each version and time stamp to a concise summary of the processing pipeline used and the samples that are included in the data file. You can find this table <a href="...">here</a>[todo].


### <a name="download_console"></a> Downloading data via the AWS console

Data that were already uploaded to the S3 bucket can be downloaded by logging in to AWS through <a href="https://911998420209.signin.aws.amazon.com/console.">this link</a> and navigating through the project's directory structure.


### <a name="download_terminal"></a> Downloading data via terminal

To read from the S3 bucket:
1. Run your user-specific credentials file to set aws keys
2. Sync your folders via `aws s3 sync <source> <target> [--options]` (where source is the aws folder). You can also use `aws s3 ls <target> [--options]` to list contents of a directory. Checkout more commands <a href= "https://docs.aws.amazon.com/cli/latest/userguide/cli-services-s3-commands.html">here</a>.



## <a name="upload"></a> Data upload

### <a name="upload_naming"></a> File naming conventions

New samples can only be uploaded as pre-aligned, gzip-compressed `fastq` files (i.e. `.fastq.gz` files). Each file must follow the following naming convention:
```
<donor_id>_<library_type>_<library_id>_<10x_reaction>_<seq_run>_<lane>_<read>.fastq.gz
```
<!--
where
* `donor_id` - the donor ID, should be consistent with the <a href="https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2">samples spreadsheet</a>; do not include a dot
* `organ` - the organ, which should be consistent with the Dictionary tab in the <a href="https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2">samples spreadsheet</a>; for example, using `SPL` for spleen and `LIV` for liver.
* `cell_type` the cell type, which should be consistent with the Dictionary tab in the <a href="https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2">samples spreadsheet</a>; for example, using `T` for T cells and `B` for B cells.
* `sequencer_output` - the standard output given by the raw data from the sequencer; for example, the output for a specific lane.
-->
The entries `<donor_id>` and `<library_id>` must be consistent with the donors sheet and the samples sheet of the <a href="https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/">IA Google spreadsheet</a>, and that library type can take a value out of the following five possible values: GEX (for gene expression), ADT (for CITE-seq), TCR (T-Cell Receptor), BCR (T-Cell Receptor), HTO (hashtag; in case a hashtag library was processed separately).

A few examples for valid file names include:
```
591C_ADT_CZI-IA10034945_0001_01_L001_R1.fastq.gz
583B_TCR_CZI-IA9924344_0001_01_L001_R1.fastq.gz
583B_GEX_CZI-IA9924327_0001_01_L001_R1.fastq.gz
```

<!--- 
Tcells-rest vs Tcells-stim, cell isolation (FAC-sorting, magnetic beads, etc)  information can go in obs field of h5ad
Use established naming convention for donors and tissues
 --->

### <a name="upload_upload"></a> Upload to S3

Once the Google spreadsheet is updated and the fastq files are ready to be uploaded, you can move forward to the actual data upload, which can be done using a designated script.

First, download the <a href="https://raw.githubusercontent.com/YosefLab/Immune-Aging-Data-Hub/main/scripts/upload.py?token=ABMIR5MFHUTH6CWXMNS4YGDASVMJW"> upload.py</a> script (open the link, right click on the screen, and save as "upload.py"). Second, make sure you have a working installation of python with the `pandas` library installed. If you do not have python with `pandas` already installed, you can simply download and install <a href="https://www.anaconda.com/products/individual">Anaconda</a>, a distribution of python that includes many popular libraries for data science, including `pandas`.

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
* Deleting files that have already been uploaded can only be done by an admin. If you believe that some raw files were uploaded by mistake and should be deleted please email Elior (erahmani@berkeley.edu) and Galen (gx2113@columbia.edu).


## <a name="visualization"></a> Data visualization

### <a name="visualization_live"></a> Live VISION Session

We provide a live VISION session of the latest version of the harmonized data.
VISION allows to... [ref]...
Every time a new sample or a batch of samples has been uploaded, the data will be processed and harmonized with all the existing samples in the project, followed by generating a new VISION session.

The most updated VISION sessions can be found here, separated by tissue:
* <a href="...">Lung</a> 
* <a href="...">Liver</a> 

### <a name="visualization_local"></a> Local VISION Session

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


## <a name="processing"></a> Data Processing

This section documents the pipeline for data processing that will be performed by the Yosef group members.

### <a name="processing_prerequisites"></a> Prerequisites

Download anaconda with python >= 3.7...make sure your anaconda is in the path etc.

...explain about the .yml files and versioning and why we use it...

For example,.. 
```
conda env create -f iadh.v1.march_2021.yml
```
This will create a new environment under `$CONDA_PREFIX"/envs/iadh.v1.march_2021"`

install cell ranger..

download reference genome, link to 10x

adding keys to the config files - (1) if adding keys that should not affect on versioning then update VARIABLE_CONFIG_KEYS in the processing files... (2) otherwise it will automatically generate a new version etc.. but need to make sure the script is backwards compatible...

### <a name="processing_sample"></a> Processing a single sample

The script `run_processing_pipeline.py` will...

For example, run:
```python
python run_processing_pipeline.py test_config_file.txt
```

This script initiates processing_pipeline.py which.. runs a job..

### <a name="processing_harmonization"></a> Data Harmonization

The script `harmonization_pipeline.py`... would make sense to run on the node with the GPU (s130)


## <a name="maintenance"></a> Data Hub maintenance


### <a name="maintenance_iams"></a> Setting up an AWS IAM

Explanations...

### <a name="maintenance_credentials"></a> Generating an AWS credentials file

Explanations...


---

## Authors & Maintainers of the Immune Aging Data Hub

Elior Rahmani <erahmani@berkeley.edu>

Galen Xing <gx2113@columbia.edu>

Allon Wagner <allonwag@berkeley.edu>

Nir Yosef <niryosef@berkeley.edu>
