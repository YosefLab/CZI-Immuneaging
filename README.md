# The Immune Aging Data Hub

This repository serves as the data hub of the Immune Aging project.
Here you can find all the resources and information needed for accessing and visualizing the currently available processed data of the project, as well as a complete documentation of the data collection process, including data upload and data processing.

This page did not answer your question? Please <a href="...">open an issue</a> and label it as a 'question'.

## Table of contents
1. [Setting up data access](#access)
    1. [Prerequisites](#access_prerequisites)
    2. [Reading and Writing to the S3](#read_write)
    3. [List of donors and samples](#access_spreadsheet)
    4. [Downloading data from the S3 website](#access_s3)
    5. [Directory structure on S3](#access_s3_structure)
2. [Data upload](#upload)
    1. [Data naming conventions](#upload_naming)
    2. [Including metadata](#upload_metadata)
    3. [Upload to S3](#upload_upload)
3. [Data Download](#download)
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

## <a name="access"></a> Setting up data access

### <a name="access_prerequisites"></a> Prerequisites

The project data is stored on an Amazon Web Services (AWS) Simple Storage Service (S3; "S3 bucket"). Whether you will be uploading or just downloading data from the S3 bucket, you will need to set up your access via the following steps:

1. Get an IAM user and credentials file:
You will be given an Identity and Access Management (IAM) user account and a user-specific credentials file (DO NOT SHARE THIS FILE) to access the S3. An IAM user is associated with a specific user in the Immune Aging S3 bucket, and the credentials file allows using the data uploading scripts and downloading data via terminal (to be described later).
<br /><br />
The Yosef group will manage the generation of both IAM users and credentials file for collaborators on the Immune Aging project. 
In order to set up an IAM user and receive credentials file please email Galen (gx2113@columbia.edu) and cc your Immune Aging PI, who will need to approve your request.

1. Install AWS CLI (Optional but highly recommended):
Uploading and downloading data is best done via the AWS Command Line Interface (CLI) from <a href="https://aws.amazon.com/cli/">here</a>. This will also let you use our helper scripts for reading and writing data.

**Granting temporary/one-off access:**
To share data with temporary members (like research assistants temporarily helping out) or for one-off access, email Galen (gx2113@columbia.edu) with the S3 folder you want to share (or upload data to) and he will generate a script that will automatically upload/download data when run. This is super easy for him to do (says Galen writing this), so don't hesitate to reach out. DO NOT SHARE YOUR USER SPECIFIC CREDENTIALS!

### <a name="read_write"></a> Reading and Writing to the S3
To read data from the S3 bucket, use the `read_data.sh` script located in scripts [TODO]. See more in the [Data Download](#download) section.
To upload data to the S3 bucket, use the `upload_data.sh` script locatied in scripts [TODO]. See more in the [Data Upload](#upload) section.

### <a name="access_spreadsheet"></a> List of donors and samples

The list and raw metadata of the project's donors and samples can be found in <a href="https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2" target="_blank">this Google spreadsheet</a>.
In order to get access to this spreadsheet, please contact...[todo].

For data uploaders, this file serves as the samples organizer and should include a description of all the data that are being generated. 
Further instructions about naming conventions to be used in this spreadsheet for new samples are provided under the Data upload section.

Note that this spreadsheet does not reflect at any given moment the existing data on S3 or which samples have already been processed. Specifically, it may include information about donors and samples that either were not sequenced yet, are being currently sequenced, being uploaded, or pending data processing.

### <a name="access_s3"></a> Downloading data from the S3 website

Data that were already uploaded to the S3 bucket can be downloaded by logging in to AWS through <a href="https://911998420209.signin.aws.amazon.com/console.">this link</a> and navigating through the project's directory structure (see below for details).

In addition, we provide a simple script for automatically downloading processed data - either specific samples or harmonized datasets with multiple samples. See Section (<a name="data_download">Data download</a>.

### <a name="access_s3_structure"></a> Directory structure on S3

Whether you will be logging in to the 3S bucket for downloading data or using the download script (<a name="data_download">Data download</a> section), it will be helpful to understand the directory structure, which shows how we maintain difference versions of the data (e.g., data that went through different processing pipelines).

The data on S3 is structured as follows:

* raw_columbia/
* raw_sanger/
* processed/
* harmonized/
* vision/

The first two directories are designated for raw data uploads (.fastq files).
The third and fourth directories, which are discussed in detail below, are processed versions of the raw data. The `vision` directory will be discussed later under Data Visualization.

The `processed` directory includes a `h5ad` file for each sample. This file can be used for downstream analysis in single-cell analysis tools such as scvi-tools [ref] and Seurat[ref]. 
Throughout the life cycle of the Immune Aging project we expect to improve our data processing pipelines. We therefore version the processed data by structuring the directories accordingly and by including a `.log` file with each processed sample that describes in details the processing pipeline used, such as the aligner used, software versions etc. The log file is important to guarantee reproducibility at the long-term, as we expect our pipeline to update from time to time, as new best practices in the field emerge.

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

For example, given two files `data_name.harmonized.v1.Feb21.h5ad, data_name.harmonized.v1.April21.h5ad`, we can tell by the file names that both used version `v1` of the harmonization pipeline, however, given the two different time stamps, we know that each file includes a different subset of samples.
The exact information will appear in the matching log files of the data files, however, we also maintain a lookup table that maps each version and time stamp to a concise summary of the processing pipeline used and the samples that are included in the data file. You can find this table <a href="...">here</a>[todo].


## <a name="upload"></a> Data upload

### <a name="upload_naming"></a> Data naming conventions

New samples can only be uploaded as pre-aligned, gzip-compressed `fastq` files (i.e. `.fastq.gz` files). Each file should follow a specific naming convention, which includes donor ID, organ, and cell type (unsorted mononuclear cells [MNCs], T cells, B cells, or myeloid cells).
The upload process (described later) will now allow uploading files that do not follow the required naming convention.

Each file name must be structured as follows
```
<donor_id>_<organ>_<cell_type>.<sequencer_output>.fastq
```
where
* `donor_id` - the donor ID, which should be consistent with the <a href="https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2">samples spreadsheet</a>; do not include a dot
* `organ` - the organ, which should be consistent with the Dictionary tab in the <a href="https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2">samples spreadsheet</a>; for example, using `SPL` for spleen and `LIV` for liver.
* `cell_type` the cell type, which should be consistent with the Dictionary tab in the <a href="https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2">samples spreadsheet</a>; for example, using `T` for T cells and `B` for B cells.
* `sequencer_output` - the standard output given by the raw data from the sequencer; for example, the output for a specific lane.

A few examples for valid file names:
```
CUIMC-457_LNG_T.S1_L001_R1_001.fastq.gz
CUIMC-457_LNG_T.S1_L002_R1_001.fastq.gz
D456_COL_MNC.S1_L001_R1_001.fastq.gz
D456_COL_MNC.S1_L001_R2_001.fastq.gz
```

### <a name="upload_metadata"></a> Including metadata

TODO... 

<!--- 
Tcells-rest vs Tcells-stim, cell isolation (FAC-sorting, magnetic beads, etc)  information can go in obs field of h5ad
Use established naming convention for donors and tissues
 --->

### <a name="upload_upload"></a> Upload to S3

The very first step in uploading data is to update the <a href=https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2>samples spreadsheet</a> and name the new `fastq.gz` files following the naming conventions provided in the previous subsection.

Once the samples spreadsheet is updated and the files are ready, you can move on to the actual data upload.
We provide the script `data_upload.sh`, which uploads data files of a single sample by running the following command in terminal:
```
data_upload.sh <sample_prefix> <institute> <credentials_file>
```

Description of arguments:

* `<sample_prefix>` - the sample prefix, which is essentially the prefix `<donor_id>_<organ>_<cell_type>` that is described in the previous subsection. In case the data files of the samples are not in the current working directory you should include a path to data files.
* `<institute>` - can either be set to `columbia` or `sanger`.
* `<credentials_file>` - path to the user-specific credentials file (see Prerequisites subsection under Data Access).

For example, for uploading the files
```
CUIMC-457_LNG_T.S1_L001_R1_001.fastq.gz
CUIMC-457_LNG_T.S1_L002_R1_001.fastq.gz
```
on behalf of the Columbia group we run:
```
data_upload.sh /path/to/CUIMC-457_LNG_T columbia /path/to/iam_credentials.sh
```
where `/path/to/CUIMC-457_LNG_T` is the path fo the data files and the prefix of the sample that we want to upload and `/path/to/ai_credentials.sh` is the path to the credentials file of the user.

Finally, after uploading new data, please notify Elior Rahmani by email (erahmani@berkeley.edu). The Yosef team will then run the new sample through the processing pipeline.

Notes:

* `data_upload.sh` will not allow uploading data files that do not properly follow the naming convention; in such a case an error will be issued.
* Uploading more than one sample can be done by running the script multiple times, once per each sample.

## <a name="download"></a> Data Download

To read from the S3 bucket:
1. Run your user-specific credentials file to set aws keys
2. Sync your folders via `aws s3 sync <source> <target> [--options]` (where source is the aws folder). You can also use `aws s3 ls <target> [--options]` to list contents of a directory. Checkout more commands <a href= "https://docs.aws.amazon.com/cli/latest/userguide/cli-services-s3-commands.html">here</a>.

It is also possible to use the online AWS S3 interface but it is recommended to use the CLI. 


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
