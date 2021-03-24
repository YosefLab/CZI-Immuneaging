# The Immune Aging Data Hub

This repository serves as the data hub of the Immune Aging project.
Here you can find resources and all the information needed for accessing and visualizing the currently available processed data of the project, as well as a complete documentation for the data collection process, including data upload, data processing, and maintaining this data hub.

This page did not answer your question? Please <a href="...">open an issue</a> and label it as a 'question'.

## Table of contents
1. [Setting up data access](#access)
    1. [Prerequisites](#access_prerequisites)
    2. [S3 Structure](#s3_structure)
2. [Data Upload](#upload)
    1. [Updating the Google Spreadsheet](#spreadsheet)
    2. [Naming Convention](#naming_convention)
    3. [Upload to S3](#upload_s3)
3. [Data Download](#download)
4. [Data visualization](#visualization)
    1. [Live VISION sessions](#live_vision)
    2. [Local VISION sessions](#local_vision)
5. [Data processing](#processing)
    1. [Prerequisites](#processing_prerequisites)
6. [Data Hub Maintainers](#maintainers)
    1. [Generating AWS Credentials File](#generating_keys)

---

## <a name="access"></a> Setting up data access

### <a name="prerequisites"></a> Prerequisites

The data of the project are stored on an Amazon Web Services (AWS) Simple Storage Service (S3; "S3 bucket"). Whether you will be uploading or just downloading data from the S3 bucket, you will need an account and a user-specific credentials file.

[Galen]
For opening an account...
The use-specific credentials file is a simple text file that will be provided by the Yosef lab. In order to get a key, please contact: ...

In addition to setting up an AWS account and getting a credentials file, you will need to download the AWS  Command Line Interface (CLI) from <a href="https://aws.amazon.com/cli/">here</a>.

Finally, uploading or downloading data requires having the scripts on your system. [explanations how to download from the github]. Note that these scripts can only be executed on Linux/Unix/MacOSx.

### <a name="spreadsheet"></a> List of samples

The list of donors and samples to be sequenced can be found in <a href=https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2>this spreadsheet</a>.
In order to get access to this spreadsheet, please contact...
Note that this spreadsheet does not reflect at any given moment the data on S3 and which samples have already been processed, therefore it may contain samples that are still not available on the S3 bucket.


### <a name="s3_structure"></a> S3 Directory Structure

* raw_columbia/
* raw_sanger/
* processed/
* harmonized/
* vision/

The first two directories are designated for raw data uploads (.fastq files).
The third and fourth directories are processed versions of the raw data. The `vision` directory will be discussed later under Data Visualization

In more detail, the `processed` directory includes a `h5ad` file for each sample. This file can be used for downstream analysis in single-cell analysis tools such as scvi-tools [ref] and Seurat[ref]. 
Throughout the life cycle of the Immune Aging project we expect to improve our data processing pipelines. We therefore version the processed data by structuring the directories accordingly and by including a `.log` file with each processed sample that describes in details the processing pipeline used, such as the aligner used, software versions etc. The log file is important to allow reproducibility at the long-term, as we expect out pipeline procedures to update from time to time, given new best practices in the field.

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
The exact differences and the complete information about each `h5ad` file can be found in the file's log. In addition, we further maintain a lookup table that provides a concise summary of each version. You can find this table <a href="...">here</a>[todo].

The `harmonized` directory includes a `h5ad` with harmonized data containing multiple samples. Since the samples in the Immune Aging project are not generated and processed at the same time, for every new incoming sample or a batch of samples we will pool all existing data at the specific point in time.

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

For example, given two files `data_name.harmonized.v1.Feb21.h5ad, data_name.harmonized.v1.April21.h5ad`, we can tell by the file names that both used version `v1` of the harmonization pipeline.

In addition to the information in the log files, we further maintain a lookup table that provides a concise summary of each version and time stamp; particularly, the table indicates which samples were harmonized in a given data file. You can find this table <a href="...">here</a>[todo].


## <a name="upload"></a> Data upload

### <a name="naming_conventions"></a> Data Naming Conventions

Need to follow naming conventions as follows...

Note that sample names should not include a dot in its name. Sample name should not share its prefix with any other sample on S3. For example,... 

Not following the above conventions will fail the upload, which is described next.


### <a name="upload_s3"></a> Upload to S3

The very first step for uploading data is to update <a href=https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2>this spreadsheet</a> with the new file names, following the naming conventions provided above.

Once the spreadsheet is updated, you can move on to the data upload.
We provide a script for uploading pre-aligned `.fastq` files.
The upload script `data_upload.sh` gets three arguments as an input as follows:
```
data_upload.sh <sample_prefix> <institute> <credentials_file>
```

Description:

* `<sample_prefix>` - the sample prefix, which is the `sample_name` part of the data files
formatted as as described in the previous subsection. In case the data files of the samples are not in the current working directory you should include a path to data files.
* `<institute>` - can either be set to `columbia` or `sanger`.
* `<credentials_file>` - path to the user-specific credentials file (see Prerequisites subsection under Data Access).

For example, in order to... , we run:
```
data_upload.sh my_sample columbia keys.sh
```
 
Finally, after uploading new data, please notify ... by email...  [so that we can manually initiate the processing pipeline]

Notes:

* `data_upload.sh` will not allow to upload data files that do not properly follow the naming convention; in such a case an error will be issued.
* Uploading more than one sample can be done by running script multiple times, once per each sample.


## <a name="download"></a> Data Download

Users with access permission can see the list of samples on the S3 bucket by logging in to...

It is also possible to use the AWS CLI...


## <a name="visualization"></a> Data visualization

### <a name="live_vision"></a> Live VISION Session

We provide a live VISION session of the latest version of the harmonized data.
VISION allows to... [ref]...
Every time a new sample or a batch of samples has been uploaded, the data will be processed and harmonized with all the existing samples in the project, followed by generating a new VISION session.

The most updated VISION sessions can be found here, separated by tissue:
* <a href="...">Lung</a> 
* <a href="...">Liver</a> 

### <a name="live_vision"></a> Local VISION Session

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


## <a name="data_processing"></a> Data Processing

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

### <a name="processing_sample"></a> Data Harmonization

The script `harmonization_pipeline.py`... would make sense to run on the node with the GPU (s130)


### <a name="generating_keys"></a> Generating AWS credentials files

Explanations...


---

## Authors & Maintainers of the Immune Aging Data Hub

Elior Rahmani <erahmani@berkeley.edu>

Galen Xing <gx2113@columbia.edu>

Allon Wagner <allonwag@berkeley.edu>

Nir Yosef <niryosef@berkeley.edu>
