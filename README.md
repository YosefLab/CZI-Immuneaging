# The Immune Aging Data Hub

This repository serves as the data hub of the Immune Aging project.
Here you can find all the resources and information needed for accessing the currently available processed data of the project, as well as a complete documentation of the data collection process, including instructions for data upload, data processing and other data management procedures.

This page did not answer your question? Please <a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/issues">open an issue</a> and label it as a 'question'.

## Table of contents
1. [Overview](#overview)
1. [Preliminaries](#preliminaries)
    1. [Setting up data access](#preliminaries_access)
    1. [The Immune Aging samples spreadsheet](#preliminaries_spreadsheet)
    1. [The Immune Aging Dashboard](#dashboard)
1. [Code organization](#code_organization)
    1. [Align library overview](#align_library_overview)
    1. [Process library overview](#process_library_overview)
    1. [Process sample overview](#process_sample_overview)
1. [Data Download](#download)
    1. [Directory structure on S3](#download_structure)
    1. [Downloading data via the AWS console](#download_console)
    1. [Downloading data via terminal](#download_terminal)
1. [Data upload](#upload)
    1. [File naming conventions](#upload_naming)
    1. [Upload to S3](#upload_upload)
1. [Data visualization](#visualization)
1. [Data Processing](#processing)
    1. [Prerequisites](#processing_prerequisites)
    1. [Aligning libraries](#aligning_libraries)
    1. [Processing libraries](#processing_libraries)
    1. [Processing samples](#processing_samples)
    1. [Integrating samples](#integrating_samples)
    1. [Updating the dashboard](#dashboard_update)
    1. [Reproducibility and version control](#version_control)
    1. [Sandbox environment](#sandbox_environment)
    1. [Job queue](#job_queue)
    1. [Generating job configs](#job_configs_generation)
    1. [Executing job queue jobs](#job_queue_execution)
---

## <a name="overview"></a> Overview

Here is an overview of the various components involved in this project:

<center><img src="figures/project_overview.svg?raw=true" width="100%;" alt="IA project overview"></center>

## <a name="preliminaries"></a> Preliminaries

### <a name="preliminaries_access"></a> Setting up data access

The project's data is stored on an Amazon Web Services (AWS) Simple Storage Service (S3) bucket. Whether you will be uploading data to the Immune Aging S3 bucket or just using it for downloading data, you will need to set up your access as follows:

1. Get an IAM user and credentials file:
You will need an Identity and Access Management (IAM) user account, which will be associated with a specific user in the Immune Aging S3 bucket. Your user account will allow you to access the S3 bucket through the <a href="https://911998420209.signin.aws.amazon.com/console">AWS console</a>.
<br /><br />
If you will be uploading data to the bucket or will be using terminal commands for downloading data (i.e., rather than using the AWS console for downloading; terminal download, which will be described later, allows to download large amount of files more conveniently), you will also need a user-specific credentials file that will allow you to access the S3 bucket via terminal.
<br /><br />
The Yosef group will manage the generation of both IAM users and credentials file for collaborators on the Immune Aging project. 
In order to set up an IAM user and receive credentials file please email Elior Rahmani (erahmani@berkeley.edu) and cc your PI for approval of your request (the project's PIs are Menna Clatworthy, Donna Farber, Muzz Haniffa, Joanne Jones, Peter Sims, Sarah Teichmann, Roser Vento, and Nir Yosef). Note that IAM accounts and credentials files will be issued with read access only for all users, except for designated data uploaders and data processors who will also get write (i.e., upload) access. In any case, **DO NOT SHARE THE CREDENTIALS FILE WITH ANYONE**.

1. Install AWS CLI (Optional for downloading data):

* Uploading raw data (.fastq files) is done via the AWS Command Line Interface (CLI), which can be downloaded from <a href="https://aws.amazon.com/cli/">here</a>. While you will not need to explicitly use the AWS CLI for data upload, the upload script (to be described later) requires it as a dependency.
* Downloading multiple data files becomes more streamlined by using the AWS CLI rather than the <a href="https://911998420209.signin.aws.amazon.com/console">AWS console</a>, which is a web-based graphical user interface. For more details see [Data Download](#download).

**Note:** Non-US based users that wish to download data via AWS CLI or upload data using the upload tool described below may need to manually configure their location prior to any download/upload. Specifically, if your download/upload fails and you receive an error "Could not connect to the endpoint URL" then log in to the <a href="https://911998420209.signin.aws.amazon.com/console">AWS console</a> with your account and properly set the location in the top right corner of the console to match your physical location (e.g., users located in London should set `eu-west-2`).
If the problem persists, run in terminal the command `aws configure`. You will then be prompted to enter your `AWS Access Key ID` and `AWS Secret Access Key` - enter those keys as they appear in your credentials file. When prompted to enter your `Default region name`, enter the region name as you specified in the AWS console (e.g., London users should enter `eu-west-2`).

<!--
**Granting temporary/one-off access:**
To share data with temporary members (like research assistants temporarily helping out) or for one-off access, email Galen (gx2113@columbia.edu) with the S3 folder you want to share (or upload data to) and he will generate a script that will automatically upload/download data when run. This is super easy for him to do (says Galen writing this), so don't hesitate to reach out. DO NOT SHARE YOUR USER SPECIFIC CREDENTIALS!
-->

### <a name="preliminaries_spreadsheet"></a> The Immune Aging samples spreadsheet

The project's donors and samples, as well as the raw metadata, can be found in the <a href="https://docs.google.com/spreadsheets/d/1XC6DnTpdLjnsTMReGIeqY4sYWXViKke_cMwHwhbdxIY/edit?usp=sharing_eip&ts=6054e1a2" target="_blank">AI samples Google Spreadsheet</a>.
You should never share this link with anyone outside the immune aging project. Also, **DO NOT MAKE EDITS IN THIS FILE**, unless you are a designated data uploader and you read all the instructions under [Data upload](#upload) - the data processing scripts rely on the data in this spreadsheet and some changes can compromise data integrity.

Note that:

* This spreadsheet does not reflect at any given moment the existing data on the S3 bucket or which samples have already been processed. Specifically, it may include information about donors and samples that either were not sequenced yet, are being currently sequenced, being uploaded, or pending data processing.
* The metadata of the donors and samples in the spreadsheet is also available in the processed data files (see below).

### <a name="dashboard"></a> The Immune Aging Dashboard

We provide this <a href= "https://docs.google.com/spreadsheets/d/1v1TT1uJHTAxemnPmORX8XvSDZCXWq27ZK0Nbfugjts0/edit?usp=sharing">dashboard</a>, which reflects the current state of the data. The dashboard includes the following three sheets:
1. Tissue Coverage - describes the coverage of tissues across donors.
2. Samples QC - provides basic quality control summaries for each processed sample, as well as detailed warnings and/or errors that were issued during processing.
3. Tissue Integration - shows basic information on the current tissue-level integrations in the project, including the sample IDs and ages of the samples that were used for generating the data for each tissue. The last two columns provide links to the data files and to figures that demonstrate the integration through multiple UMAPs colored by different metadata and variables of the experimental design. Note that these links are accessible only after logging into the AWS console ([see instructions here](#preliminaries_access)).

## <a name="code_organization"></a> Code organization
There are 4 main classes of scripts that comprise our codebase:
- Library alignment
- Library processing
- Sample processing
- Sample integration
- Bash script and config file generation scripts
- Utilities

### <a name="align_library_overview"></a> Align library overview
- Align and count using cellranger count for GEX and GEX+ADT and GEX+HTO and cellranger vdj for BCR and TCR

### <a name="process_library_overview"></a> Process library overview
- If GEX:
    - Start with the aligned library h5ad adata file
    - Add library alignment metrics (from cellranger count) as .obs columns in the adata
    - Check the library against a list of known poor quality libraries and if so terminate execution.
    - Append library name to the cell barcodes to unique-ify them.
    - For each tissue, download the csv file containing its non-immune cell barcodes and remove them from adata.
    - Move any features in adata.X that are HTO counts to adata.obs.
    - Move any features in adata.X that are protein counts to adata.obsm.
    - Switch protein names to their internal names defined in the protein panels.
    - Save control and non-control proteins in different obsm objects.
    - Filter cells that dont have less than 600 genes expressed.
    - Compute for each cell the percentage of total counts that come from mitochondrial genes (start with “MT-”) and discard cells that have >20%.
    - Compute for each cell the percentage of total counts that come from ribosomal genes (start with “RPS” or “RPL”). We currently dont filter cells based on this metric.
    - Exclude any mitochondrial genes and the “MALAT1” gene.
    - If applicable, use hashsolo to classify cells as Doublet or assign them to the correct sample, based on their HTO counts.
    - Save and upload the results to AWS.
- If BCR or TCR:
    - Start with the filtered_contig_annotations.csv file
    - Use scirpy.io.read_10x_vdj to read it as an adata object.
    - Add library alignment metrics (from cellranger vdj) as .obs columns in the adata.
    - Filter out cells that were called with low_confidence by cellranger vdj.
    - Filter out cells that are considered multichain or ambiguous by scipry (are likely doublets).
    - Validate that there are no non-IR cells, and that cells' receptor_type match the library type being processed. Else terminate execution.
    - Save and upload the results to AWS.

### <a name="process_sample_overview"></a> Process sample overview
- Download h5ad files of processed libraries (GEX, BCR, and TCR) from AWS.
- For each GEX library that is multiplexed (has Classification obs column), select cells that correspond to this sample, if any.
- Skip a GEX library if it has less than N cells for this sample. N = 50 if this is a jejunum sample otherwise else 200.
- This is something we will need later: For each library take all genes that are expressed in at least 5 cells then take the intersection of all such genes across all libraries. Call these solo_genes.
- Concatenate cells from all GEX libraries, call this adata_gex.
- Skip any BCR/TCR library whose corresponding GEX library was not selected above.
- Concatenate cells from all BCR libraries. Same for TCR. Make sure the concatenation appends the same “-x” suffix to cells that it did for the cells of the corresponding GEX lib.
- Filter out any cells that are both in the concatenated BCR and TCR libs.
- Concatenate the concatenated BCR and TCR libs together in an object called adata_ir.
- Merge adata_gex with adata_ir, discarding any cells in the latter that dont have a match in the former. Call this adata.
- Compute the per-cell protein coverage (total number of protein reads per cell across control and non-control proteins). If the median per-cell protein coverage across all cells is 0, disregard protein counts altogether.
- Run decontX for estimating contamination levels from ambient RNA. Run a Rscript command for this. If there is more than one GEX library, use the “batch” column from the GEX concatenation done earlier. decontX returns contamination_levels per cells (save it as an obs in adata) and decontaminated counts (save it in a layer in adata).
- Filter out all cells from adata whose total decontaminated counts (after rounding), is less than 100.
- Make a copy of adata, call it adata_copy. Have adata_copy.X be the rounded decontaminated counts.
- Keep only the genes in adata_copy that are in solo_genes.
- Filter out any genes from adata_copy that are categorized as VDJ genes (the full list of this is kept on AWS). This set of is expected to express high donor-level variability (per the property of the adaptive immune repertoire) that is not interesting to us in the context of this study. Thus we filter these out prior to HVG selection in adata_copy.
- Detect the top 3k highly variable genes in adata_copy using seurat_v3 flavor.
- Download all celltypist models specified in the config file. Run prediction on adata_copy using every specified celltypist model and store the results (predicted labels and majority voting) in adata.obs.
- Using a custom-made celltypist model specifically trained to detect red blood cells, detect any erythrocytes (“RBC” predicted label) in adata and filter them out.
- If we have proteins, train a TotalVI model on adata_copy. Store the resulting latent representation in adata.obxm[“X_totalVI”]. Use it to compute UMAP coordinates and store it in adata.obsm[“X_umap_totalvi”].
- Train a SCVI model on adata_copy. Store the resulting latent representation in adata.obxm[“X_scVI”]. Use it to compute UMAP coordinates and store it in adata.obsm[“X_umap_scVI”].
- Train a SOLO model (one per-batch, if needed) from the pre-trained SCVI model and use it to filter cells in adata to only those predicted by SOLO as “singlet”.
- Run PCA on adata_copy. Use it to compute UMAP coordinates and store it in adata.obsm[“X_umap_pca”].

## <a name="download"></a> Data Download

### <a name="download_structure"></a> Directory structure on S3

Whether you will be downloading data from the S3 bucket [directly from the AWS consol](#download_console) or [via terminal](#download_terminal), you will need to know the bucket's directory structure.

The root directory of the Immune Aging S3 bucket includes the following sub-directories:

📂 raw_columbia/  
👉 contain raw sequencing data from NY  
📂 raw_sanger/  
👉 contain raw sequencing data from UK  
📂 aligned_libraries/  
👉 contains the output of library alignment for each library in the project. An intermediate output of the pipeline.   
📂 processed_libraries/  
👉 contains the output of library processing for each library in the project. An intermediate output of the pipeline.  
📂 processed_samples/   
👉 contains a processed data file for each sample in the project. An intermediate output of the pipeline    
📂 **integrated_samples**/   
👉 includes integration of samples. This is the directory that contains the final output of the pipeline  
📂 scanvi_integrated_samples/    
👉 includes integration of samples using SCANVI  
📂 cell_filtering/  
👉 for each tissue, contains a list of cells that we decided to filter out due to being non immune cell types. We download and use these files during library processing. This also contains the list of poor quality libraries to exclude (poor_quality_libs.csv).  
📂 combined_lib_alignment_metrics/  
👉 contains, for each donor (and also one file for all_donors), some high level metrics for all processed GEX/IR libraries for that donor. These are pulled from cellranger output and aggregated here.  
📂 per-compartment-barcodes/  
👉 contains the list of manually curated barcodes for each compartment. Under archived/ I store the previous versions before replacing them with new ones.  
📂 annotated_objects/  
👉 used to store annotated data objects  
📂 scratch/  
👉 used as an unstructured "dropbox" for testing or temporary artifacts  
📂 unpublished_celltypist_models/  
👉 contains a celltypist model we use for RBC detection (not officially published by celltypist). During sample processing we download the model from this location.  
📂 vdj_genes/  
👉 contains a csv file with the list of VDJ genes to exclude. During library processing we download the csv from this location.  
📂 job_queue/  
👉 a job queue that is used during processing (more details under [Data Processing](#processing))  

Here is more detail about some of the main directories listed above that contain outputs from the final and intermediate steps of the pipeline:

The `aligned_libraries` directory stores data of aligned libraries and supplementary files related to the alignment. Particularly, it includes a `h5ad` file for each sequenced library (a file format which can be used for downstream analysis in single-cell analysis tools such as <a href="https://scvi-tools.org/">scvi-tools</a> and <a href="https://satijalab.org/seurat/">Seurat</a>).
Throughout the life cycle of the Immune Aging project we expect to improve our data processing pipelines. We therefore version each intermediate step of processed data by structuring the directories accordingly and by including a `.log` file with each processed data file that describes in details the processing pipeline used. The log file is important to guarantee reproducibility at the long-term, as we expect our pipeline to update from time to time, as new best practices in the field emerge.

The structure of the `aligned_libraries` directory is as follows:

* aligned_libraries/
    * configs/
        * align_library.v1.txt
        * align_library.v2.txt
        * ...
    * v1/
        * library1_gex
            * library1.cellranger.cloupe.cloupe
            * library1.cellranger.feature_ref.csv
            * library1.cellranger.libraries.csv
            * library1.cellranger.metrics_summary.csv
            * library1.cellranger.web_summary.html
            * library1.v1.h5ad
            * align_library.library1.v1.log
        * library2_ir
            * library2.cellranger.airr_rearrangement.tsv
            * library2.cellranger.all_contig_annotations.csv
            * library2.cellranger.all_contig_annotations.json
            * library2.cellranger.all_contig.fasta
            * library2.cellranger.filtered_contig_annotations.v1.csv
            * library2.cellranger.filtered_contig.fasta
            * library2.cellranger.metrics_summary.csv
            * library2.cellranger.vloupe.vloupe
            * library2.cellranger.web_summary.html
            * align_library.library2.v1.log
        * ...
    * v2/
        * ...
    * ...

The `configs` directory includes versioned configuration files for the alignment pipeline. For each version of the aligned data, one designated directory (e.g., directory `v1` for version 1) includes the alignment outputs for each library. Particularly, it includes a .h5ad file (or filtered_contig_annotations.csv file), output files from cellranger, and a .log file documenting the execution of the pipeline on the library. Note: libraries are not named arbitrarily as library1, library2 etc. but rather take the following naming convention: `<donor_id>_<seq_run>_<library_type>_<library_id>`.

The `processed_libraries` directory stores data of processed libraries (i.e., beyond alignment), an intermediate product before the sample-level processing; its structure is as follows:

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

Here, files with the prefix `process_library.configs.` include the configurations that were used in the execution of the library processing pipeline, and files with a `.log` suffix are documentation of the execution of the pipeline. As in the `aligned_libraries` directory, the actual library names follow the naming convention `<donor_id>_<seq_run>_<library_type>_<library_id>`.

The `processed_samples` directory stores data of processed samples (i.e., sample-level integration across different libraries); its structure is as follows:

* processed_samples/
    * sample1/
        * v1/
            * sample1.processed.v1.h5ad
            * sample1.processed.v1.decontx_model.RData
            * sample1.processed.v1.scvi_model.zip
            * process_sample.sample1.v1.log
            * process_sample.configs.sample1.v1.txt
        * v2/
            * sample1.processed.v2.h5ad
            * sample1.processed.v2.decontx_model.RData
            * sample1.processed.v2.scvi_model.zip
            * process_sample.sample1.v2.log
            * process_sample.configs.sample11.v2.txt
        * ...
    * sample2/
        * ...
    * ...

Here, files with the prefix `process_sample.configs.` include the configurations that were used in the execution of the sample processing pipeline, and files with a `.log` suffix are documentation of the execution of the pipeline. 
Actual sample names follow the naming convention `<sample_id>_<data_type>`.

**Note**:
In some cases, the sub-directories under `processed_libraries` and `processed_samples` may not contain data files (i.e., h5ad files) due to a failure of the processing pipeline to run on a specific library or sample. In such cases, there will still be a log file that should include some error message at the end.

The `integrated_samples` directory stores integrated data. Particularly, integration of multiple samples for every tissue in the project. The structure of the directory is as follows:

* integrated_samples/
    * tissue_level/
        * tissue1/
            * v1/
                * tissue1.v1.h5ad
                * tissue1.v1.gene_stats.csv
                * tissue1.v1.celltypist_dotplots.zip
                * tissue1.v1.scvi_model_batch_key_donor_id.zip
                * tissue1.v1.totalvi_model_batch_key_donor_id.zip
                * integrate_samples.tissue1.v1.log
                * integrate_samples.configs.tissue1.v1.txt
                * [same as above for stim and unstim segments of the tissue if applicable]
        * tissue2/
        * ...
    * compartment_level/
        * compartment1/
            * v1/
                * compartment1.v1.h5ad
                * compartment1.v1.gene_stats.csv
                * compartment1.v1.celltypist_dotplots.zip
                * compartment1.v1.scvi_model_batch_key_donor_id.zip
                * compartment1.v1.totalvi_model_batch_key_donor_id.zip
                * compartment1.v1.scvi_model_batch_key_donor_id+tissue.zip
                * compartment1.v1.totalvi_model_batch_key_donor_id+tissue.zip
                * integrate_samples.compartment1.v1.log
                * integrate_samples.configs.compartment1.v1.txt
                * [same as above for stim and unstim segments of the compartment if applicable]
        * compartment2/
        * ...

The `.h5ad` file is the integrated data file, the `.zip` files incudes the scVI and totalVI models, the `.txt` file includes the configurations that were used in the execution of the sample integration (i.e., `integrate_samples.py`), and the `.log` file provides a documentation of the execution of the integration. 

### <a name="download_console"></a> Downloading data via the AWS console

Data can be downloaded from the Immune Aging S3 bucket by logging in to AWS through <a href="https://911998420209.signin.aws.amazon.com/console.">this link</a> and navigating through the project's directory structure.

### <a name="download_terminal"></a> Downloading data via terminal

To read from the S3 bucket via terminal:
1. Run your user-specific credentials file to set the aws keys (note that the credential file is a shell file, which can be executed in a linux/mac environment).
2. Sync directories to your local machine via `aws s3 sync <source> <target> [--options]` (where source is the aws folder). You can also use `aws s3 ls <target> [--options]` to list contents of a directory. Check out more commands <a href= "https://docs.aws.amazon.com/cli/latest/userguide/cli-services-s3-commands.html">here</a>.

## <a name="upload"></a> Data upload

### <a name="upload_naming"></a> File naming conventions

New data can only be uploaded as pre-aligned, gzip-compressed `fastq` files (i.e., `.fastq.gz` files). Each file must follow the following naming convention:
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

First, download the <a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/blob/main/data_management/scripts/upload.py"> upload.py</a> script (open the link, right click on the screen, and save as "upload.py"). Please make sure you always use the most updated version of the script. Second, make sure you have a working installation of python with the `pandas` library installed. If you do not have python with `pandas` already installed, you can simply download and install <a href="https://www.anaconda.com/products/individual">Anaconda</a>, a distribution of python that includes many popular libraries for data science, including `pandas`.

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

Processed data can be easily visualized using <a href="https://chanzuckerberg.github.io/cellxgene/">cellxgene</a>, which provides a browser-based user-interface for basic data exploration. Currently, only our sample-level data files can be visualized by cellxgene (i.e., .h5ad files under the `processed_samples` directory on the S3 bucket).

In order to run cellxgene, first install it via terminal by running:

``pip install cellxgene``

Then, after downloading a processed .h5ad file from the S3 bucket (see [Data Download](#download)) you can open a cellxgene session with the data by running:

``cellxgene launch filename.h5ad --open``

## <a name="data_processing"></a> Data Processing

The processing scripts can be found in this repository, <a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/tree/main/data_processing/scripts">here</a>. It is advised that you clone the repository and make sure to pull updates before working on new data.

Our pipeline is largely composed of a body of scripts, each run in sequence depending on the stage of the pipeline that is being executed. As a reminder, here are the 4 stages involved:
- Library alignment
- Library processing
- Sample processing
- Sample integration
- [Recommended] Updating the dasboard

### <a name="processing_prerequisites"></a> Prerequisites

First, Download anaconda with python >= 3.7.
Second, get the .yml file of the latest version of the python environment from <a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/tree/main/envs">here</a> (i.e., immune_aging.py_env.*.yml where * is the latest version available), and use it to set up a new conda environment via terminal.

For example:
```
conda env create -f immune_aging.py_env.v1.yml
```
This will create a new environment under `$CONDA_PREFIX"/envs/immune_aging.py_env.v1.yml"`. This can take a while (1hr+).

Activate the environment by typing:
```
source activate immune_aging.py_env.v1
```
and then complete the setup of the environment by opening an R session and running the commands in the latest version of the R setup script (`immune_aging.R_setup.v*.R`), which can be found <a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/tree/main/envs">here</a>.

Note that every time you wish to run the processing scripts you need to activate the environment. At the end of the execution you can deactivate the environment by typing
```
conda deactivate
```

### <a name="aligning_libraries"></a> Aligning libraries

Once the fastq files of new data from a donor have been uploaded to the S3 bucket and the Google Spreadsheet has been properly updated to include all metadata about the donor and the associated samples we can start the library alignment. The script `generate_library_alignment_script.py` generates a shell script that can be used to execute the library alignment script `align_library.py` on all the libraries associated with the donor.

Note that the script `generate_library_alignment_script.py` requires a configuration file, which will be used for the execution of `align_library.py`.
<a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/blob/main/data_processing/configs_templates/align_library.configs_file.example.txt">Here</a> you can find a template for generating such a configuration file, and <a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/blob/main/data_processing/configs_templates/align_library_configs.md">here</a> you can find a description of each of the configuration fields.

⚙️ How to align libraries  
You don't have to be on a machine with  a GPU for this stage.
1. Start a tmux session
2. Activate your IA environment
3. Create alignment configs text file
4. Generate alignment script (bash script) by running:  
`python generate_library_alignment_script.py <donor_id> <seq_run> <path to your data_processing/scripts code> <path to your configs.txt file> <output dir where to place alignment results>`
5. Run the bash script

### <a name="processing_libraries"></a> Processing libraries

The script `process_library.py` runs on a single aligned library. All the parameters for library processing are defined in a configuration file. Once you set up a configuration file you can run the library processing script as follows:

```python
python process_library.py configs.txt
```
<a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/blob/main/data_processing/configs_templates/process_library.configs_file.example.txt">Here</a> you can find a template for generating configuration files for `process_library.py`, and <a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/blob/main/data_processing/configs_templates/process_library_configs.md">here</a> you can find a description of each of the fields in the configuration file.

⚙️ How to process libraries  
You don't have to be on a machine with  a GPU for this stage.
1. Locally create a directory called process_configs
1. Generate processing library config files -- for output destination use the absolute path to the process_configs dir above  
`python generate_processing_config_files.py library <path to your data_processing/scripts code>  <path to your process_configs directory> <donor_id> <seq_run> <path to your AWS_immuneaging_credentials_file.sh> v1 v1`
1. Make sure on AWS the job queue folders are all empty
1. Upload the results of generate configs script to AWS
1. Locally create a directory called processing_results
1. Generate process library bash script  
`python generate_processing_scripts.py <path to your AWS_immuneaging_credentials_file.sh> <path to your processing_results directory> <path to your data_processing/scripts code> <output dir where to place the generated bash file>`
1. Make sure you've done the needful to avoid overwriting the existing version AWS (by having a different config file than the on on AWS) if that is your intention
1. Start a tmux session
1. Run the generated bash script
1. Once done, run the digest process library logs script against the generated process library logs:  
`python digest_logs.py print_digest library <donor_id> <seq_run> aws <version of the processed libs> <working dir for this script> <path to your AWS_immuneaging_credentials_file.sh>`

### <a name="processing_samples"></a> Processing samples

The script `process_sample.py` processes data for a single sample (i.e., specific donor, tissue, and stimulation/non-stimulation). For samples that were sequenced using multiple libraries this script integrates the data from the different libraries. All the parameters for the different steps are defined in a configuration file. Once you set up a configuration file you can run the sample processing script as follows:

```python
python process_sample.py configs.txt
```
<a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/blob/main/data_processing/configs_templates/process_sample.configs_file.example.txt">Here</a> you can find a template for generating configuration files for `process_sample.py`, and <a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/blob/main/data_processing/configs_templates/process_sample_configs.md">here</a> you can find a description of each of the fields in the configuration file.

⚙️ How to process samples  
You DO have to be on a machine with  a GPU for this stage.  
1. Locally create a directory called process_configs
1. Generate processing sample config files -- for output destination use the absolute path to the process_configs dir above  
`python generate_processing_config_files.py sample <path to your data_processing/scripts code> <path to your process_configs directory> <donor_id> <seq_run> <path to your AWS_immuneaging_credentials_file.sh> <latest processed gex lib version for this donor> <latest processed ir lib version for this donor>`
1. Make sure on AWS the job queue folders are all empty
1. Upload the results of generate configs script to AWS
1. Locally create a directory called processing_results
1. Generate process sample bash script  
`python generate_processing_scripts.py <path to your AWS_immuneaging_credentials_file.sh> <path to your processing_results directory> <path to your data_processing/scripts code> <output dir where to place the generated bash file>`
1. Make sure you've done the needful to avoid overwriting the existing version on AWS (by having a different config file than the on on AWS) if that is your intention
1. Start a tmux session
1. Run the generated bash script
1. Once done, run the digest process sample logs script against the generated process sample logs:
`python digest_logs.py print_digest sample <donor_id> <seq_run> aws <version of the processed samples> <working dir for this script> <path to your AWS_immuneaging_credentials_file.sh>`
1. Update the dashboard using the output of the below (you still need to copy and paste it manually in the dashboard):  
`python digest_logs.py get_csv sample <donor_id> <seq_run> aws <version of the processed samples> <working dir for this script> <path to your AWS_immuneaging_credentials_file.sh>`

### <a name="integrating_samples"></a> Integrating samples

The script `integrate_samples.py` can be used for integrating a specified list of processed samples. More specifically, we use it for integrating all samples of a given tissue or compartment, thus creating tissue-level or compartment-level integration of the data. This script requires a configuration file - <a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/blob/main/data_processing/configs_templates/integrate_samples.configs_file.example.txt">here</a> you can find a template for generating such a configuration file, and <a href="https://github.com/YosefLab/Immune-Aging-Data-Hub/blob/main/data_processing/configs_templates/integrate_samples_configs.md">here</a> you can find a description of each of the configuration fields.

In addition, the script `generate_integration_config_files_and_scripts.py` automatically generates configuration files for `integrate_samples.py` - one per tissue or compartment - by collecting for each tissue/compartment the list of available processed samples. Note that for a compartment, we simply grab all samples at this stage. Later, during the execution of the integrate_samples script, we only select cells from each sample that belong to that compartment.

⚙️ How to integrate samples compartment-wise   
You DO have to be on a machine with  a GPU for this stage.
1. Locally create a directory called process_configs
1. Generate integration config files and script (though we won't use the script below) -- for output destination use the absolute path to the process_configs dir above  
`python generate_integration_config_files_and_script.py <path to your data_processing/scripts code> <path to your process_configs directory> <path to your AWS_immuneaging_credentials_file.sh> compartment`
1. Start a tmux session
1. Integrate all compartments. For example for the B compartment:  
`python integrate_samples.py process_configs/integrate_samples.B.configs.txt`

⚙️ How to integrate samples tissue-wise   
You DO have to be on a machine with  a GPU for this stage.
1. Locally create a directory called process_configs
1. Generate integration config files and script (though we won't use the script below) -- for output destination use the absolute path to the process_configs dir above  
`python generate_integration_config_files_and_script.py <path to your data_processing/scripts code> <path to your process_configs directory> <path to your AWS_immuneaging_credentials_file.sh> tissue`
1. Start a tmux session
1. Integrate all tissues. For example for the spleen tissue:   
`python integrate_samples.py /data/yosef2/scratch/immuneaging/vv_runs/process_configs/integrate_samples.SPL.configs.txt`
1. Update the dashboard using the output of the below (you still need to copy and paste it manually in the dashboard):   
`python dashboard_utils.py tissue_integration_results <working dir for this script> <path to your AWS_immuneaging_credentials_file.sh> <version of integrated tissues> True`

### <a name="version_control"></a> Reproducibility and version control

Our pipeline is versioned via:
- built-in versioning on the AWS S3 bucket
- an ad-hoc versioning mechanism we implemented as follows: Each stage of the pipeline (library alignment, processing, sample integration, etc.) is associated with a configuration file that drives the parameters of that step. A snapshot of this config file is uploaded and saved on S3. The "version" of a run is essentially the version of the config file for that run (for each part of the pipeline). Outputs of that stage are also stamped with the version of the run (for each run).

In some scripts you will also see a variable called "VARIABLE_CONFIG_KEYS". This indicates parameters that can be changed in the config file without causing the system to create a new version (for example name of the person who ran the script).

Adding keys to config files:
If adding keys that should not affect versioning then update VARIABLE_CONFIG_KEYS in the processing files. Also, to ensure backwards compatiblity, make sure that when you reference the new key in the code, you account for the scenario where it may not exist.

### <a name="sandbox_environment"></a> Sandbox environment

The configuration files for running `process_library.py` and `process_sample.py` include a `sandbox_mode` argument. Setting this argument to `"True"` indicates that outputs should not be uploaded to the S3 bucket. The sandbox environment allows data owners to experiment with different configurations. Once the configurations were tuned and reported by the data owner as appropriate (see next subsection), a system admin can run the processing while setting `sandbox_mode` to `"False"`.

It is likely that during data curation data owners will find bugs and/or will have suggestions for implementing additional/different logics in the processing scripts. You can either post an issue or make a pull request with your suggestions. If you are suggesting any changes please bear in mind that any updates will have to maintain backwards compatibility in order to guarantee reproducibility of previous versions of the processed data.

## <a name="job_queue"></a> Job queue

Once a data owner makes a final decision about configurations for the processing of specific libraries and samples, the final configuration files should be uploaded to the S3 bucket through the <a href="https://911998420209.signin.aws.amazon.com/console">AWS console</a>. Specifically, configurations for processing libraries should be uploaded to `s3://immuneaging/job_queue/process_library/` and configurations for processing samples should be uploaded to `s3://immuneaging/job_queue/process_sample/`. Once configuration files are uploaded to these directories they are considered as jobs to be executed, and the outputs of the processing will be saved, stamped with a version, and become viewable via the S3 bucket to everyone with data access in the project.

Once an admin starts executing the jobs in the queue (see [Executing job queue jobs](#admins_job_queue_execution)), the configuration files will be removed from their original directories and will be moved to 
`s3://immuneaging/job_queue/process_library.running/` and `s3://immuneaging/job_queue/process_sample.running/` to indicate their status.

**NOTE**:
At the moment, the system does not notify the admins about new jobs in the queue. If you upload new jobs please notify Elior Rahmani by email (erahmani@berkeley.edu).

### <a name="job_configs_generation"></a> Generating job configs
A job can be of two types: library processing or sample processing. Each job is characterized by a config file that describes how the job needs to run. The script `generate_processing_config_files` can be used to generate config files for all process_library and process_sample jobs associated with a given donor. It can be run as follows:

```
python generate_processing_config_files.py <config_type> <code_path> <output_destination> <donor_id> <seq_run>
```

where config_type can be one of "library", "sample" or "all".

After generating the config files, we can run the following commands to upload them to the job queue on AWS (after setting the AWS credentials as env variables):

```
aws s3 sync config_files s3://immuneaging/job_queue/process_library/ --exclude "*" --include "process_library*.configs.txt"
```

```
aws s3 sync config_files s3://immuneaging/job_queue/process_sample/ --exclude "*" --include "process_sample.configs*.txt"
```

assuming config_files is the directory containing the config files.

Once these jobs are generated and uploaded to AWS, we can proceed to executing them, see [Executing job queue jobs](#admins_job_queue_execution).

### <a name="job_queue_execution"></a> Executing job queue jobs
Once the job queue on AWS (see [Job queue](#job_queue)) has jobs to run, we need to generate processing scripts for each type of job (sample processing or library processing) and execute them. The script `generate_processing_scripts.sh` helps automate this process. It can be run as follows:

```
python generate_processing_scripts.py <aws_credentials_file> <output_dir> <code_path> <output_path>
```

The script syncs the job queue from AWS down to the local machine. It then crawls the queued jobs (each characterized by their corresponding config files) for each type of job: process_library or process_sample. For each job (each config file), it adds the commands needed to execute that job to a list of commands keyed by donor. In the end, we end up with a set of shell script files each specific to a [donor, job_type] tuple. Each script file contains all of the commands needed to execute the jobs associated with that donor and job_type. These commands are, for example: activate the specified conda environment, execute the process_sample python script at the given code_path, deactivate the conda environment, etc. The script also updates the remote job queue on AWS to indicate which jobs are running.

Once all the shell script files are generated, you can execute them in the shell to kick off the execution of each corresponding job. Make sure to execute all process_library jobs before executing any process_sample jobs.

**Note:** `output_dir` specifies the location for the processing outputs and `output_path` specifies the location to which the shell file with the run commands will be saved.

---

## Admins of the Immune Aging Data Hub

Elior Rahmani <erahmani@berkeley.edu>
Valeh Valiollah Pour Amiri <valehvpa@berkeley.edu>

## Authors of the Immune Aging Data Hub

Elior Rahmani, Valeh Valiollah Pour Amiri, Galen Xing, Allon Wagner, and Nir Yosef
