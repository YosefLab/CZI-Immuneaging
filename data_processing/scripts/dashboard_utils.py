# This script can be used to generate various metrics that can be used to produce our data dashboard.
# - To get tissue coverage csv, run as follows:
#   python dashboard_utils.py tissue_coverage
# - To get tissue integration results csv, run as follows:
#   python dashboard_utils.py tissue_integration_results <working_dir> <s3_access_file> <rm_working_dir>
#   where <rm_working_dir> is "True" or "False" for indicating whether the working dir should be removed at the end; this argument is optional (default value is "True").

import sys
import io
import csv
import numpy as np
import pandas as pd
import os
import traceback
import scanpy as sc
import zipfile
import shutil

import utils
from logger import RichLogger

sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
logger = RichLogger()

def get_tissue_coverage_csv():
    logger.add_to_log("Downloading the Donors sheet from Google drive...")
    donors_csv = utils.read_immune_aging_sheet("Donors")
    donors = donors_csv["Donor ID"]
    logger.add_to_log("Downloading the Samples sheet from Google drive...")
    samples = utils.read_immune_aging_sheet("Samples")
    tissues = np.unique(samples["Organ"])

    csv_rows = []
    for donor in donors:
        index = samples["Donor ID"] == donor
        organs = samples["Organ"][index] # organ and tissue are synonym
        stimulations = samples["Stimulation"][index]

        # build the dictionary of tissue to a list of tissue stim status
        # a tissue can have more than one stim status (e.g. it can be nonstim and stim)
        tissue_stimulation = {tissue: [] for tissue in tissues}
        for idx,_ in organs.items():
            stim_status = tissue_stimulation[organs[idx]]
            # if for this donor, this tissue was so far found to have no data then set it to
            # stimulations[i]. If it was already reported to be stimulated with stimulations[i]
            # then no need to do anything. Otherwise, append stimulations[i] to whatever we have
            # so far
            if stim_status == [] or stimulations[idx] not in stim_status:
                stim_status.append(stimulations[idx])
            tissue_stimulation[organs[idx]] = stim_status

        # fill in the row and add it to the list of rows
        csv_row = {}
        for tissue in tissues:
            tissue_stimulation[tissue].sort()
            stim_status = tissue_stimulation[tissue]
            csv_row[tissue] = "No Data" if stim_status == [] else ";".join(stim_status)
        csv_rows.append(csv_row)

    # write the csv
    csv_file = io.StringIO()
    writer = csv.DictWriter(csv_file, fieldnames=tissues)
    writer.writeheader()
    writer.writerows(csv_rows)
    print(csv_file.getvalue())
    csv_file.close()

def generate_tissue_integration_figures(adata, tissue, version, working_dir, base_aws_url, base_s3_url, base_s3_dir):
    
    # add two age fields for the umaps - age as a quantitative variable and as a categorical variable using 3 age groups
    ages_str = np.unique(adata.obs["age"])
    ages = {}
    for a in ages_str:
        try:
            ages[a] = float(a)
        except ValueError:
            # this is a range and not a number
            b = sorted([float(j) for j in a.split("-")])
            ages[a] = b[0]+(b[1]-b[0])/2 # take the middle of the range
    adata.obs["age_continuous"] = np.array([ages[j] for j in adata.obs["age"].values])
    age_categorical = []
    grp1_label = "{}-39 y/o".format(round(np.min(adata.obs["age_continuous"].values)))
    grp2_label = "40-59 y/o"
    grp3_label = "60-{} y/o".format(round(np.max(adata.obs["age_continuous"].values)))
    for j in adata.obs["age_continuous"].values:
        if j <= 40:
            age_categorical.append(grp1_label)
        elif j <= 60:
            age_categorical.append(grp2_label)
        else:
            age_categorical.append(grp3_label)
    adata.obs["age_categorical"] = pd.Categorical(np.array(age_categorical))

    fields_to_plot = ["site", "donor_id", "sample_id", "total_counts", "sex", "age_continuous", "age_categorical", "stimulation"]
    vmins = ["p0", "p0", "p0", "p10", "p0", "p0", "p0", "p0",] # minimal value for sc.pl.umap (p0 means using the minimal value in the data; p10 means the 10th percentile)
    vmaxs = ["p100", "p100", "p100", "p90", "p100", "p100", "p100", "p100",] # maximal value for sc.pl.umap (p100 means using the maximal value in the data)
    min_cell_type_frac = 0.001 # for generating annotation figures (cells annotated as coming from cell types that demonstrate less than min_cell_type_frac fraction of cells in the data will be removed form the plots)
    def extract_abundant_cell_types(adata, min_cell_type_frac, obs_key):
        # for extracting only cells that were annotated as a cell type that was used for at least min_num_cells_per_ct cells
        cell_types = np.unique(adata.obs[obs_key])
        keep = []
        for cell_type in cell_types:
            if np.sum(adata.obs[obs_key] == cell_type) >= round(min_cell_type_frac*adata.n_obs):
                keep.append(cell_type)    
        return(adata.obs[obs_key].isin(keep))

    logger.add_to_log("Generating figures for tissue: {}, version: {}\n".format(tissue,version))
    prefix = "{}.{}".format(tissue, version)
    
    figures_urls = []

    # plot umaps and color by metadata; plot all umap versions (i.e. based on all dimensionality reductions we have)
    umap_keys = ["X_umap_scvi_integrated","X_umap_pca"]
    integration_models = ["scvi","pca"]
    if "X_umap_totalvi_integrated" in adata.obsm:
        umap_keys.append("X_umap_totalvi_integrated")
        integration_models.append("totalvi")
    
    for mdl in range(len(umap_keys)):
        adata.obsm["X_umap"] = adata.obsm[umap_keys[mdl]]
        integration_model = integration_models[mdl]
        figures_dir = os.path.join(working_dir,"{}.{}.figures".format(prefix,integration_model))
        os.system("mkdir -p " + figures_dir)
        for field_index in range(len(fields_to_plot)):
            field = fields_to_plot[field_index]
            # looks like sc.pl.umap cannot take an absolute path; move files from the output dir of sc.pl.umap to the working directory
            filename = ".{}.{}.{}.pdf".format(prefix,integration_model,field)
            sc.pl.umap(adata, color=[field], save = filename, vmin = vmins[field_index], vmax = vmaxs[field_index])
            filename = "umap" + filename
            shutil.move("figures/" + filename, os.path.join(figures_dir, filename))

        # for each celltypist model used plot a umap colored by the predicted cell types; remove lowly abundant cell types
        celltypist_mdl = 1
        while "celltypist_model.{}".format(str(celltypist_mdl)) in adata.obs:
            celltypist_mdl_name = "celltypist_{}".format(adata.obs["celltypist_model." + str(celltypist_mdl)][1].split("/")[-1].split(".")[0])
            filename = ".{}.{}.{}.pdf".format(prefix,integration_model,celltypist_mdl_name)
            obs_key = "celltypist_predicted_labels.{}".format(str(celltypist_mdl))
            abundant_cell_types = extract_abundant_cell_types(adata, min_cell_type_frac, obs_key)
            sc.pl.umap(adata[abundant_cell_types,], color=[obs_key], save = filename)
            filename = "umap" + filename
            shutil.move("figures/" + filename, os.path.join(figures_dir, filename))
            celltypist_mdl += 1

        logger.add_to_log("Saving .zip file with the figures...")
        zip_filename = "{}.{}.figures.zip".format(prefix,integration_model)
        figures_file_path = os.path.join(working_dir, zip_filename)
        zipf = zipfile.ZipFile(figures_file_path, 'w', zipfile.ZIP_DEFLATED)
        utils.zipdir(figures_dir, zipf)
        zipf.close()
        
        logger.add_to_log("Uploading {} to aws...".format(zip_filename))
        aws_destination = "{}/{}/{}/{}".format(base_s3_url, base_s3_dir, tissue, version)
        sync_cmd = 'aws s3 sync --no-progress {} {} --exclude "*" --include {}'.format(working_dir, aws_destination, zip_filename)
        logger.add_to_log("sync_cmd: {}".format(sync_cmd))
        logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
        
        shutil.rmtree(figures_dir)
        figures_urls.append("{}{}/{}/{}/{}".format(base_aws_url, base_s3_dir, tissue, version, zip_filename))
    
    shutil.rmtree("figures")
    # concatenate urls
    figures_url = "[{}]".format(" ".join(figures_urls))
    return figures_url

def get_tissue_integration_results_csv(working_dir: str, s3_access_file: str, rm_working_dir: bool):
    logger.add_to_log("Downloading the Samples sheet from Google drive...")
    samples = utils.read_immune_aging_sheet("Samples")
    tissues = np.unique(samples["Organ"])

    # define the csv headers
    CSV_HEADER_TISSUE: str = "Tissue"
    CSV_HEADER_DONOR_COUNT: str = "# Donors"
    CSV_HEADER_DONORS: str = "Donors"
    CSV_HEADER_DONORS_AGE: str = "Donors age"
    CSV_HEADER_CELL_COUNT: str = "# Cells"
    CSV_HEADER_SAMPLE_COUNT: str = "# Samples"
    CSV_HEADER_SAMPLES: str = "Samples"
    CSV_HEADER_ANNDATA: str = "Anndata URL"
    CSV_HEADER_FIGURES: str = "Figures URL"

    # prepare the needful for downloading data from aws
    utils.set_access_keys(s3_access_file)
    #BASE_AWS_URL = "https://immuneaging.s3.us-west-1.amazonaws.com"
    BASE_AWS_URL = "https://s3.console.aws.amazon.com/s3/object/immuneaging?prefix="
    BASE_S3_URL = "s3://immuneaging"
    BASE_S3_DIR = "integrated_samples/tissue_level"
    version = "v1" # TODO add ability to get the latest version (github issue #32)

    csv_rows = []
    for tissue in tissues:
        # create a csv row and initialize it
        csv_row = {}
        csv_row[CSV_HEADER_TISSUE] = tissue
        csv_row[CSV_HEADER_DONOR_COUNT] = -1
        csv_row[CSV_HEADER_DONORS] = []
        csv_row[CSV_HEADER_DONORS_AGE] = []
        csv_row[CSV_HEADER_CELL_COUNT] = 0
        csv_row[CSV_HEADER_SAMPLE_COUNT] = -1
        csv_row[CSV_HEADER_SAMPLES] = []
        csv_row[CSV_HEADER_ANNDATA] = "Pending processing"
        csv_row[CSV_HEADER_FIGURES] = "Pending processing"

        # download the integrated anndata from aws
        try:
            file_name = "{}.{}.h5ad".format(tissue, version)

            # first check to see if we have an anndata for this tissue
            ls_cmd = "aws s3 ls {}/{}/{}/{}/".format(BASE_S3_URL, BASE_S3_DIR, tissue, version)
            files = os.popen(ls_cmd).read()
            logger.add_to_log("aws response: {}\n".format(files))
            found = False
            for f in files.rstrip().split('\n'):
                if f.split(' ')[-1] == file_name:
                    found = True
                    break

            if not found:
                logger.add_to_log("Integrated annotated data not found for tissue {}".format(tissue))
            else:
                csv_row[CSV_HEADER_ANNDATA] = "{}{}/{}/{}/{}".format(BASE_AWS_URL, BASE_S3_DIR, tissue, version, file_name)
                dir_name = "{}/{}/{}/{}/".format(BASE_S3_URL, BASE_S3_DIR, tissue, version)
                sync_cmd = 'aws s3 sync --no-progress {} {} --exclude "*" --include {}'.format(dir_name, working_dir, file_name)
                logger.add_to_log("sync_cmd: {}".format(sync_cmd))
                logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
                
                # now extract info from the adata
                h5ad_file = os.path.join(working_dir, file_name)
                adata = sc.read_h5ad(h5ad_file)
                csv_row[CSV_HEADER_CELL_COUNT] = adata.n_obs
                csv_row[CSV_HEADER_SAMPLES] = np.unique(adata.obs["sample_id"])
                csv_row[CSV_HEADER_SAMPLE_COUNT] = len(csv_row[CSV_HEADER_SAMPLES])
                donor_ids = np.unique(adata.obs["donor_id"])
                csv_row[CSV_HEADER_DONORS] = donor_ids
                csv_row[CSV_HEADER_DONOR_COUNT] = len(csv_row[CSV_HEADER_DONORS])
                ages = []
                for donor_id in donor_ids:
                    ages.append(np.unique(adata[adata.obs["donor_id"] == donor_id,].obs["age"])[0])
                csv_row[CSV_HEADER_DONORS_AGE] = np.array(ages)
                
                # generate figures
                csv_row[CSV_HEADER_FIGURES] = generate_tissue_integration_figures(adata, tissue, version, working_dir, BASE_AWS_URL, BASE_S3_URL, BASE_S3_DIR)

                if rm_working_dir:
                    # clean up
                    os.system("rm -r {}/*".format(working_dir))
        except Exception:
            logger.add_to_log("Execution failed for tissue {} with the following error:\n{}".format(tissue, traceback.format_exc()), "critical")
            if rm_working_dir:
                os.system("rm -r {}/*".format(working_dir))
            logger.add_to_log("Continuing execution for other tissues...")
            continue
        if found:
            csv_rows.append(csv_row)

    # write the csv
    csv_file = io.StringIO()
    field_names = [
        CSV_HEADER_TISSUE,
        CSV_HEADER_DONOR_COUNT,
        CSV_HEADER_DONORS,
        CSV_HEADER_DONORS_AGE,
        CSV_HEADER_CELL_COUNT,
        CSV_HEADER_SAMPLE_COUNT,
        CSV_HEADER_SAMPLES,
        CSV_HEADER_ANNDATA,
        CSV_HEADER_FIGURES,
    ]
    writer = csv.DictWriter(csv_file, fieldnames=field_names)
    writer.writeheader()
    writer.writerows(csv_rows)
    print(csv_file.getvalue())
    csv_file.close()


action = sys.argv[1]
assert(action in ["tissue_coverage", "tissue_integration_results"])
if action == "tissue_coverage":
    get_tissue_coverage_csv()
else:
    working_dir = sys.argv[2] # must be the absolute path to the local directory where we will place artifacts downloaded from aws
    s3_access_file = sys.argv[3] # must be the absolute path to the aws credentials file
    assert(os.path.isdir(working_dir))
    assert(os.path.isfile(s3_access_file))
    rm_working_dir = True
    if len(sys.argv) > 4:
        assert(sys.argv[4] == "True" or sys.argv[4] == "False")
        rm_working_dir = sys.argv[4] == "True"
    get_tissue_integration_results_csv(working_dir, s3_access_file, rm_working_dir)
