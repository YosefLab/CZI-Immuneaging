# This script can be used to generate various metrics that can be used to produce our data dashboard.
# - To get tissue coverage csv, run as follows:
#   python dashboard_utils.py tissue_coverage
# - To get tissue integration results csv, run as follows:
#   python dashboard_utils.py tissue_integration_results <working_dir> <s3_access_file>

import sys
import io
import csv
import numpy as np
import os
import traceback
import scanpy as sc

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

def get_tissue_integration_results_csv(working_dir: str, s3_access_file: str):
    logger.add_to_log("Downloading the Samples sheet from Google drive...")
    samples = utils.read_immune_aging_sheet("Samples")
    tissues = np.unique(samples["Organ"])

    # define the csv headers
    CSV_HEADER_TISSUE: str = "Tissue"
    CSV_HEADER_DONOR_COUNT: str = "# Donors"
    CSV_HEADER_DONORS: str = "Donors"
    CSV_HEADER_CELL_COUNT: str = "# Cells"
    CSV_HEADER_SAMPLE_COUNT: str = "# Samples"
    CSV_HEADER_SAMPLES: str = "Samples"
    CSV_HEADER_ANNDATA: str = "Anndata URL"

    # prepare the needful for downloading data from aws
    utils.set_access_keys(s3_access_file)
    BASE_AWS_URL = "https://immuneaging.s3.us-west-1.amazonaws.com"
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
        csv_row[CSV_HEADER_CELL_COUNT] = -1
        csv_row[CSV_HEADER_SAMPLE_COUNT] = -1
        csv_row[CSV_HEADER_SAMPLES] = []
        csv_row[CSV_HEADER_ANNDATA] = "Pending processing"

        # download the integrated anndata from aws
        try:
            file_name = "{}.{}.h5ad".format(tissue, version)

            # first check to see if we have an anndata for this tissue
            ls_cmd = "aws s3 ls {}/{}/{}/{}".format(BASE_S3_URL, BASE_S3_DIR, tissue, version)
            files = os.popen(ls_cmd).read()
            logger.add_to_log("aws response: {}\n".format(files))
            found = False
            for f in files.rstrip().split('\n'):
                if f.split('/')[-1] == file_name:
                    found = True
                    break

            # TODO remove once we have more or all tissues
            if not found:
                logger.add_to_log("Integrated annotated data not found for issue {}".format(tissue))
            else:
                csv_row[CSV_HEADER_ANNDATA] = "{}/{}/{}/{}/{}".format(BASE_AWS_URL, BASE_S3_DIR, tissue, version, file_name)
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
                csv_row[CSV_HEADER_DONORS] = np.unique(adata.obs["donor_id"])
                csv_row[CSV_HEADER_DONOR_COUNT] = len(csv_row[CSV_HEADER_DONORS])
                # clean up
                os.system("rm {}/*".format(working_dir))
        except Exception:
            logger.add_to_log("Execution failed for tissue {} with the following error:\n{}".format(tissue, traceback.format_exc()), "critical")
            os.system("rm {}/*".format(working_dir))
            logger.add_to_log("Continuing execution for other tissues...")
            continue
        csv_rows.append(csv_row)

    # write the csv
    csv_file = io.StringIO()
    field_names = [
        CSV_HEADER_TISSUE,
        CSV_HEADER_DONOR_COUNT,
        CSV_HEADER_DONORS,
        CSV_HEADER_CELL_COUNT,
        CSV_HEADER_ANNDATA,
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
    get_tissue_integration_results_csv(working_dir, s3_access_file)