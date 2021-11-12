# This script can be used to generate various metrics that can be used to produce our data dashboard.
# - To get tissue coverage csv, run as follows:
#   python dashboard_utils.py tissue_coverage

import sys
import io
import csv
import numpy as np

import utils
from logger import RichLogger

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


action = sys.argv[1]
assert(action in ["tissue_coverage"])
get_tissue_coverage_csv()
