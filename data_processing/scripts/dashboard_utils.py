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
    donors = utils.read_immune_aging_sheet("Donors")
    logger.add_to_log("Downloading the Samples sheet from Google drive...")
    samples = utils.read_immune_aging_sheet("Samples")
    tissues = np.unique(samples["Organ"])

    csv_rows = []
    for donor in donors:
        csv_row = {tissue: "No Data" for tissue in tissues}
        for tissue in tissues:
            index = samples["Donor ID"] == donor
            organs = samples["Organ"][index] # organ and tissue are synonym
            stimulations = samples["Stimulation"][index]
            for i in range(len(organs)):
                tissue_stimulation = []
                # if for this donor, this tissue was found to have no data then set it to
                # stimulations[i]. If it was already reported to be stimulated with
                # stimulations[i] then no need to do anything. Otherwise, append
                # stimulations[i] to whatever was reported before
                if tissue_stimulation == []:
                    tissue_stimulation = stimulations[i]
                elif stimulations[i] not in tissue_stimulation:
                    tissue_stimulation.append(stimulations[i])
                csv_row[tissue] = tissue_stimulation.sort().join(" ➕ ")
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
