# This script can be used to generate various metrics that can be used to produce our data dashboard.
# - To get tissue coverage csv, run as follows:
#   python dashboard_utils.py tissue_coverage

import sys
import io
import csv

import utils
from logger import RichLogger

logger = RichLogger()

def get_tissue_coverage_csv():
    # logger.add_to_log("Downloading the Donors sheet from from Google drive...")
    # donors = utils.read_immune_aging_sheet("Tissue Coverage (in progress)")

    logger.add_to_log("Downloading the Tissue Coverage sheet from Google drive...")
    tissues = utils.read_immune_aging_dashboard_sheet("Tissue Coverage (in progress)")

    logger.add_to_log("Downloading the Samples sheet from Google drive...")
    samples = utils.read_immune_aging_sheet("Samples")

    csv_rows = []
    for donor in donors:
        csv_row = {tissue: "No Data" for tissue in tissues}
        for tissue in tissues:
            index = samples["Donor ID"] == donor
            organs = samples["Organ"][index]
            stimulations = samples["Stimulation"][index]
            for i in range(len(organs)):
                stimulation_status = []
                # if for this donor, this tissue was found to have no data then set it to
                # stimulations[i]. If it was already reported to be stimulated with
                # stimulations[i] then no need to do anything. Otherwise, append it
                # to whatever was reported before
                if stimulation_status == []:
                    stimulation_status = stimulations[i]
                elif stimulations[i] not in stimulation_status:
                    stimulation_status.append(stimulations[i])
                csv_row[tissue] = stimulation_status.sort().join(" ➕ ")
                # if stimulations[i] == "anti-CD3/anti-CD28":
                #     # if for this donor, this tissue was found to have no data then set it to
                #     # "anti-CD3/anti-CD28". If it was already reported to be stimulated with
                #     # "anti-CD3/anti-CD28" then no need to do anything. Otherwise, append it
                #     # to whatever was reported before
                #     if csv_row[tissue] == "No Data":
                #         csv_row[tissue] = "anti-CD3/anti-CD28"
                #     elif not csv_row[tissue].contains("anti-CD3/anti-CD28"):
                #         csv_row[tissue] += "➕ " + "anti-CD3/anti-CD28"
                # elif stimulations[i] == "PMA+I":
                #     # same as above
                #     if csv_row[tissue] == "No Data":
                #         csv_row[tissue] = "PMA+I"
                #     elif not csv_row[tissue].contains("PMA+I"):
                #         csv_row[tissue] += "➕ " + "PMA+I"
                # elif stimulations[i] == "Nonstim":
                #     # same as above
                #     if csv_row[tissue] == "No Data":
                #         csv_row[tissue] = "PMA+I"
                #     elif not csv_row[tissue].contains("PMA+I"):
                #         csv_row[tissue] += "➕ " + "PMA+I"
        csv_rows.append(csv_row)

    # write the csv
    csv_file = io.StringIO()
    writer = csv.DictWriter(csv_file, fieldnames=tissues)
    writer.writerows(csv_rows)
    print(csv_file.getvalue())
    csv_file.close()


action = sys.argv[1]
assert(action in ["tissue_coverage"])
get_tissue_coverage_csv()