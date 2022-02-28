## Run as follows: python gather_lib_alignment_qcs.py <code_path> <output_destination> <s3_access_file> <task_type>

import re
import sys
import os
from typing import List
import pandas as pd
import seaborn as sns
import scanpy as sc
import matplotlib.pyplot as plt
import logging

from logger import RichLogger

logging.getLogger('matplotlib').setLevel(logging.WARNING)

code_path = sys.argv[1]
output_destination = sys.argv[2]
s3_access_file = sys.argv[3]
task_type = sys.argv[4]

assert task_type in ["csv", "adata"]

sys.path.append(code_path)
from utils import *

set_access_keys(s3_access_file)

logger = RichLogger()

CSV_FIELDS_FOR_GEX = {
    "Estimated Number of Cells",
    "Mean Reads per Cell",
    "Median Genes per Cell",
    "Number of Reads",
    "Sequencing Saturation",
    "Median UMI Counts per Cell",
    # the below are only present if we have CITE or HTO data
    "Antibody: Number of Reads",
    "Antibody: Mean Reads per Cell",
    "Antibody: Sequencing Saturation",
    "Antibody: Median UMIs per Cell (summed over all recognized antibody barcodes)",
}

CSV_FIELDS_FOR_IR = {
    "Estimated Number of Cells",
    "Mean Read Pairs per Cell",
    "Number of Read Pairs",
    "Reads Mapped to Any V(D)J Gene",
    # BCR
    "Median IGH UMIs per Cell",
    "Median IGK UMIs per Cell",
    "Median IGL UMIs per Cell",
    # TCR
    "Median TRA UMIs per Cell",
    "Median TRB UMIs per Cell",
}

def combine_data(samples, donor_id, seq_run, site):
    indices = (samples["Donor ID"] == donor_id) & (samples["Seq run"] == float(seq_run))

    def add_lib(lib_type: str, all_libs: set) -> None:
        if lib_type not in ["GEX", "BCR", "TCR"]:
            raise ValueError("Unsupported lib_type: {}. Must be one of: GEX, BCR, TCR".format(lib_type))
        column_name = "{} lib".format(lib_type)
        libs_all = samples[indices][column_name]
        failed_libs = set()
        ls_cmd = "aws s3 ls s3://immuneaging/aligned_libraries --recursive"
        ls  = os.popen(ls_cmd).read()
        for i in range(len(libs_all)):
            if libs_all.iloc[i] is np.nan:
                continue
            libs = libs_all.iloc[i].split(",")
            for lib in libs:
                # find the latest aligned_lib_version
                latest_version = -1
                if len(ls) != 0:
                    filenames = ls.split("\n")
                    if lib_type == "GEX":
                        # search for patterns of <lib id>.vX.h5ad. If there is a match, group
                        # one is "<lib id>.v" and group 2 is "X" (X can be any integer >0)
                        pattern = "({}\.v)(\d+)\.h5ad$".format(lib)
                    else:
                        pattern = "({}_{}\.cellranger\.filtered_contig_annotations\.v)(\d+)\.csv".format(lib_type, lib)
                    for filename in filenames:
                        m = re.search(pattern, filename)
                        if bool(m):
                            version = int(m[2])
                            if latest_version < version:
                                latest_version = version
                if latest_version == -1:
                    failed_libs.add(lib)
                    # skip
                    continue
                aligned_lib_version = "v" + str(latest_version)
                all_libs.add((lib,lib_type,aligned_lib_version))
        for fl in failed_libs:
            logger.add_to_log("No aligned libraries found on AWS for lib id {} lib type {}. Skipping.".format(fl, lib_type), level="warning")

    def combine_metrics_for_lib(libs: set, generic_lib_type: str, data_dir: str):
        lib_data_dir = os.path.join(data_dir, generic_lib_type)
        os.system("mkdir -p " + lib_data_dir)
        all_metrics_files = []
        for lib in libs:
            lib_id = lib[0]
            lib_type = lib[1]
            aligned_lib_version = lib[2]
            logger.add_to_log("Downloading metrics.csv file for lib id {}, lib type {} from S3...".format(lib_id, lib_type))
            metrics_csv_file_name = "{}_{}_{}_{}.cellranger.metrics_summary.csv".format(donor_id, seq_run, lib_type, lib_id)
            sync_cmd = 'aws s3 sync --no-progress s3://immuneaging/aligned_libraries/{}/{}_{}_{}_{}/ {} --exclude "*" --include {}'.format(
                aligned_lib_version, donor_id, seq_run, lib_type, lib_id, lib_data_dir, metrics_csv_file_name
            )
            logger.add_to_log("sync_cmd: {}".format(sync_cmd))
            logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
            metrics_csv_file = os.path.join(lib_data_dir, metrics_csv_file_name)
            if not os.path.isfile(metrics_csv_file):
                msg = "Failed to download file {} from S3.".format(metrics_csv_file)
                logger.add_to_log(msg, level="error")
                raise ValueError(msg)
            all_metrics_files.append((metrics_csv_file,lib_id,lib_type,aligned_lib_version))

        # create the combined csv metrics file
        all_dfs = []
        for f in all_metrics_files:
            df = pd.read_csv(f[0])
            fields_to_keep = CSV_FIELDS_FOR_GEX if generic_lib_type == "GEX" else CSV_FIELDS_FOR_IR
            for c in df.columns:
                if c not in fields_to_keep:
                    del df[c]
            df["Donor ID"] = donor_id
            df["Site (UK/NY)"] = site
            df["Lib ID"] = f[1]
            df["Lib Type"] = f[2]
            df["Aligned Lib Version"] = f[3]
            all_dfs.append(df)
        if len(all_dfs) > 0:
            combined_df = pd.concat(all_dfs, ignore_index=True)
            combined_metrics = os.path.join(data_dir, "{}_{}_all_{}_metrics.csv".format(donor_id,seq_run,generic_lib_type))
            with open(combined_metrics, 'w') as f:
                combined_df.to_csv(f)
            # upload the combined csv file to AWS
            logger.add_to_log("Uploading combined metrics file {} to S3...".format(combined_metrics.split("/")[-1]))
            sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/combined_lib_alignment_metrics/{}/ --exclude "*" --include {}'.format(data_dir, generic_lib_type, combined_metrics.split("/")[-1])
            logger.add_to_log("sync_cmd: {}".format(sync_cmd))
            logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
        else:
            combined_metrics = ""

        return combined_metrics

    def combine_adatas_for_lib(libs: set, generic_lib_type: str, data_dir: str):
        if generic_lib_type != "GEX":
            raise NotImplementedError
        lib_data_dir = os.path.join(data_dir, generic_lib_type)
        os.system("mkdir -p " + lib_data_dir)
        all_adata_dfs = []
        for lib in libs:
            lib_id = lib[0]
            lib_type = lib[1]
            aligned_lib_version = lib[2]
            logger.add_to_log("Downloading aligned h5ad file for lib id {}, lib type {} from S3...".format(lib_id, lib_type))
            aligned_h5ad_file_name = "{}_{}.{}.{}.h5ad".format(donor_id, seq_run, lib_id, aligned_lib_version)
            sync_cmd = 'aws s3 sync --no-progress s3://immuneaging/aligned_libraries/{}/{}_{}_{}_{}/ {} --exclude "*" --include {}'.format(
                aligned_lib_version, donor_id, seq_run, lib_type, lib_id, lib_data_dir, aligned_h5ad_file_name
            )
            logger.add_to_log("sync_cmd: {}".format(sync_cmd))
            logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
            aligned_h5ad_file = os.path.join(lib_data_dir, aligned_h5ad_file_name)
            if not os.path.isfile(aligned_h5ad_file):
                msg = "Failed to download file {} from S3.".format(aligned_h5ad_file)
                logger.add_to_log(msg, level="error")
                raise ValueError(msg)
            logger.add_to_log("Reading aligned h5ad file...")
            adata = sc.read_h5ad(aligned_h5ad_file)
            a = adata.X.toarray()
            umi_counts = a.sum(axis=-1)
            gene_counts = np.count_nonzero(a, axis=-1)
            df = pd.DataFrame(data={'Gene Counts': gene_counts, "UMI Counts": umi_counts}, index=adata.obs.index)
            all_adata_dfs.append((df,lib_id,lib_type,aligned_lib_version))

        return all_adata_dfs

    gex_libs = set()
    ir_libs = set()
    add_lib("GEX", gex_libs)
    add_lib("BCR", ir_libs)
    add_lib("TCR", ir_libs)

    # create a new directory for the outputs
    data_dir = os.path.join(output_destination, "_".join([donor_id, seq_run]))
    os.system("mkdir -p " + data_dir)

    if task_type == "csv":
        combined_gex = combine_metrics_for_lib(gex_libs, "GEX", data_dir)
        combined_ir = combine_metrics_for_lib(ir_libs, "IR", data_dir)
    else:
        combined_gex = combine_adatas_for_lib(gex_libs, "GEX", data_dir)
        # no IR for now
        combined_ir = []

    return combined_gex, combined_ir

def combine_csv_all_donors(lib_type: str, per_donor_files: List[str]):
    per_donor_csv = [pd.read_csv(f) for f in per_donor_files if os.path.isfile(f)]
    combined_df = pd.concat(per_donor_csv, ignore_index=True)
    all_donors_metrics = os.path.join(output_destination, "all_donors_{}_metrics.csv".format(lib_type))
    with open(all_donors_metrics, 'w') as f:
        combined_df.to_csv(f)
    logger.add_to_log("☑ Uploading combined metrics across all donors for lib type {} to S3...".format(lib_type))
    sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/combined_lib_alignment_metrics/{}/ --exclude "*" --include {}'.format(output_destination, lib_type, all_donors_metrics.split("/")[-1])
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

def plot_data_all_donors(lib_type: str, per_donor_data: List[pd.DataFrame]):
    if lib_type != "GEX":
        raise NotImplementedError
    concat_dfs = []
    for l in per_donor_data :
        # l (one per donor) is a list of Dataframes, each of which has data for a given lib id that we need to plot
        for lib in l:
            df = lib[0]
            lib_id = lib[1]
            # cut out outliers
            q = df["Gene Counts"].quantile(0.99)
            df = df[df["Gene Counts"] < q]
            q = df["UMI Counts"].quantile(0.99)
            df = df[df["UMI Counts"] < q]
            # melt
            df = df.melt(var_name="Count type", value_name="Counts")
            # add lib id
            df["Lib id"] = lib_id
            # add it to the list
            concat_dfs.append(df)
            # sns.boxplot(x="Count type", y="Counts", data=df).set_title(lib_id)
    d = pd.concat(concat_dfs)
    d_file = os.path.join(output_destination, "all_donors_per_{}_lib_counts_data.csv".format(lib_type))
    with open(d_file, 'w') as f:
        d.to_csv(f)
    logger.add_to_log("☑ Uploading combined lib data across all donors for lib type {} to S3...".format(lib_type))
    sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/combined_lib_alignment_metrics/{}/ --exclude "*" --include {}'.format(output_destination, lib_type, d_file.split("/")[-1])
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))
    # now plot
    g = sns.catplot(x="Count type", y="Counts", col="Lib id", data=d, kind="violin", col_wrap=4).set(xlabel=None, ylabel=None)
    plt.show()
    fig_path = os.path.join(output_destination, "all_donors_per_{}_lib_counts.svg".format(lib_type))
    g.fig.savefig(fig_path, dpi=100)
    logger.add_to_log("☑ Uploading combined lib plots across all donors for lib type {} to S3...".format(lib_type))
    sync_cmd = 'aws s3 sync --no-progress {} s3://immuneaging/combined_lib_alignment_metrics/{}/ --exclude "*" --include {}'.format(output_destination, lib_type, fig_path.split("/")[-1])
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}\n".format(os.popen(sync_cmd).read()))

donors = read_immune_aging_sheet("Donors")
samples = read_immune_aging_sheet("Samples")
per_donor_gex_combined = []
per_donor_ir_combined = []
for donor in np.unique(donors["Donor ID"]):
    indices = donors["Donor ID"] == donor
    sites = np.unique(donors[indices]["Site (UK/ NY)"])
    for site in sites:
        indices = samples["Donor ID"] == donor
        seq_runs = np.unique(samples[indices]["Seq run"])
        for seq_run in seq_runs:
            seq_run = "00" + str(seq_run)
            combined_gex, combined_ir = combine_data(samples, donor, seq_run, site)
            per_donor_gex_combined.append(combined_gex)
            per_donor_ir_combined.append(combined_ir)
if task_type == "csv":
    combine_csv_all_donors("GEX", per_donor_gex_combined)
    combine_csv_all_donors("IR", per_donor_ir_combined)
else:
    plot_data_all_donors("GEX", per_donor_gex_combined)
