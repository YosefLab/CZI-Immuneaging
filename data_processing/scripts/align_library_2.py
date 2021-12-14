import sys
import os
import time
import scanpy as sc
import pandas as pd
import numpy as np
import re
import hashlib

align_lib_script = sys.argv[0]
configs_file = sys.argv[1]
code_path = sys.argv[2]
lib_ids = sys.argv[3]

# we don't need this, but we use it to name the output h5ad file, which is
# used downstream and expects the "GEX" in the name, so keep it around for
# now
lib_type = "GEX"

sys.path.append(code_path)

from logger import SimpleLogger
from utils import *

VARIABLE_CONFIG_KEYS = ["donor",
"seq_run",
"output_destination",
"aligner_software_path",
"alignment_ref_genome_path",
"berkeley_user",
"s3_access_file",
"code_path"]

def get_aligner_cmd(aligner, donor_id, seq_run, data_dir, data_dir_fastq, samples, cite_key, chemistry, GEX_lib = None, ADT_lib = None, HTO_lib = None, TCR_lib = None, BCR_lib = None, protein_panel = None):
    assert aligner == "cellranger" # no other option is currently implemented
    assert GEX_lib is not None
    # set chemistry
    assert chemistry == "5'v1" or chemistry == "5'v2" or chemistry == "3'v2" or chemistry == "3'v3"
    if chemistry == "5'v1" or chemistry == "5'v2":
        chem = "fiveprime"
    if chemistry == "3'v2":
        chem = "SC3Pv2"
    if chemistry == "3'v3":
        chem = "SC3Pv3"
    libs = []
    GEX_lib_name = "_".join([donor_id, seq_run, "GEX", GEX_lib])
    outputs_to_save = [os.path.join(data_dir,GEX_lib_name,"outs/web_summary.html"),
        os.path.join(data_dir,GEX_lib_name,"outs/metrics_summary.csv"),
        os.path.join(data_dir,GEX_lib_name,"outs/cloupe.cloupe")]
    aligned_data_dir = os.path.join(data_dir, GEX_lib_name,"outs/filtered_feature_bc_matrix/")
    if (ADT_lib is not None) or (HTO_lib is not None):
        # prepare libraries csv file - see https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis#libraries-csv
        libs.append([os.path.join(data_dir_fastq, GEX_lib), GEX_lib_name, "Gene Expression"])
        if ADT_lib is not None:
            ADT_lib_name = "_".join([donor_id, seq_run, "ADT", ADT_lib])
            libs.append([os.path.join(data_dir_fastq, ADT_lib), ADT_lib_name, "Antibody Capture"])
        if HTO_lib is not None:
            HTO_lib_name = "_".join([donor_id, seq_run, "HTO", HTO_lib])
            libs.append([os.path.join(data_dir_fastq, HTO_lib), HTO_lib_name, "Antibody Capture"])
        libraries_df = pd.DataFrame(libs, columns = ['fastqs', 'sample', 'library_type'])
        libraries_csv_file = os.path.join(data_dir,"libraries.csv")
        libraries_df.to_csv(libraries_csv_file, sep=',', header=True, index=False)	
        outputs_to_save.append(libraries_csv_file)
        # prepare feature reference csv file from HTO and/or ADT.
        # first, prep feature information for the hashtags
        sample_indices = []
        GEX_libs = samples["GEX lib"].values
        for i in range(len(GEX_libs)):
            if not pd.isna(GEX_libs[i]) and GEX_lib in GEX_libs[i]:
                sample_indices.append(i)
        feature_ref = []
        hashtags = samples["Hashtag"][sample_indices].values
        hashtag_ids = samples["Sample_ID"][sample_indices].values
        counter = 0
        for hashtag in hashtags:
            for i in range(cite_key.shape[0]):
                if hashtag in cite_key.loc[i].values:
                    feature_ref.append([hashtag_ids[counter], hashtag_ids[counter], cite_key.loc[i].values[6], cite_key.loc[i].values[5], cite_key.loc[i].values[3], "Antibody Capture"])
                    counter += 1
                    continue
        feature_ref_df = pd.DataFrame(feature_ref, columns = ['id', 'name', 'read', 'pattern', 'sequence', 'feature_type'])
        # second, incorporate protein panel if applicable
        if ADT_lib is not None:
            assert protein_panel is not None
            feature_ref_df = pd.concat([feature_ref_df, protein_panel])
        feature_ref_csv_file = os.path.join(data_dir,"feature_ref.csv")
        feature_ref_df.to_csv(feature_ref_csv_file, sep=',', header=True, index=False)
        outputs_to_save.append(feature_ref_csv_file)
        if BCR_lib is not None or TCR_lib is not None:
            # prepare config csv file for cellranger multi - see https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/multi#config
            multi_config_csv_file = os.path.join(data_dir,"multi_config.csv")
            outputs_to_save.append(multi_config_csv_file)
            aligner_cmd = "{0} multi --id={1} --csv={2}".format(aligner_software_path, GEX_lib_name, multi_config_csv_file)
        else:
            aligner_cmd = "{0} count --id={1} --libraries={2} --transcriptome={3} --feature-ref={4} --chemistry={5}".format(aligner_software_path, GEX_lib_name, libraries_csv_file, aligner_genome_file, feature_ref_csv_file, chem)
    #elif BCR_lib is not None or TCR_lib is not None:
        # TODO
        #pass
    else:
        # GEX with no hashtags, ADT, BCR or TCR
        aligner_cmd = "{0} count --id={1} --transcriptome={2} --fastqs={3} --chemistry={4}".format(aligner_software_path, 
            GEX_lib_name, aligner_genome_file, os.path.join(data_dir_fastq, GEX_lib), chem)
    return (aligner_cmd, aligned_data_dir, outputs_to_save)

configs_dir_remote = "s3://immuneaging/aligned_libraries/configs/"
configs_file_remote_prefix = "align_library."

lib_ids = lib_ids.split(',')
configs = load_configs(configs_file)
set_access_keys(configs["s3_access_file"])
donor_id = configs["donor"]
seq_run = configs["seq_run"]
donor_run = "_".join([donor_id,seq_run])
aligner = configs["aligner"]
aligner_software_path = configs["aligner_software_path"]
aligner_genome_file = os.path.join(configs["alignment_ref_genome_path"],configs["alignment_ref_genome_file"])

data_dir = os.path.join(configs["output_destination"],"S3",donor_run,lib_ids[0])

configs_version = get_configs_version_alignment(configs, data_dir, configs_dir_remote, configs_file_remote_prefix, VARIABLE_CONFIG_KEYS)

# TODO added for debugging, delete
print("*** DEBUG: configs_version: {}".format(configs_version))
if configs_version != "v2":
    sys.exit()

h5ad_file = "{}.{}.{}.h5ad".format(donor_run, lib_ids[0], configs_version)
logger_file = os.path.join("align_library.{}.{}.{}.log".format(donor_run,lib_ids[0],configs_version))
h5ad_file_exists = False
logger_file_exists = False

# check if aligned files are already on the server
ls_cmd = 'aws s3 ls s3://immuneaging/aligned_libraries --recursive'
aligned_files = os.popen(ls_cmd).read()
for f in aligned_files.rstrip().split('\n'):
    j = f.split('/')[-1]
    if j == h5ad_file:
        h5ad_file_exists = True
    if j == logger_file:
        logger_file_exists = True

if logger_file_exists and h5ad_file_exists:
    print("Files {0} and {1} already exist on the S3 bucket; exiting align_library.py".format(h5ad_file, logger_file))
    sys.exit()

logger_file_path = os.path.join(data_dir, logger_file)
if os.path.isfile(logger_file_path):
    os.remove(logger_file_path)
logger = SimpleLogger(filename = logger_file_path)

logger.add_to_log("Running align_library.py...")
logger.add_to_log("Starting time: {}".format(get_current_time()))
with open(align_lib_script, "r") as f:
    logger.add_to_log("align_library.py md5 checksum: {}".format(hashlib.md5(bytes(f.read(), 'utf-8')).hexdigest()))

logger.add_to_log("donor: {}\nSeq run: {}\nlib_ids: {}".format(donor_id, seq_run, lib_ids))
logger.add_to_log("using the following configurations:\n{}".format(str(configs)))

logger.add_to_log("Loading metadata from the Google spreadsheet...")
samples = read_immune_aging_sheet("Samples")
donors = read_immune_aging_sheet("Donors")
cite_key = read_immune_aging_sheet("CITE key")

logger.add_to_log("Extracting site information...")
sites = samples["Site"][(samples["Donor ID"] == donor_id) & (samples["Seq run"] == float(seq_run))].values
assert all(sites == "UK") or all(sites == "NY")
site = sites[0]
site_s3_dir = "s3://immuneaging/raw_"
if site == "UK":
    site_s3_dir = site_s3_dir+"sanger/"
else:
    site_s3_dir = site_s3_dir+"columbia/"

logger.add_to_log("detected site_s3_dir = {}".format(site_s3_dir))

logger.add_to_log("Extracting chemistry information...")
fields = ["GEX lib","GEX chem"]
if lib_ids[1] != "none":
    fields.append("CITE chem")

if lib_ids[2] != "none":
    fields.append("HTO chem")

if lib_ids[3] != "none":
    fields.append("BCR chem")

if lib_ids[4] != "none":
    fields.append("TCR chem")

libs_chem = samples[fields]
lib_chems = []
for i in range(libs_chem.shape[0]):
    if not pd.isna(libs_chem.iloc[i,0]) and lib_ids[0] in libs_chem.iloc[i,0]:
        lib_chems.append(libs_chem.iloc[i,1:].values)

# all libraries must have the same chemistry
lib_chems_flat = np.array(lib_chems).flatten()
assert(all(lib_chems_flat == lib_chems_flat[0]))
chemistry = lib_chems_flat[0]
logger.add_to_log("detected chemistry = {}".format(chemistry))

logger.add_to_log("Downloading fastq files from S3...")
# don't rely on 'aws s3 sync' commands which seem to be buggy in some cases; copy file if it doesn't exist on local.
ls_cmd = 'aws s3 ls {0}'.format(site_s3_dir)
ls = os.popen(ls_cmd).read()
data_dir_fastq = os.path.join(data_dir, "fastq")
os.system("mkdir -p " + data_dir_fastq)
for lib_id in lib_ids:
    if lib_id != "none":
        logger.add_to_log("Downloading fastq files from S3 for lib {}...".format(lib_id))
        data_dir_lib = os.path.join(data_dir_fastq, lib_id)
        os.system("mkdir -p " + data_dir_lib)
        lib_pattern = "{}_{}.*{}.*.fastq.gz".format(donor_id, seq_run, lib_id)
        for i in re.split(' |\n',ls.rstrip()):
            if re.search(lib_pattern, i):
                f = os.path.join(data_dir_lib,i)
                if os.path.isfile(f):
                    logger.add_to_log("file {} is already in local.".format(f))
                else:
                    cp_cmd = 'aws s3 cp --no-progress {0} {1}'.format(os.path.join(site_s3_dir,i), data_dir_lib)
                    logger.add_to_log("cp_cmd: {}".format(cp_cmd))
                    logger.add_to_log("aws response: {}".format(os.popen(cp_cmd).read()))

logger.add_to_log("Preparing alignment command...")
GEX_lib = lib_ids[0]
ADT_lib = None if lib_ids[1] == "none" else lib_ids[1]
HTO_lib = None if lib_ids[2] == "none" else lib_ids[2]
BCR_lib = None if lib_ids[3] == "none" else lib_ids[3]
TCR_lib = None if lib_ids[4] == "none" else lib_ids[4]
if lib_ids[1] == "none":
    protein_panel = None
else:
    # get protein panel
    GEX_libs = samples["GEX lib"].values
    for i in range(len(GEX_libs)):
        if not pd.isna(GEX_libs[i]) and GEX_lib in GEX_libs[i]:
            protein_panel = read_immune_aging_sheet(samples["Protein panel"].loc[i])
            protein_panel = protein_panel[["id", "name", "read", "pattern", "sequence", "feature_type"]]
            break
alignment_cmd, aligned_data_dir, aligner_outputs_to_save = get_aligner_cmd(aligner, donor_id, seq_run, data_dir, data_dir_fastq, samples, cite_key, chemistry, GEX_lib, ADT_lib, HTO_lib, TCR_lib, BCR_lib, protein_panel)

alignment_exists = dir_and_files_exist(aligned_data_dir, aligner_outputs_to_save)
prefix = "_".join([donor_id, seq_run, lib_type, GEX_lib])
if alignment_exists:
    logger.add_to_log("Alignment outputs for the following command already exist:")
    logger.add_to_log("alignment_cmd:\n{}".format(alignment_cmd))
    logger.add_to_log("Skipping alignment.")
else:
    logger.add_to_log("Running the following alignment command:\n{}".format(alignment_cmd))
    alignment_output = os.popen(alignment_cmd).read()
    alignment_exists = dir_and_files_exist(aligned_data_dir, aligner_outputs_to_save)
    if not alignment_exists:
        # remove the output directory, which is required in order to prevent errors in a following execution of cellranger		
        os.system("rm -r {}".format(os.path.join(data_dir, prefix)))
        if "We detected an unsupported chemistry combination (SC5P-R2, SC5P-PE)" in alignment_output:
            logger.add_to_log("Alignment failed due to: an unsupported chemistry combination (SC5P-R2, SC5P-PE).", level="error")
            logger.add_to_log("Rerunning after changing chemistry argument...")
            alignment_cmd = alignment_cmd[0:alignment_cmd.index("--chemistry=")] + "--chemistry=SC5P-R2"
            logger.add_to_log("alignment_cmd:\n{}".format(alignment_cmd))
            logger.add_to_log("Output from aligner:\n" + os.popen(alignment_cmd).read())
            alignment_exists = dir_and_files_exist(aligned_data_dir, aligner_outputs_to_save)
        else:
            logger.add_to_log("Alignment failed. Alignment output:\n{}".format(alignment_output), level="error")

if not alignment_exists:
    msg = "Not all alignment outputs were generated. Terminating execution."
    logger.add_to_log(msg, level="error")
    print(msg)
    sys.exit()

logger.add_to_log("Uploading aligner outputs to S3...")
for out in aligner_outputs_to_save:
    l = out.split('/')
    out_file = "{0}.{1}.{2}".format(prefix,aligner,l[-1])
    if len(l) == 1:
        out_dir = "."
    else:
        out_dir = "/".join(l[0:-1])
    cp_cmd = "cp {0} {1}".format(out, os.path.join(data_dir, out_file))
    logger.add_to_log("copying file in local...")
    logger.add_to_log("cp_cmd: {}".format(cp_cmd))
    logger.add_to_log("cp_cmd result: {}".format(os.popen(cp_cmd).read()))
    sync_cmd = 'aws s3 sync --no-progress {0} s3://immuneaging/aligned_libraries/{1}/{2} --exclude "*" --include {3}'.format(data_dir, configs_version, prefix, out_file)
    logger.add_to_log("Uploading aligner output {}...".format(out_file))
    logger.add_to_log("sync_cmd: {}".format(sync_cmd))
    logger.add_to_log("aws response: {}".format(os.popen(sync_cmd).read()))

logger.add_to_log("Converting aligned data to h5ad...")
adata = sc.read_10x_mtx(aligned_data_dir, gex_only = False)
logger.add_to_log("Saving file as {}...".format(os.path.join(data_dir, h5ad_file)))
adata.write(os.path.join(data_dir, h5ad_file), compression="lzf")

logger.add_to_log("Uploading h5ad file to S3...")
sync_cmd = 'aws s3 sync --no-progress {0} s3://immuneaging/aligned_libraries/{1}/{2} --exclude "*" --include {3}'.format(data_dir, configs_version, prefix, h5ad_file)
logger.add_to_log("sync_cmd: {}".format(sync_cmd))
logger.add_to_log("aws response: {}".format(os.popen(sync_cmd).read()))

msg = "Done aligning library {}".format(GEX_lib)
logger.add_to_log(msg)
print(msg)

# upload log file to S3
cmd = 'aws s3 sync --no-progress {0} s3://immuneaging/aligned_libraries/{1}/{2} --exclude "*" --include {3}'.format(data_dir, configs_version, prefix, logger_file.split('/')[-1])
os.system(cmd)

## remove fastq files
# os.system("rm -r {}".format(data_dir_fastq))