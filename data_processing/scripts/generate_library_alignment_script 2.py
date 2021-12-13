"""
This script creates a bash script that can be used for executing the library alignment pipeline for a single donor.
The script gets a configuration file (for a specific donor), path to the immune aging code base, and a filename for the bash output file.

**NOTE**: the path to the code base and path to the configs file must be absolute paths.

Example:
python generate_library_alignment_script.py donor_id seq_run /path/to/code /path/to/align_libraries.configs.582C_001.txt python_env_version /path/to/output_dir
"""

import os
import sys
import pandas as pd

donor_id = sys.argv[1]
seq_run = sys.argv[2]
code_path = sys.argv[3]
configs_file = sys.argv[4]
python_env_version = sys.argv[5]
output_dir = sys.argv[6]

sys.path.append(code_path)
from utils import *

samples = read_immune_aging_sheet(sheet="Samples")
sample_indices = samples.index[(samples["Donor ID"] == donor_id) & (samples["Seq run"] == float(seq_run))]
runs = []
for i in range(len(sample_indices)):
    gex = samples["GEX lib"][sample_indices[i]].split(',')
    assert len(gex)>0 # all samples must have a GEX lib (at least one)

    def get_feature_lib_if_exists(lib_name: str) -> List[str]:
        if pd.isna(samples[lib_name][sample_indices[i]]):
            lib = ['none' for j in range(len(gex))]
        else:
            lib = samples[lib_name][sample_indices[i]].split(',')
            # TODO Q: This seems to hold empirically but experimentally what
            # guarantees that the number of adt and hto libs is the same as
            # the number of gex libs? Same for bcr/tcr libs
            assert len(lib) == len(gex)
        return lib
    adt = get_feature_lib_if_exists("ADT lib")
    hto = get_feature_lib_if_exists("HTO lib")
    bcr = get_feature_lib_if_exists("BCR lib")
    tcr = get_feature_lib_if_exists("TCR lib")
    for j in range(len(gex)):
        r = (gex[j], adt[j], hto[j], bcr[j], tcr[j])
        if r not in runs:
            runs.append(r)

l_cd_mkdir = []
l_align_lib = []
l_align_msg = []
donor_run = "_".join([donor_id,seq_run])
for r in runs:
    l_cd_mkdir.append(os.path.join(output_dir, "S3", donor_run, r[0])) # name the dir after the GEX lib's name
    l_align_lib.append("python {0}/align_library.py {1} {0} {2}".format(code_path, configs_file, ",".join(r)))
    l_align_msg.append("echo \"Execution of align_library.py on libs {} is complete.\"".format(",".join(r)))

l1 = ["source activate {}".format(python_env_version),
    "mkdir -p " + os.path.join(output_dir, "S3"),
    "mkdir -p " + os.path.join(output_dir, "S3", donor_run)]
l2 = ["conda deactivate"] 

alignment_script = "{}_align_libraries.sh".format(donor_run)
f = open(alignment_script,'w')
f.write("\n".join(l1) + "\n")
for i in range(len(l_cd_mkdir)):
    f.write("mkdir -p " + l_cd_mkdir[i] + "\n")
    f.write("cd " + l_cd_mkdir[i] + "\n")
    f.write(l_align_lib[i] + "\n")
    f.write(l_align_msg[i] + "\n")
f.write("\n".join(l2) + "\n")
f.close()
print("generated file " + alignment_script)
