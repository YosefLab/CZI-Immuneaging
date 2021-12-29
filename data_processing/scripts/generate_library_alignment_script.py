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
GEX_runs = []
BCR_runs = []
TCR_runs = []
for i in range(len(sample_indices)):
    gex = samples["GEX lib"][sample_indices[i]].split(',')
    assert len(gex)>0 # all samples must have a GEX lib (at least one)

    def get_feature_lib_if_exists(lib_name: str) -> List[str]:
        if pd.isna(samples[lib_name][sample_indices[i]]):
            lib = ['none' * len(gex)]
        else:
            lib = samples[lib_name][sample_indices[i]].split(',')
            assert len(lib) == len(gex)
        return lib

    adt = get_feature_lib_if_exists("ADT lib")
    hto = get_feature_lib_if_exists("HTO lib")
    for j in range(len(gex)):
        r = (gex[j], adt[j], hto[j])
        if r not in GEX_runs:
            GEX_runs.append(r)

    def get_bcr_tcr_libs(lib_type: str, lib_runs: List[str]) -> None:
        assert lib_type in ["BCR", "TCR"]
        if not pd.isna(samples["{} lib".format(lib_type)][sample_indices[i]]):
            libs = samples["{} lib".format(lib_type)][sample_indices[i]].split(',')
            assert len(libs) == len(gex)
            for l in libs:
                if l not in lib_runs:
                    lib_runs.append(l)

    get_bcr_tcr_libs("BCR", BCR_runs)
    get_bcr_tcr_libs("TCR", TCR_runs)

l_cd_mkdir = []
l_align_lib = []
l_align_msg = []
donor_run = "_".join([donor_id,seq_run])
for r in GEX_runs:
    # name the dir after the GEX lib's name
    # often the same lib (same lib id) is used for GEX, BCR and TCR thus distinguish the lib type in the dir name
    l_cd_mkdir.append(os.path.join(output_dir, "S3", donor_run, "{}-{}".format(r[0],"GEX"))) 
    l_align_lib.append("python {0}/align_library.py {1} {0} GEX {2}".format(code_path, configs_file, ",".join(r)))
    l_align_msg.append("echo \"Execution of align_library.py on GEX {} is complete.\"".format(",".join(r)))
for r in BCR_runs:
    l_cd_mkdir.append(os.path.join(output_dir, "S3", donor_run, "{}-{}".format(r,"BCR")))
    l_align_lib.append("python {0}/align_library.py {1} {0} BCR {2}".format(code_path, configs_file, r))
    l_align_msg.append("echo \"Execution of align_library.py on BCR {} is complete.\"".format(r))
for r in TCR_runs:
    l_cd_mkdir.append(os.path.join(output_dir, "S3", donor_run, "{}-{}".format(r,"TCR")))
    l_align_lib.append("python {0}/align_library.py {1} {0} TCR {2}".format(code_path, configs_file, r))
    l_align_msg.append("echo \"Execution of align_library.py on TCR {} is complete.\"".format(r))

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
