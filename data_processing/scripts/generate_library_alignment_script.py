"""
This script creates a bash script that can be used for executing the library alignment pipeline for a single donor.
The script gets a configuration file (for a specific donor), path to the immune aging code base, and a filename for the bash output file.

**NOTE**: the path to the code base and path to the configs file must be absolute paths.

Example:
python generate_library_alignment_script.py /path/to/code /path/to/align_libraries.configs.582C_001.txt 582C_001_align_libraries.sh
"""

import os
import sys
import pandas as pd

code_path = sys.argv[1]
configs_file = sys.argv[2]
alignment_script = sys.argv[3]

sys.path.append(code_path)
from utils import *

configs = load_configs(configs_file)
donor_id = configs["donor"]
output_dir = configs["output_destination"]
code_path = configs["code_path"]
python_env_version = configs["python_env_version"]
r_env_version = configs["r_env_version"]
seq_run = configs["seq_run"]
donor_run = "_".join([donor_id,seq_run])

samples = read_immune_aging_sheet(sheet="Samples")
sample_indices = samples.index[(samples["Donor ID"] == donor_id) & (samples["Seq run"] == float(seq_run))]
GEX_runs = []
BCR_runs = []
TCR_runs = []
for i in range(len(sample_indices)):
    gex = samples["GEX lib"][sample_indices[i]].split(',')
    assert len(gex)>0
    if pd.isna(samples["ADT lib"][sample_indices[i]]):
        adt = ['none' for j in range(len(gex))]
    else:
        adt = samples["ADT lib"][sample_indices[i]].split(',')
    if pd.isna(samples["HTO lib"][sample_indices[i]]):
        hto = ['none' for j in range(len(gex))]
    else:
        hto = samples["HTO lib"][sample_indices[i]].split(',')
    for j in range(len(gex)):
        r = (gex[j], adt[j], hto[j])
        if r not in GEX_runs:
            GEX_runs.append(r)
    if not pd.isna(samples["BCR lib"][sample_indices[i]]):
        bcr = samples["BCR lib"][sample_indices[i]].split(',')
        for j in bcr:
            if j not in BCR_runs:
                BCR_runs.append(j)
    if not pd.isna(samples["TCR lib"][sample_indices[i]]):
        tcr = samples["TCR lib"][sample_indices[i]].split(',')
        for j in tcr:
            if j not in TCR_runs:
                TCR_runs.append(j)

l_cd_mkdir = []
l_align_lib = []
l_align_msg = []
for r in GEX_runs:
    l_cd_mkdir.append(os.path.join(output_dir, "S3", donor_run, r[0]))
    l_align_lib.append("python {0}/align_library.py {1} {0} GEX {2}".format(code_path, configs_file, ",".join(r)))
    l_align_msg.append("echo \"Execution of align_library.py on GEX {} is complete.\"".format(",".join(r)))

# uncomment the following lines in order to allow aligning of BCR and TCR libraries
"""
for r in BCR_runs:
    l_cd_mkdir.append("mkdir -p " + os.path.join(output_dir, "S3", donor_id, r))
    l_align_lib.append("python {0}/align_library.py {1} {0} BCR {2}".format(code_path, configs_file, ",".join(r)))
for r in TCR_runs:
    l_cd_mkdir.append("mkdir -p " + os.path.join(output_dir, "S3", donor_id, r))
    l_align_lib.append("python {0}/align_library.py {1} {0} TCR {2}".format(code_path, configs_file, ",".join(r)))
"""

l1 = ["source activate {}".format(python_env_version),
    "mkdir -p " + os.path.join(output_dir, "S3"),
    "mkdir -p " + os.path.join(output_dir, "S3", donor_run)]
l2 = ["conda deactivate ",
    "source activate {}".format(r_env_version),
    #"Rscript ...", ## TODO complete the Rscript command
    "conda deactivate"] 

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
