## This script generates configuration files for integrate_samples.py (one configuration file per tissue - for integrating all currently available samples from the given tissue type)
## run as follows (use absolute paths): python generate_tissue_integration_config_files.py <code_path> <output_destination> <s3_access_file>

import sys
import os
import json
import numpy as np
import re 

code_path = sys.argv[1]
output_destination = sys.argv[2]
s3_access_file = sys.argv[3]

sys.path.append(code_path)
from utils import *

python_env = "immune_aging.py_env.v4"

pilot_donors = ["CUIMC-457","CUIMC-471","390C"]

celltypist_model_urls = "https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_Low.pkl,https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_High.pkl,https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_Conde/v1/Immune_All_PIP.pkl"
leiden_resolutions = "1.0,2.0,3.0,5.0,7.0"

integration_configs = {
        "sandbox_mode": "False",
        "data_owner": "erahmani",
        "code_path": code_path,
        "output_destination": output_destination,
        "s3_access_file": s3_access_file,
        "integration_level": "tissue",
        "protein_levels_max_sds": 5,
        "n_highly_variable_genes": 3000,
        "highly_variable_genes_flavor": "seurat_v3",
        "batch_key": "donor_id",
        "scvi_max_epochs": 400,
        "totalvi_max_epochs": 400,
        "empirical_protein_background_prior": "False",
        # The following non-default configurations of scvi and totalvi can be used for speed up in case of very large numbers of cells.
        #"use_layer_norm": "none",
        #"use_batch_norm": "both",
        #"early_stopping": "True",
        #"early_stopping_patience": 45,
        #"batch_size": 1024,
        #"limit_train_batches": 20,
        "neighborhood_graph_n_neighbors": 15,
        "umap_min_dist": 0.5,
        "umap_spread": 1.0,
        "umap_n_components": 2,
        "celltypist_model_urls": celltypist_model_urls,
        "celltypist_dotplot_min_frac": 0.005,
        "leiden_resolutions": leiden_resolutions,
        "vdj_genes": "s3://immuneaging/vdj_genes/vdj_gene_list_v1.csv",
        "python_env_version": python_env,
        "r_setup_version": "immune_aging.R_setup.v2",
        "pipeline_version": "v3",
    }

set_access_keys(s3_access_file)

samples = read_immune_aging_sheet("Samples")
tissues = np.unique(samples["Organ"][np.logical_not(pd.isnull(samples["Organ"]))])

outfile = open(os.path.join(code_path,"integrate_samples_runs.sh"),'w') 
outfile.write("source activate {}\n".format(python_env))

no_integration = []

for tissue in tissues:
    # get all samples from the requested tissue as appears in the google spreadsheet
    final_sample_ids = []
    versions = []
    # consider all samples from the current tissue except for samples coming from pilot donor
    indices = (samples["Organ"] == tissue) & (~samples["Donor ID"].isin(pilot_donors))
    sample_ids = samples["Sample_ID"][indices]
    for sample_id in sample_ids:
        ls_cmd = "aws s3 ls s3://immuneaging/processed_samples/{}_GEX --recursive".format(sample_id)
        ls  = os.popen(ls_cmd).read()
        if len(ls) == 0:
            continue
        # find the latest version available
        filenames = ls.split("\n")
        latest_version = -1
        for filename in filenames:
            p = re.search("(.)(\d)+(.)log$", filename)
            if bool(p):
                version = int(filename[p.span()[0]+1:-4])
                if latest_version<version:
                    latest_version = version
        if latest_version>-1:
            version = "v"+str(latest_version)
            # check if h5ad file is available
            ls_cmd = "aws s3 ls s3://immuneaging/processed_samples/{0}_GEX/{1}/{0}_GEX.processed.{1}.h5ad".format(sample_id,version)
            if len(os.popen(ls_cmd).read())>0:
                versions.append(version)
                final_sample_ids.append(sample_id)
    if len(final_sample_ids)>1:
        tissue_integration_configs = integration_configs
        tissue_integration_configs["output_prefix"] = tissue
        tissue_integration_configs["sample_ids"] = ",".join(final_sample_ids)
        tissue_integration_configs["processed_sample_configs_version"] = ",".join([str(i) for i in versions])
        filename = os.path.join(output_destination,"integrate_samples.{}.configs.txt".format(tissue))
        with open(filename, 'w') as f:
            json.dump(tissue_integration_configs, f)
        print("generated configs file " + filename)
        outfile.write("python integrate_samples.py {}\n".format(filename))
    else:
        no_integration.append((tissue, len(final_sample_ids)))
        

outfile.write("conda deactivate")
outfile.close()

if len(no_integration) > 0:
    print("Integration config file was not generated for the following tissues that do not have more than one processed sample:")
    print("\n".join(["{}, number of processed samples:{}".format(no_integration[i][0],no_integration[i][1]) for i in range(len(no_integration))]))
