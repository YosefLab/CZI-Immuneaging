## This script generates configuration files for integrate_samples.py.and
## There are currently two levels of integration: tissue-level and compartment-level.
## In tissue-level integration, this script generates one configuration file per tissue - for integrating all currently available samples from the given tissue type.
## In compartment-level integration, this script generates one file per compartment - for integrating cells belonging to that compartment from all currently available samples.
## run as follows (use absolute paths): python generate_integration_config_files_and_scripts.py <code_path> <output_destination> <s3_access_file> <integration_level>

import sys
import os
import json
import numpy as np
import re 

code_path = sys.argv[1]
output_destination = sys.argv[2]
s3_access_file = sys.argv[3]
integration_level = sys.argv[4]
apply_filtering = sys.argv[5]
assert integration_level in ["tissue", "compartment", "all"]

sys.path.append(code_path)
from utils import *

python_env = "immune_aging.py_env.v4"

pilot_donors = ["CUIMC-457","CUIMC-471","390C"]
bad_sample_id = ["647C-LIV-142", '759B-SKN-1', '768B-SKN-286', '768B-THY-287'] # 2xhigh ambient, 2xlow hashtag confidence due to high protein aggregates

celltypist_model_urls = "https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_Low.pkl,https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_High.pkl"
leiden_resolutions = "10.0" if integration_level=="all" else "3.0, 10.0"

integration_configs = {
        "folder_local_files": "processing_results_newest/process_integration/All",
        "compartment_barcode_csv_file": "Top_level_annotation_consensus.csv",
        "sandbox_mode": "False",
        "data_owner": "cane11",
        "code_path": code_path,
        "output_destination": output_destination,
        "s3_access_file": s3_access_file,
        "integration_level": integration_level,
        "protein_levels_max_sds": None,
        "n_highly_variable_genes": 10000,
        "highly_variable_genes_flavor": "seurat_v3",
        "batch_key": "donor_id",
        "empirical_protein_background_prior": "True",
        "n_layers": 2,
        "gene_likelihood": "nb",
        "scvi_max_epochs": 50 if integration_level!="tissue" else 200,
        "totalvi_max_epochs": 200,
        "early_stopping": True,
        "batch_size": 1024 if integration_level!="tissue" else 256,
        "reduce_lr_on_plateau": False,
        "n_epochs_kl_warmup": 10,
        "reduce_lr_on_plateau": False,
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
        "pipeline_version": "qc_230227_seq_batch",
        "include_stim": False,
        "filtering": {
            "apply_filtering": apply_filtering,
            "filter_name": "230310",
            "percolation_score_median":{
                "total_counts_median_cluster_scores": 0.5
            },
            "sum_percolation_score_mean_cluster": {
                "BLO":0.3,
                "BMA":0.3,
                "LLN":0.5,
                "JEJEPI":0.7,
                "JEJLP":0.3,
                "SPL":0.3,
                "LNG":0.3,
                "BAL":0.3,
                "ILN":0.3,
                "MLN":0.3,
                "COLEPI":0.3,
                "COLLP":0.3,
                "THY":0.3,
                "SKN":0.3,
                "LIV":0.3
                },
            "celltypes_passing_filtering": {
                "all": ['Mast cells',  'HSC/MPP'],
                "BMA": ['Mast cells',  'HSC/MPP', "Erythroid", 'Megakaryocytes/Platelets'],
                "SPL": ['Mast cells',  'Macrophages'],
                "LIV": ['Mast cells',  'Monocytes'],
                "MLN": ['Mast cells',  'Macrophages'],
                "LLN": ['Mast cells',  'Macrophages'],
                "ILN": ['Mast cells',  'Macrophages'],
            }
        }
    }

set_access_keys(s3_access_file)

samples = read_immune_aging_sheet("Samples")
tissues_or_compartments = []
if integration_level == "tissue":
    tissues_or_compartments = np.unique(samples["Organ"][np.logical_not(pd.isnull(samples["Organ"]))])
elif integration_level == "compartment":
    tissues_or_compartments = ["T", "B", "M"]
else:
    tissues_or_compartments = ["All"]
    
if integration_configs["filtering"]["apply_filtering"]=="True":
    filtered = "filtered_"
else:
    filtered = ""

outfile = open(os.path.join(output_destination,f"integrate_samples_runs_{filtered}{integration_level}.sh"),'w')

no_integration = []

done = False
for tissue_or_compartment in tissues_or_compartments:
    # get all samples from the requested tissue or compartment as appears in the google spreadsheet
    # consider all samples from the current tissue or compartment except for samples coming from pilot donor
    indices = ~samples["Donor ID"].isin(pilot_donors)
    if integration_level=="tissue":
        indices = indices & (samples["Organ"] == tissue_or_compartment)
    if integration_level == "tissue" or not done:
        final_sample_ids = []
        versions = []
        sample_ids = samples["Sample_ID"][indices]
        for sample_id in sample_ids:
            if sample_id in bad_sample_id:
                print(f'Removed {sample_id} from integration as it was defined in bad_sample_id.')
                continue
            ls_cmd = "aws s3 ls s3://immuneaging/processed_samples/{}_GEX --recursive | grep .log".format(sample_id)
            ls  = os.popen(ls_cmd).read()
            if len(ls) == 0:
                continue
            # find the latest version available
            filenames = ls.split("\n")[:-1]
            if len(filenames):
                for latest_version in sorted([int(i.split('.')[-2][1:]) for i in filenames], reverse=True):
                    version = "v"+str(latest_version)
                    # check if h5ad file is available
                    ls_cmd = "aws s3 ls s3://immuneaging/processed_samples/{0}_GEX/{1}/{0}_GEX.processed.{1}.h5ad".format(sample_id,version)
                    if len(os.popen(ls_cmd).read())>0:
                        versions.append(version)
                        final_sample_ids.append(sample_id)
                        break
        done = True
    final_integration_configs = integration_configs
    
    final_integration_configs["output_prefix"] = tissue_or_compartment if integration_configs['filtering']['apply_filtering']=='False' else tissue_or_compartment + '_filtered'
        
    if len(final_sample_ids):
        final_integration_configs["sample_ids"] = ",".join(final_sample_ids)
        final_integration_configs["processed_sample_configs_version"] = ",".join([str(i) for i in versions])
        if final_integration_configs["filtering"]["apply_filtering"]=="True":
            file_save = "integrate_samples.{}_filtered.configs.txt".format(tissue_or_compartment)
        else:
            file_save = "integrate_samples.{}.configs.txt".format(tissue_or_compartment)
        filename = os.path.join(output_destination, file_save)
        with open(filename, 'w') as f:
            json.dump(final_integration_configs, f)
        print("generated configs file " + filename)
        outfile.write("{}\n".format(filename))
    else:
        assert integration_level=="tissue" 
        no_integration.append((tissue_or_compartment, len(final_sample_ids)))

outfile.close()

if len(no_integration) > 0:
    print("Integration config file was not generated for the following tissues that do not have more than one processed sample:")
    print("\n".join(["{}, number of processed samples:{}".format(no_integration[i][0],no_integration[i][1]) for i in range(len(no_integration))]))
