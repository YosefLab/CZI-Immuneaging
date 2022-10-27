## This script generates configuration files and bash scripts for running the integrate_using_scanvi.py script.
## There are currently two levels of integration: tissue-level and compartment-level.
## In tissue-level integration, this script generates one configuration file per tissue - for integrating all currently available samples from the given tissue type.
## In compartment-level integration, this script generates one file per compartment - for integrating cells belonging to that compartment from all currently available samples.
## run as follows (use absolute paths): python generate_scanvi_integration_config_files_and_scripts.py <code_path> <output_destination> <s3_access_file> <integration_level>

import sys
import os
import json

code_path = sys.argv[1]
output_destination = sys.argv[2]
s3_access_file = sys.argv[3]
integration_level = sys.argv[4]
assert integration_level in ["tissue", "compartment"]

sys.path.append(code_path)
from utils import *

python_env = "immune_aging.py_env.v4"
celltypist_model_urls = "https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_Low.pkl,https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_High.pkl"
leiden_resolutions = "1.0,2.0,3.0,5.0,7.0,10.0"
# these tissues are ignored since we only have a single non-pilot sample for them
skip_tissues = ["COL", "JEJ"]

integration_configs = {
        "sandbox_mode": "False",
        "data_owner": "erahmani",
        "code_path": code_path,
        "output_destination": output_destination,
        "s3_access_file": s3_access_file,
        "integration_level": integration_level,
        "batch_key": "donor_id" if integration_level == "tissue" else "donor_id,donor_id+tissue",
        "labels_key": "manual_labels",
        "unlabeled_category": "Unknown",
        "n_samples_per_label": 100,
        "n_latent": 30,
        "neighborhood_graph_n_neighbors": 15,
        "umap_min_dist": 0.5,
        "umap_spread": 1.0,
        "umap_n_components": 2,
        "celltypist_model_urls": celltypist_model_urls,
        "celltypist_dotplot_min_frac": 0.005,
        "leiden_resolutions": leiden_resolutions,
        "python_env_version": python_env,
        "r_setup_version": "immune_aging.R_setup.v2",
        "pipeline_version": "v4",
    }

# add the path to the latest annotated anndata object in s3. For now this is fairly simple but we can update it depending on
# what we decide on the nomenclature, organization and format of these files
s3_ann_path = f"s3://immuneaging/annotated_objects"
ann_version = get_latest_object_version(s3_access_file, s3_ann_path)
if ann_version == "v-1":
    print("Failed to find the latest annotated object. Terminating...")
    sys.exit()
integration_configs["latest_annotated_object_version"] = ann_version

outfile = open(os.path.join(code_path,"integrate_using_scanvi_runs.sh"),'w') 
outfile.write("source activate {}\n".format(python_env))

tissues_or_compartments = get_tissues_or_compartments(s3_access_file, integration_level, skip_tissues)
not_found = []
for tissue_or_compartment in tissues_or_compartments:
    s3_integrated_path = f"s3://immuneaging/integrated_samples/{integration_level}_level/{tissue_or_compartment}"
    version = get_latest_object_version(s3_access_file, s3_integrated_path)
    if version != "v-1":
        final_integration_configs = integration_configs
        final_integration_configs["latest_integrated_object_version"] = version
        final_integration_configs["output_prefix"] = tissue_or_compartment
        filename = os.path.join(output_destination,"integrate_using_scanvi.{}.configs.txt".format(tissue_or_compartment))
        with open(filename, 'w') as f:
            json.dump(final_integration_configs, f)
        print("generated configs file " + filename)
        outfile.write("python integrate_using_scanvi.py {}\n".format(filename))
    else:
        not_found.append(tissue_or_compartment)

outfile.write("conda deactivate")
outfile.close()

if len(not_found) > 0:
    print("Integration config file was not generated for the following tissues or compartments:")
    print("\n".join(not_found))
