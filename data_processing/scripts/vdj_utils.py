import io
import os
import pandas as pd
import numpy as np
import anndata
from anndata import AnnData
import csv
from typing import List
import gc
from utils import *
import uuid
from logger import SimpleLogger

def get_vdj_lib_to_gex_lib_mapping(samples=None):
    # Returns a mapping of all vdj libraries to their corresponding gex libraries
    # in the form of two dictionaries, one for bcr libs and one for tcr libs
    if samples is None:
        samples = read_immune_aging_sheet("Samples")

    def add_lib(lib_type: str, all_libs: dict) -> None:
        if lib_type not in ["BCR", "TCR"]:
            raise ValueError("Unsupported lib_type: {}. Must be one of: BCR, TCR".format(lib_type))
        column_name = "{} lib".format(lib_type)
        libs_all = samples[column_name]
        gex_libs_all = samples["GEX lib"]
        for i in range(len(libs_all)):
            if libs_all.iloc[i] is np.nan:
                continue
            libs = libs_all.iloc[i].split(",")
            gex_libs = gex_libs_all.iloc[i].split(",")
            for j in range(len(libs)):
                lib = libs[j]
                corresponding_gex_lib = gex_libs[j]
                all_libs[lib] = corresponding_gex_lib

    all_bcr_libs = {}
    add_lib("BCR", all_bcr_libs)

    all_tcr_libs = {}
    add_lib("TCR", all_tcr_libs)

    return all_bcr_libs, all_tcr_libs

def get_gex_lib_to_vdj_lib_mapping():
    bcr_to_gex, tcr_to_gex = get_vdj_lib_to_gex_lib_mapping()

    # make sure there is a 1:1 mapping between gex and bcr, and same for tcr
    gex = list(bcr_to_gex.values())
    if len(gex) != len(np.unique(gex)):
        raise ValueError("There is not a 1:1 mapping between gex libs and bcr libs")
    gex = list(tcr_to_gex.values())
    if len(gex) != len(np.unique(gex)):
        raise ValueError("There is not a 1:1 mapping between gex libs and tcr libs")

    gex_to_bcr = {v:k for k,v in bcr_to_gex.items()}
    gex_to_tcr = {v:k for k,v in tcr_to_gex.items()}

    return gex_to_bcr, gex_to_tcr

def add_vdj_lib_ids_to_adata(adata: AnnData, gex_to_bcr, gex_to_tcr):
    adata.obs["bcr_library_id"] = np.array([gex_to_bcr[l] if l in gex_to_bcr else "NA" for l in adata.obs["library_id"].values])
    adata.obs["tcr_library_id"] = np.array([gex_to_tcr[l] if l in gex_to_tcr else "NA" for l in adata.obs["library_id"].values])

def add_vdj_lib_ids_to_integrated_data(tissues_dir: str):
    gex_to_bcr, gex_to_tcr = get_gex_lib_to_vdj_lib_mapping()

    files = os.listdir(tissues_dir)
    files = [f for f in files if f.endswith(".h5ad")]
    for file in files:
        file_path = os.path.join(tissues_dir, file)
        adata = anndata.read_h5ad(file_path)
        add_vdj_lib_ids_to_adata(adata, gex_to_bcr, gex_to_tcr)
        adata.write(file_path, compression="lzf")
        del adata
        gc.collect()

def report_vdj_vs_cell_label_metrics_all_libs(ir_lib_type: str, tissue_adatas: List[AnnData]):
    if ir_lib_type not in ["BCR", "TCR"]:
        raise ValueError("Unsupported lib_type: {}. Must be one of: BCR, TCR".format(ir_lib_type))
    is_b = ir_lib_type == "BCR"
    ir_libs = list(get_vdj_lib_to_gex_lib_mapping()[0 if is_b else 1].keys())
    samples = read_immune_aging_sheet("Samples")
    first = True
    for lib in ir_libs:
        # get some other metadata associated with this lib
        col_name = "{} lib".format(ir_lib_type)
        samples["Pick?"] = samples[col_name].fillna('').apply(lambda x: "Yes" if lib in x.split(",") else "No")
        idx = samples["Pick?"] == "Yes"
        sites = samples[idx]["Site"]
        donors = samples[idx]["Donor ID"]
        organs = set(samples[idx]["Organ"])
        if len(set(sites)) > 1:
            raise ValueError("More than one site was found for lib id {} of type {}".format(lib, ir_lib_type))
        if len(set(donors)) > 1:
            raise ValueError("More than one donor was found for lib id {} of type {}".format(lib, ir_lib_type))
        site, donor = set(sites).pop(), set(donors).pop()
        # get adatas
        adatas = []
        for tdata in tissue_adatas:
            lib_id_col = "{}cr_library_id".format("b" if is_b else "t")
            if lib_id_col not in tdata.obs.columns.values:
                raise ValueError("Make sure all tissue adata objects have {}".format(lib_id_col))
            adata_s = tdata[tdata.obs[lib_id_col] == lib, :].copy()
            if len(adata_s) > 0:
                adatas.append(adata_s)
        if len(adatas) == 0:
            print("No cells found for library {}".format(lib))
            continue
        elif len(adatas) == 1:
            adata = adatas[0]
        else:
            adata = adatas[0].concatenate(adatas[1:], join="outer")
        report_vdj_vs_cell_label_metrics(adata, lib, ir_lib_type, site, donor, organs, get_csv=True, skip_header=not first, skip_debug_print=True)
        if first:
            first = False
        del adata
        del adatas
        gc.collect()

def report_vdj_vs_cell_label_metrics(
    adata: AnnData,
    ir_lib_id: str,
    ir_lib_type: str,
    site: str,
    donor: str,
    organs,
    get_csv: bool = False,
    skip_header: bool = False,
    skip_debug_print: bool = False
):
    def debug_print(msg: str):
        if not skip_debug_print:
            print(msg)

    ct_key_high = "celltypist_majority_voting.Immune_All_High.totalvi.leiden_resolution_2.0"
    ct_key_low = "celltypist_majority_voting.Immune_All_Low.totalvi.leiden_resolution_2.0"

    # Define known cell label categories
    # high hierarchy aka low resolution
    known_b_cells_high = ["B cells", "B-cell lineage", "Plasma cells"]
    known_t_cells_high = ["Double-negative thymocytes", "Double-positive thymocytes", "ETP", "T cells"]
    known_non_b_non_t_cells_high = [
        "DC",
        "DC precursor",
        "Endothelial cells",
        "Epithelial cells",
        "Erythrocytes",
        "Erythroid",
        "Fibroblasts",
        "Granulocytes",
        "ILC",
        "Monocyte precursor",
        "Monocytes",
        "HSC/MPP",
        "Macrophages",
        "Mast cells",
        "Megakaryocyte precursor",
        "Megakaryocytes/platelets",
        "Early MK",
        "MNP",
        "Mono-mac",
        "Myelocytes",
        "pDC",
        "pDC precursor",
        "Promyelocyte"
    ]
    known_non_b_cells_high = known_t_cells_high + known_non_b_non_t_cells_high
    known_non_t_cells_high = known_b_cells_high + known_non_b_non_t_cells_high
    # low hierarchy aka high resolution
    known_b_cells_low = ["Cycling B cells"]
    known_t_cells_low = ["Cycling gamma-delta T cells", "Cycling T cells", "Early lymphoid/T lymphoid"]
    known_non_b_non_t_cells_low = [
        "Cycling DCs",
        "Cycling monocytes",
        "Cycling NK cells"
    ]
    known_non_b_cells_low = known_t_cells_low + known_non_b_non_t_cells_low
    known_non_t_cells_low = known_b_cells_low + known_non_b_non_t_cells_low

    def run_analysis(lib_id: str, receptor_type: str):
        if receptor_type not in ["BCR", "TCR"]:
            raise ValueError("Unknown receptor_type passed: {}".format(receptor_type))
        cell_type = "b" if receptor_type == "BCR" else "t"
        is_b = receptor_type == "BCR"
        debug_print("Running for {} of type {}...".format(lib_id, receptor_type))

        # Look at labels of cells that have receptors of the given type 
        debug_print("Looking at labels of cells that have {} cell receptors...".format(cell_type))
        ir_cells_by_receptor = adata[adata.obs["{}-has_ir".format(receptor_type)] == "True", :].copy()
        num_ir_cells_by_receptor_1 = len(ir_cells_by_receptor)
        debug_print("Examining high hierarchy cell type labels...")
        ir_cells_by_type = ir_cells_by_receptor.obs[ct_key_high]
        # known b/t cell types high
        indices = ir_cells_by_type.isin(known_b_cells_high if is_b else known_t_cells_high)
        num_known_ir_cells_high = indices.sum()
        ir_cells_by_type = ir_cells_by_type[~indices]
        # known non b/t cell types high
        indices = ir_cells_by_type.isin(known_non_b_cells_high if is_b else known_non_t_cells_high)
        num_known_non_ir_cells_high = indices.sum()
        ir_cells_by_type = ir_cells_by_type[~indices]
        if len(ir_cells_by_type) == 0:
            debug_print("No cells left to examine in low hierarchy")
            # continue anyway - all the steps below will just return 0 in this case
        else:
            debug_print("Examining low hierarchy cell type labels...")
        # low hierarchy
        ir_cells_by_type_low = ir_cells_by_receptor[ir_cells_by_type.index, :].obs[ct_key_low]
        # known b/t cell types low
        indices = ir_cells_by_type_low.isin(known_b_cells_low if is_b else known_t_cells_low)
        num_known_ir_cells_low = indices.sum()
        ir_cells_by_type_low = ir_cells_by_type_low[~indices]
        # known non b/t cell types low
        indices = ir_cells_by_type_low.isin(known_non_b_cells_low if is_b else known_non_t_cells_low)
        num_known_non_ir_cells_low = indices.sum()
        if num_ir_cells_by_receptor_1 == 0:
            pct_known_ir_cells_high_n_low, pct_known_non_ir_cells_high_n_low = -1, -1
        else:
            pct_known_ir_cells_high_n_low = ((num_known_ir_cells_high + num_known_ir_cells_low) * 100) / num_ir_cells_by_receptor_1
            pct_known_non_ir_cells_high_n_low = ((num_known_non_ir_cells_high + num_known_non_ir_cells_low) * 100) / num_ir_cells_by_receptor_1

        draw_separator_line()

        # Look at presence of b/t cell receptors for cells that have b/t cell labels (will give us a hint about V(D)J library seq. saturation)
        debug_print("Looking at presence of {}'s for cells that have {} cell labels ...".format(receptor_type, cell_type))
        # known b/t cell types high and low
        indices_high = adata.obs[ct_key_high].isin(known_b_cells_high if is_b else known_t_cells_high)
        indices_low = adata.obs[ct_key_low].isin(known_b_cells_low if is_b else known_t_cells_low)
        indices = indices_high | indices_low
        ir_cells_by_type = adata[indices, :].copy()
        num_ir_cells_by_type = len(ir_cells_by_type)
        # look at b/t cell receptors
        ir_cells_by_receptor = ir_cells_by_type[ir_cells_by_type.obs["{}-has_ir".format(receptor_type)] == "True", :].copy()
        num_ir_cells_by_receptor_2 = len(ir_cells_by_receptor)
        # num_ir_cells_by_type can be 0 e.g. in the case of SKIN
        pct_ir_cells_by_receptor_2 = -1 if num_ir_cells_by_type == 0 else (num_ir_cells_by_receptor_2 * 100) / num_ir_cells_by_type
        # look for any cell receptors of the "opposite" cell type
        key = "TCR-has_ir" if is_b else "BCR-has_ir"
        opposite_cells_by_receptor = ir_cells_by_type[ir_cells_by_type.obs[key] == "True", :].copy()
        num_opposite_cells_by_receptor = len(opposite_cells_by_receptor)
        pct_opposite_cells_by_receptor = -1 if num_ir_cells_by_type == 0 else (num_opposite_cells_by_receptor * 100) / num_ir_cells_by_type

        # Look at presence of b/t cell receptors for cells that DO NOT have b/t cell labels  (will give us a hint about false positive receptors)
        debug_print("Looking at presence of {}'s for cells that DO NOT have {} cell labels ...".format(receptor_type, cell_type))
        # known non b/t cell types high and low
        indices_high = adata.obs[ct_key_high].isin(known_b_cells_high if is_b else known_t_cells_high)
        indices_low = adata.obs[ct_key_low].isin(known_b_cells_low if is_b else known_t_cells_low)
        indices = indices_high | indices_low
        non_ir_cells_by_type = adata[~indices, :].copy()
        num_non_ir_cells_by_type = len(non_ir_cells_by_type)
        # look at b/t cell receptors
        ir_cells_by_receptor = non_ir_cells_by_type[non_ir_cells_by_type.obs["{}-has_ir".format(receptor_type)] == "True", :].copy()
        num_ir_cells_by_receptor_3 = len(ir_cells_by_receptor)
        pct_ir_cells_by_receptor_3 = -1 if num_non_ir_cells_by_type == 0 else (num_ir_cells_by_receptor_3 * 100) / num_non_ir_cells_by_type

        draw_separator_line()

        debug_print("Results:")
        if not get_csv:
            print("Of {} {} cells by receptor:".format(num_ir_cells_by_receptor_1, cell_type))
            print("\t{} cells are of {} type in high hierarchy".format(num_known_ir_cells_high, cell_type))
            print("\t{} cells are not of {} cell type in high hierarchy".format(num_known_non_ir_cells_high, cell_type))
            print("\t{} cells are of {} type in low hierarchy".format(num_known_ir_cells_low, cell_type))
            print("\t{} cells are not of {} cell type in low hierarchy".format(num_known_non_ir_cells_low, cell_type))

            print("Of {} {} cells by type (high or low):".format(num_ir_cells_by_type, cell_type))
            print("\t{} cells have {} cell receptors".format(num_ir_cells_by_receptor_2, cell_type))
            print("\t{} cells have {} cell receptors".format(num_opposite_cells_by_receptor, "t" if is_b else "b"))

            print("Of {} non {} cells by type (high or low):".format(num_non_ir_cells_by_type, cell_type))
            print("\t{} cells have {} cell receptors".format(num_ir_cells_by_receptor_3, cell_type))

            print("cell label keys used: {}, {}".format(ct_key_high, ct_key_low))
        elif get_csv:
            csv_rows = [
                {
                    "Lib Id": lib_id,
                    "Site": site,
                    "Donor ID": donor,
                    "Organ(s)": ",".join([str(o) for o in organs]),

                    "# {} cells by receptor".format(cell_type): num_ir_cells_by_receptor_1,
                    "Out of {} cells by receptor: # {} cells (high)".format(cell_type, cell_type): num_known_ir_cells_high,
                    "Out of {} cells by receptor: # {} cells (low)".format(cell_type, cell_type): num_known_ir_cells_low,
                    "Out of {} cells by receptor: pct {} cells (high&low)".format(cell_type, cell_type): pct_known_ir_cells_high_n_low,

                    "Out of {} cells by receptor: # non {} cells (high)".format(cell_type, cell_type): num_known_non_ir_cells_high,
                    "Out of {} cells by receptor: # non {} cells (low)".format(cell_type, cell_type): num_known_non_ir_cells_low,
                    "Out of {} cells by receptor: pct non {} cells (high&low)".format(cell_type, cell_type): pct_known_non_ir_cells_high_n_low,

                    "# {} cells by type".format(cell_type): num_ir_cells_by_type,
                    "Out of {} cells by type: # cells w/ {} cr".format(cell_type, cell_type): num_ir_cells_by_receptor_2,
                    "Out of {} cells by type: pct cells w/ {} cr (aka seq. saturation)".format(cell_type, cell_type): pct_ir_cells_by_receptor_2,
                    "Out of {} cells by type: # cells w/ {} cr".format(cell_type, "t" if is_b else "b"): num_opposite_cells_by_receptor,
                    "Out of {} cells by type: pct cells w/ {} cr".format(cell_type, "t" if is_b else "b"): pct_opposite_cells_by_receptor,

                    "# non {} cells by type".format(cell_type): num_non_ir_cells_by_type,
                    "Out of non {} cells by type: # cells w/ {} cr".format(cell_type, cell_type): num_ir_cells_by_receptor_3,
                    "Out of non {} cells by type: pct cells w/ {} cr (aka false positive)".format(cell_type, cell_type): pct_ir_cells_by_receptor_3,
                }
            ]
            csv_file = io.StringIO()
            field_names = list(csv_rows[0].keys())
            writer = csv.DictWriter(csv_file, fieldnames=field_names)
            if not skip_header:
                writer.writeheader()
            writer.writerows(csv_rows)
            print(csv_file.getvalue())
            csv_file.close()

    run_analysis(ir_lib_id, ir_lib_type)

def gather_extra_info_for_ir_libs(
    lib_type: str,
    libs_csv_path: str,
    all_ir_lib_metrics_csv_path: str,
    all_gex_lib_metrics_csv_path: str,
    working_dir: str,
    s3_access_file: str,
):
    if lib_type not in ["BCR", "TCR"]:
        raise ValueError(f"Unknown lib type: {lib_type}")

    libs_df = pd.read_csv(libs_csv_path, header=None)
    libs = libs_df.values.squeeze().tolist()

    ir_lib_metrics_df = pd.read_csv(all_ir_lib_metrics_csv_path) 
    ir_lib_metrics_df = ir_lib_metrics_df[ir_lib_metrics_df["Lib Type"] == lib_type]
    ir_lib_metrics_df = ir_lib_metrics_df.set_index("Lib ID")

    gex_lib_metrics_df = pd.read_csv(all_gex_lib_metrics_csv_path) 
    gex_lib_metrics_df = gex_lib_metrics_df.set_index("Lib ID")

    ir_to_gex = get_vdj_lib_to_gex_lib_mapping()[0 if lib_type == "BCR" else 1]

    # define the csv headers
    CSV_HEADER_LIB_ID: str = "Lib Id"
    CSV_HEADER_ESTIMATED_CELLS_CELLRANGER: str = "Estimated # of Cells (cellranger)"
    CSV_HEADER_MEAN_RPPC_CELLRANGER: str = "Mean Read Pairs per Cell (cellranger)"
    CSV_HEADER_NUM_RP_CELLRANGER: str = "Number of Read Pairs (cellranger)"
    CSV_HEADER_MEAN_URPPC_CELLRANGER: str = "Mean Used Read Pairs per Cell (cellranger)"
    CSV_HEADER_GEX_LIB_ID: str = "Gex Lib Id"
    CSV_HEADER_ESTIMATED_CELLS_IN_GEX_LIB_CELLRANGER: str = "Estimated # of Cells in Gex Lib (cellranger)"
    CSV_HEADER_NUM_CELLS_POST_QC: str = "# of Cells Post QC"
    CSV_HEADER_NUM_GEX_CELLS_POST_QC: str = "# of Cells Post QC in Gex Lib"
    csv_fields = [
        CSV_HEADER_LIB_ID,
        CSV_HEADER_ESTIMATED_CELLS_CELLRANGER,
        CSV_HEADER_MEAN_RPPC_CELLRANGER,
        CSV_HEADER_NUM_RP_CELLRANGER,
        CSV_HEADER_NUM_RP_CELLRANGER,
        CSV_HEADER_MEAN_URPPC_CELLRANGER,
        CSV_HEADER_GEX_LIB_ID,
        CSV_HEADER_ESTIMATED_CELLS_IN_GEX_LIB_CELLRANGER,
        CSV_HEADER_NUM_CELLS_POST_QC,
        CSV_HEADER_NUM_GEX_CELLS_POST_QC,
    ]

    csv_rows = []
    for lib in libs:
        # it is expected that a few libs are not presented in our
        # ir_lib_metrics_df, since they were sequenced for vdj
        # recently and we are yet to extract their sequencing info
        # just skip those for now. we'll deal with them separately
        if lib not in ir_lib_metrics_df.index:
            print("lib {} not found in ir_lib_metrics_df".format(lib))
            row = {f: lib if f == CSV_HEADER_LIB_ID else "" for f in csv_fields}
            csv_rows.append(row)
            continue

        gex_lib = ir_to_gex[lib]
        row = {
            CSV_HEADER_LIB_ID: lib,

            CSV_HEADER_ESTIMATED_CELLS_CELLRANGER: ir_lib_metrics_df["Estimated Number of Cells"].loc[lib],
            CSV_HEADER_MEAN_RPPC_CELLRANGER: ir_lib_metrics_df["Mean Read Pairs per Cell"].loc[lib],
            CSV_HEADER_NUM_RP_CELLRANGER: ir_lib_metrics_df["Number of Read Pairs"].loc[lib],
            CSV_HEADER_MEAN_URPPC_CELLRANGER: ir_lib_metrics_df["Mean Used Read Pairs per Cell"].loc[lib],

            CSV_HEADER_GEX_LIB_ID: gex_lib,
            CSV_HEADER_ESTIMATED_CELLS_IN_GEX_LIB_CELLRANGER: gex_lib_metrics_df["Estimated Number of Cells"].loc[gex_lib],
        }

        donor_id = ir_lib_metrics_df["Donor ID"].loc[lib]
        def add_info_from_processed_lib_object(library_type, library_id):
            # really hacky way to work around the fact that we don't have seq runs
            # at this layer. Almost all donors have 001 as the seq run, but a few
            # have 002 or 003 so try for those too
            for seq_run in ["001", "002", "003"]:
                file_name_partial = "{}_{}_{}_{}".format(donor_id, seq_run, library_type, library_id)
                s3_processed_lib_path = "s3://immuneaging/processed_libraries/{}".format(file_name_partial)
                version = get_latest_object_version(s3_access_file, s3_processed_lib_path)
                file_name = "{}.processed.{}.h5ad".format(file_name_partial, version)
                sync_cmd = 'aws s3 sync --no-progress {}/{}/ {} --exclude "*" --include {}'.format(
                    s3_processed_lib_path, version, working_dir, file_name
                )

                print("sync_cmd: {}".format(sync_cmd))
                print("aws response: {}\n".format(os.popen(sync_cmd).read()))
                adata_file = os.path.join(working_dir, file_name)
                if not os.path.isfile(adata_file):
                    print("Failed to download file {} from S3 using seq_run: {}".format(file_name, seq_run))
                else:
                    break # don't need to try the other seq runs
            
            col_name = CSV_HEADER_NUM_CELLS_POST_QC if library_type != "GEX" else CSV_HEADER_NUM_GEX_CELLS_POST_QC
            if not os.path.isfile(adata_file):
                print("Failed to download file {} from S3.".format(file_name, seq_run))
                row[col_name] = -1
            else:
                adata = anndata.read_h5ad(adata_file)
                row[col_name] = adata.n_obs
                os.remove(adata_file)

        add_info_from_processed_lib_object(lib_type, lib)
        add_info_from_processed_lib_object("GEX", gex_lib)

        csv_rows.append(row)

    csv_file = io.StringIO()
    field_names = list(csv_rows[0].keys())
    writer = csv.DictWriter(csv_file, fieldnames=field_names)
    writer.writeheader()
    writer.writerows(csv_rows)
    print(csv_file.getvalue())
    csv_file.close()

def get_ir_gex_intersection(
    ir_type,
    s3_access_file,
    working_dir,
    save_csv_dir,
    log_file_dir,
    only_donors: Optional[List[str]] = None
):
    unique_artifact_prefix = f"{ir_type}_{str(uuid.uuid1())}"
    log_file_path = f"{log_file_dir}/vdj_seq_sat_logs_{unique_artifact_prefix}.log"
    logger = SimpleLogger(filename = log_file_path)

    irs = get_vdj_lib_to_gex_lib_mapping()[0 if ir_type == "BCR" else 1]
    df = pd.DataFrame(columns=["donor_id", "ir_type", "ir_lib", "gex_lib", "ir_gex_diff_pct", "ir_gex_pre_qc_diff_pct"])
    samples = read_immune_aging_sheet("Samples")

    for ir_id,gex_id in irs.items():
        donor_id = get_donor_id_for_lib(ir_type, ir_id, samples)
        if only_donors is not None and donor_id not in only_donors:
            continue
        adata_ir = read_library(ir_type, ir_id, s3_access_file, working_dir, "processed", logger, remove_adata=False, donor_id=donor_id)
        adata_gex = read_library("GEX", gex_id, s3_access_file, working_dir, "processed", logger, remove_adata=False, donor_id=donor_id)
        adata_gex_pre_qc = read_library(
            "GEX",
            gex_id,
            s3_access_file,
            working_dir,
            "aligned",
            logger,
            remove_adata=False,
            donor_id=donor_id
        )

        if adata_ir is None or adata_gex is None or adata_gex_pre_qc is None:
            logger.add_to_log(f"❌❌ error. ir lib: {ir_id}, gex lib: {gex_id}")
            new_row = {
                "donor_id": donor_id,
                "ir_type": ir_type,
                "ir_lib": ir_id,
                "gex_lib": gex_id,
                "ir_gex_diff_pct": "-1",
                "ir_gex_pre_qc_diff_pct": "-1",
            }
            df = df.append(new_row, ignore_index=True)
            continue
        # add the gex library ID to the cell barcode name for the aligned lib
        adata_gex_pre_qc.obs_names = adata_gex_pre_qc.obs_names + "_" + gex_id
        ir_gex_diff = len(np.setdiff1d(adata_ir.obs.index, adata_gex.obs.index))
        ir_gex_diff_pct = (ir_gex_diff/len(adata_ir.obs.index)) * 100
        # same but pre gex qc
        ir_gex_pre_qc_diff = len(np.setdiff1d(adata_ir.obs.index, adata_gex_pre_qc.obs.index))
        ir_gex_pre_qc_diff_pct = (ir_gex_pre_qc_diff/len(adata_ir.obs.index)) * 100
        logger.add_to_log(f"👉👉 ir lib: {ir_id}, ir type: {ir_type}, gex lib: {gex_id}, ir_gex_diff_pct: {ir_gex_diff_pct:.2f}, ir_gex_pre_qc_diff_pct: {ir_gex_pre_qc_diff_pct:.2f}")
        new_row = {
            "donor_id": donor_id,
            "ir_type": ir_type,
            "ir_lib": ir_id,
            "gex_lib": gex_id,
            "ir_gex_diff_pct": f"{ir_gex_diff_pct:.2f}",
            "ir_gex_pre_qc_diff_pct": f"{ir_gex_pre_qc_diff_pct:.2f}",
        }
        df = df.append(new_row, ignore_index=True)
    csv_path = f"{save_csv_dir}/vdj_seq_sat_results_{unique_artifact_prefix}.csv"
    df.to_csv(csv_path)
    logger.add_to_log(f"✅✅ results: {csv_path}")
    return df

def report_vdj_lib_ss_and_fp_metrics_for_all_libs(
        ir_lib_type: str,
        b_adata_path: str,
        t_adata_path: str,
        myeloid_adata_path: str,
        other_adata_path: str,
        csv_file_path: Optional[str] = None,
        only_donors: Optional[List[str]] = None,
    ):
    assert ir_lib_type in ["BCR", "TCR"]
    cols = ['BCR-has_ir', 'TCR-has_ir', 'bcr_library_id', 'tcr_library_id']
    b_cells = anndata.read_h5ad(b_adata_path, backed=True).obs[cols]
    t_cells = anndata.read_h5ad(t_adata_path, backed=True).obs[cols]
    myeloid_cells = anndata.read_h5ad(myeloid_adata_path, backed=True).obs[cols]
    other_cells = anndata.read_h5ad(other_adata_path, backed=True).obs[cols]

    def run_for_lib(
        lib_id: str,
        ir_lib_type: str,
        site: str,
        donor: str,
        organs: List[str],
        skip_header: bool = False,
    ):
        cell_type = "b" if ir_lib_type == "BCR" else "t"
        print("Running for {} of type {}...".format(lib_id, ir_lib_type))
        lib_id_str = f"{cell_type}cr_library_id"
        has_ir_str = f"{ir_lib_type}-has_ir"
        # get the subset of the cells that have the given library ID
        b_cells_lib = b_cells[b_cells[lib_id_str] == lib_id]
        t_cells_lib = t_cells[t_cells[lib_id_str] == lib_id]
        m_cells_lib = myeloid_cells[myeloid_cells[lib_id_str] == lib_id]
        o_cells_lib = other_cells[other_cells[lib_id_str] == lib_id]

        # Look at the presence of b cell receptors for b cells (will give us a hint about BCR library seq. saturation)
        # same for t depending on the ir_lib_type
        print("Looking at the presence of {}'s for {} cells...".format(ir_lib_type, cell_type))
        cells = b_cells_lib if ir_lib_type == "BCR" else t_cells_lib
        num_cells = len(cells)
        receptors = cells[cells[has_ir_str] == "True"]
        num_receptors = len(receptors)
        # num_cells can be 0 e.g. in the case of SKIN
        # it could also be that the corresponding GEX library failed
        pct_receptors = -1 if num_cells == 0 else (num_receptors * 100) / num_cells
        # look for any cell receptors of the "opposite" cell type
        opposite_ir_lib_type = "TCR" if ir_lib_type == "BCR" else "BCR"
        opposite_receptors = cells[cells["{}-has_ir".format(opposite_ir_lib_type)] == "True"]
        num_opposite_receptors = len(opposite_receptors)
        pct_opposite_receptors = -1 if num_cells == 0 else (num_opposite_receptors * 100) / num_cells

        # Look at the presence of b cell receptors for non-b cells (will give us a hint about false positive receptors)
        # same for t depending on the ir_lib_type
        print("Looking at the presence of {}'s for non-{} cells...".format(ir_lib_type, cell_type))
        opposite_type_cells = t_cells_lib if ir_lib_type == "BCR" else b_cells_lib
        num_cells_other_types = len(opposite_type_cells) + len(m_cells_lib) + len(o_cells_lib)
        num_receptors_other_types = len(opposite_type_cells[opposite_type_cells[has_ir_str] == "True"])
        num_receptors_other_types += len(m_cells_lib[m_cells_lib[has_ir_str] == "True"])
        num_receptors_other_types += len(o_cells_lib[o_cells_lib[has_ir_str] == "True"])
        # num_cells_other_types could be 0 if for instance the corresponding GEX library failed
        pct_receptors_other_types = -1 if num_cells_other_types == 0 else (num_receptors_other_types * 100) / num_cells_other_types

        print("✅✅ Results:")
        csv_rows = [
            {
                "Lib Id": lib_id,
                "Site": site,
                "Donor ID": donor,
                "Organ(s)": ",".join([str(o) for o in organs]),

                f"# {cell_type} cells": num_cells,
                f"# {cell_type} cells with {ir_lib_type}": num_receptors,
                f"(# {cell_type} cells with {ir_lib_type} * 100)/# {cell_type} cells (aka seq. saturation)": pct_receptors,

                f"# {cell_type} cells with {opposite_ir_lib_type}": num_opposite_receptors,
                f"(# {cell_type} cells with {opposite_ir_lib_type} * 100)/# {cell_type} cells": pct_opposite_receptors,

                f"# non-{cell_type} cells": num_cells_other_types,
                f"# non-{cell_type} cells with {ir_lib_type}": num_receptors_other_types,
                f"(# non-{cell_type} cells with {ir_lib_type} * 100)/# non-{cell_type} cells (aka false positive)": pct_receptors_other_types,
            }
        ]
        if csv_file_path is not None:
            field_names = list(csv_rows[0].keys())
            with open(csv_file_path, "a") as f:
                writer = csv.DictWriter(f, fieldnames=field_names)
                if not skip_header:
                    writer.writeheader()
                writer.writerows(csv_rows)
        else:
            csv_file = io.StringIO()
            field_names = list(csv_rows[0].keys())
            writer = csv.DictWriter(csv_file, fieldnames=field_names)
            if not skip_header:
                writer.writeheader()
            writer.writerows(csv_rows)
            print(csv_file.getvalue())
            csv_file.close()

    ir_libs = list(get_vdj_lib_to_gex_lib_mapping()[0 if ir_lib_type == "BCR" else 1].keys())
    samples = read_immune_aging_sheet("Samples")
    first = True
    for lib in ir_libs:
        # get some other metadata associated with this lib
        col_name = "{} lib".format(ir_lib_type)
        samples["Pick?"] = samples[col_name].fillna('').apply(lambda x: "Yes" if lib in x.split(",") else "No")
        idx = samples["Pick?"] == "Yes"
        sites = samples[idx]["Site"]
        donors = samples[idx]["Donor ID"]
        organs = set(samples[idx]["Organ"])
        if len(set(sites)) > 1:
            raise ValueError("More than one site was found for lib id {} of type {}".format(lib, ir_lib_type))
        if len(set(donors)) > 1:
            raise ValueError("More than one donor was found for lib id {} of type {}".format(lib, ir_lib_type))
        site, donor = set(sites).pop(), set(donors).pop()
        if only_donors is not None and donor not in only_donors:
            continue
        # run for lib
        run_for_lib(lib, ir_lib_type, site, donor, organs, skip_header=not first)
        if first:
            first = False