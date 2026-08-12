from __future__ import absolute_import

from pmultiqc.modules.common.logging import get_logger
from pmultiqc.modules.common.plots.general import (
    stat_pep_intensity
)
from pmultiqc.modules.common.common_utils import (
    cal_miss_cleavages
)


# Initialise the module logger via centralized logger
log = get_logger("pmultiqc.modules.qpx.qpx_utils")


def calculate_run_stat(psm_df_sub, proteins):
    """Calculate statistics for a specific run."""
    peptides = set(psm_df_sub["peptidoform"])

    # TODO need protein_accessions in psm.parquet
    unique_peptides = set()

    # 'modifications' is optional: writers may omit it, and the reader only loads
    # columns the file actually has.
    if "modifications" in psm_df_sub.columns:
        modified_pep = set(
            psm_df_sub[psm_df_sub["modifications"].notna()]["peptidoform"]
        )
    else:
        modified_pep = set()

    stat_run = {
        "protein_num": len(proteins),
        "peptide_num": len(peptides),
        "unique_peptide_num": len(unique_peptides),
        "modified_peptide_num": len(modified_pep)
    }

    data_per_run = {
        "proteins": proteins,
        "peptides": peptides,
        "unique_peptides": unique_peptides,
        "modified_peps": modified_pep
    }

    return stat_run, data_per_run


def get_unimod_mod_qpx(modifis, unimod_data):

    if modifis is None or len(modifis) == 0:
        return "Unmodified"

    mod_list = []

    for mod_dict in modifis:
        acc = mod_dict.get('accession')
        if not acc:
            continue

        unimod_entry = unimod_data.get_by_accession(acc.upper())

        if unimod_entry is not None:
            mod_list.append(unimod_entry.get_name())
        else:
            mod_list.append(acc.upper())

    return ",".join(set(mod_list)) if mod_list else "Unmodified"


def get_pep_intensity(pep_table, sdrf_file_df):

    pep_intensity_by_run = dict()
    pep_intensity_by_sample = dict()

    pep_table = pep_table.dropna(subset=["intensities"]).explode("intensities")

    pep_table["channel"] = pep_table["intensities"].str.get("label")
    pep_table["intensity"] = pep_table["intensities"].str.get("intensity").astype(float)

    pep_table = pep_table[pep_table["intensity"] > 0].drop(columns=["intensities"])

    tmt_to_label = {
        'TMT126': 1, 'TMT127N': 2, 'TMT127C': 3, 'TMT128N': 4, 'TMT128C': 5, 
        'TMT129N': 6, 'TMT129C': 7, 'TMT130N': 8, 'TMT130C': 9, 'TMT131': 10,
        'TMT131C': 11, 'TMT132N': 12, 'TMT132C': 13, 'TMT133N': 14, 'TMT133C': 15, 'TMT134N': 16
    }

    pep_table["Label"] = pep_table["channel"].map(tmt_to_label)
    pep_table["Label"] = pep_table["Label"].fillna(1).astype(int).astype(str)

    if sdrf_file_df.empty:
        log.info("No SDRF found, displaying pep_intensity according to run")
    else:
        sdrf_subset = sdrf_file_df[["Label", "Run", "Sample"]].drop_duplicates()
        sdrf_subset["Label"] = sdrf_subset["Label"].astype(str)

        pep_table = pep_table.merge(
            sdrf_subset,
            left_on=["run", "Label"],
            right_on=["Run", "Label"],
            how="left"
        )

        for sample_id, group in pep_table.groupby("Sample", sort=True):
            sample_name_str = f"Sample {str(sample_id)}"

            pep_intensity_by_sample[sample_name_str] = stat_pep_intensity(group["intensity"])

    for run_name, group in pep_table.groupby("run"):
        run_name = str(run_name)
        pep_intensity_by_run[run_name] = stat_pep_intensity(group["intensity"])

    return [pep_intensity_by_run, pep_intensity_by_sample]


def get_missed_cleavages(psm_df, run_df, sdrf_file_df):

    enzyme = _get_enzyme_name(run_df)
    psm_df["missed_cleavages"] = psm_df["sequence"].apply(
        lambda seq: cal_miss_cleavages(seq, enzyme)
    )

    mc_by_run = {}
    mc_by_sample = {}

    for name, group in psm_df.groupby("run"):
        sc = group["missed_cleavages"].value_counts()
        mc_by_run[str(name)] = sc.to_dict()

    if sdrf_file_df is not None and not sdrf_file_df.empty:
        merged_df = psm_df.merge(
            right=sdrf_file_df[["Sample", "Run"]].drop_duplicates(),
            left_on="run",
            right_on="Run",
            how="left"
        )

        merged_df = merged_df.dropna(subset=["Sample"])

        for name, group in merged_df.groupby("Sample", sort=True):
            sc = group["missed_cleavages"].value_counts()
            mc_by_sample[f"Sample {str(name)}"] = sc.to_dict()

    return {
        "sdrf_samples": mc_by_sample,
        "ms_runs": mc_by_run,
    }


def _get_enzyme_name(df):

    enzyme_name = "Trypsin"

    if df is None or (hasattr(df, "empty") and df.empty):
        return enzyme_name

    col_name = "enzymes"
    if hasattr(df, "columns") and col_name in df.columns:
        valid_data = df[col_name].dropna()

        if not valid_data.empty:
            name = _first_enzyme_value(valid_data.iloc[0])
            if name:
                enzyme_name = name

    return enzyme_name


def _first_enzyme_value(cell):
    """Extract the enzyme name from an 'enzymes' cell.

    The column may hold a plain string, a list/array of names, or a struct such as
    {"name": "Trypsin", ...}, so unwrap those shapes instead of blindly indexing [0]
    (which turns the string "Trypsin" into "T").
    """
    if cell is None:
        return None

    if isinstance(cell, str):
        return cell.strip() or None

    if isinstance(cell, dict):
        for key in ("name", "enzyme", "enzyme_name"):
            if cell.get(key):
                return _first_enzyme_value(cell[key])
        return next((_first_enzyme_value(v) for v in cell.values() if v), None)

    # list, tuple, numpy array, pyarrow list scalar, ...
    try:
        if len(cell) > 0:
            return _first_enzyme_value(cell[0])
    except TypeError:
        return str(cell).strip() or None

    return None
