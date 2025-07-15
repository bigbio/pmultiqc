from pathlib import Path
from typing import List, Dict, Any, Union

import pandas as pd
import re
import numpy as np
import os

from pandas._typing import ReadCsvBuffer
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from ..common.file_utils import get_filename
from ..common.calc_utils import qualUniform
from ...logging import get_logger, Timer

# Initialize logger for this module
logger = get_logger("pmultiqc.modules.maxquant")


def read(
    file: Path | ReadCsvBuffer[bytes] | ReadCsvBuffer[str] | str,
    file_type: str = None,
    filter_type: str = None,
):
    """
    Read MaxQuant output files and perform initial processing.

    Args:
        file: Path to the csv file containing the information
        file_type: Type of file (protgroup, summary, etc.)
        filter_type: Type of filtering to apply (e.g., "Reverse")

    Returns:
        Processed pandas DataFrame
    """
    file_name = get_filename(file)
    with Timer(logger, f"Reading file {os.path.basename(file_name)}"):
        mq_data = pd.read_csv(file, sep="\t", low_memory=False)
        logger.info(f"Loaded {len(mq_data)} rows from {os.path.basename(file_name)}")

    # columns to lower
    mq_data.columns = mq_data.columns.str.lower()

    # Field correction for different versions of MaxQuant
    if "contaminant" in mq_data.columns:
        logger.debug("Renaming 'contaminant' column to 'potential contaminant'")
        mq_data = mq_data.rename(columns={"contaminant": "potential contaminant"})

    if filter_type == "Reverse":
        # Reverse
        # When marked with ‘+’, this particular protein group contains no protein,
        # made up of at least 50% of the peptides of the leading protein, with a
        # peptide derived from the reversed part of the decoy database.
        # These should be removed for further data analysis.
        # The 50% rule is in place to prevent spurious protein hits to erroneously flag the protein group as reverse.
        if "reverse" in mq_data.columns:
            reverse_count = mq_data[mq_data["reverse"] == "+"].shape[0]
            mq_data = mq_data[mq_data["reverse"] != "+"].reset_index(drop=True)
            logger.info(
                f"Filtered out {reverse_count} reverse entries, {len(mq_data)} entries remaining"
            )

    # proteinGroups.txt
    if file_type == "protgroup":
        logger.debug("Processing proteinGroups file")
        if "mol. weight [kda]" in mq_data.columns:
            intensity_cols = [
                col for col in mq_data.columns if re.search(r"intensity", col, re.IGNORECASE)
            ]
            logger.debug(f"Found {len(intensity_cols)} intensity columns")

            new_intensity_cols = [
                re.sub(r"intensity", "abind", col, flags=re.IGNORECASE) for col in intensity_cols
            ]

            for col, new_col in zip(intensity_cols, new_intensity_cols):
                mq_data[new_col] = mq_data[col] / mq_data["mol. weight [kda]"]
                mq_data = mq_data.copy()

            logger.debug(f"Created {len(new_intensity_cols)} abundance index columns")
        else:
            logger.error("The column 'Mol. weight [kDa]' could not be found!")
            raise ValueError("The column 'Mol. weight [kDa]' could not be found!")

    # summary.txt
    if file_type == "summary":
        logger.debug("Processing summary file")
        row_index = mq_data.notna().sum(axis=1) > 2
        row_index.iloc[-1] = False

        original_len = len(mq_data)
        mq_data = mq_data[row_index]
        logger.debug(f"Filtered summary file from {original_len} to {len(mq_data)} rows")

    return mq_data


# 1. proteinGroups.txt
def get_protegroups(file_path: str) -> Dict[str, dict]:
    """
    Process proteinGroups.txt file and extract various metrics.

    Args:
        file_path: Path to the proteinGroups.txt file

    Returns:
        Dictionary containing various metrics extracted from the proteinGroups file
    """
    with Timer(logger, "Processing proteinGroups data"):
        mq_data = read(file_path, "protgroup", filter_type="Reverse")
        logger.info(f"Processing proteinGroups data with {len(mq_data)} entries")

        # 'intensity'
        intensity_hlm_exp_cols = [
            col for col in mq_data.columns if re.match(r"^intensity [hlm] \S+$", col)
        ]
        intensity_exp_cols = [col for col in mq_data.columns if re.match(r"^intensity \S+$", col)]
        intensity_cols = [col for col in mq_data.columns if re.match(r"intensity$", col)]

        intensity_columns = []

        if len(intensity_hlm_exp_cols):
            intensity_columns = intensity_hlm_exp_cols
            logger.debug(f"Using {len(intensity_hlm_exp_cols)} HLM intensity columns")
        elif len(intensity_exp_cols):
            intensity_columns = intensity_exp_cols
            logger.debug(f"Using {len(intensity_exp_cols)} experiment intensity columns")
        elif len(intensity_cols):
            intensity_columns = intensity_cols
            logger.debug(f"Using {len(intensity_cols)} basic intensity columns")

        logger.debug(f"Found {len(intensity_columns)} total intensity columns")

        # 1: PG: ~Contaminants
        if "potential contaminant" in mq_data.columns:
            logger.debug("Processing contaminants")
            contaminant_percent_dict = pg_contaminants(mq_data, intensity_columns)
        else:
            logger.debug("No contaminant column found")
            contaminant_percent_dict = None

        # 2. PG: ~Intensity distribution
        logger.debug("Calculating intensity distribution")
        intensity_distr_dict = pg_intensity_distr(mq_data, intensity_columns)

        # 'LFQ intensity'
        logger.debug("Looking for LFQ intensity columns")
        lfq_intensity_hlm_exp_cols = [
            col for col in mq_data.columns if re.match(r"^lfq intensity [hlm] \S+$", col)
        ]
        lfq_intensity_exp_cols = [
            col for col in mq_data.columns if re.match(r"^lfq intensity \S+$", col)
        ]

        lfq_intensity_columns = []

        if len(lfq_intensity_hlm_exp_cols):
            lfq_intensity_columns = lfq_intensity_hlm_exp_cols
            logger.debug(f"Using {len(lfq_intensity_hlm_exp_cols)} HLM LFQ intensity columns")
        elif len(lfq_intensity_exp_cols):
            lfq_intensity_columns = lfq_intensity_exp_cols
            logger.debug(f"Using {len(lfq_intensity_exp_cols)} experiment LFQ intensity columns")

        logger.debug(f"Found {len(lfq_intensity_columns)} total LFQ intensity columns")

        # LFQ
        if lfq_intensity_columns:
            logger.debug("Calculating LFQ intensity distribution")
            lfq_intensity_distr = pg_intensity_distr(mq_data, lfq_intensity_columns)
        else:
            logger.debug("No LFQ intensity columns found")
            lfq_intensity_distr = None

        # Summary Table
        summary_table_dict = pg_summary_table(mq_data, lfq_intensity_columns, intensity_columns)

        # Number of Peptides identified Per Protein
        num_pep_per_protein_dict = peptide_per_protein(mq_data)

        # PCA
        if len(intensity_columns) > 1:
            logger.debug("Calculating raw intensity PCA")
            raw_intensity_pca = pg_pca(mq_data, intensity_columns)
        else:
            logger.debug("Not enough intensity columns for PCA")
            raw_intensity_pca = None

        if len(lfq_intensity_columns) > 1:
            logger.debug("Calculating LFQ intensity PCA")
            lfq_intensity_pca = pg_pca(mq_data, lfq_intensity_columns)
        else:
            logger.debug("Not enough LFQ intensity columns for PCA")
            lfq_intensity_pca = None

        result = {
            "pg_contaminant": contaminant_percent_dict,
            "pg_intensity_distri": intensity_distr_dict,
            # 'intensity_cols': intensity_columns,
            # 'lfq_intensity_cols': lfq_intensity_columns,
            # 'pg_data': mq_data,
            "pg_lfq_intensity_distri": lfq_intensity_distr,
            "raw_intensity_pca": raw_intensity_pca,
            "lfq_intensity_pca": lfq_intensity_pca,
            "protein_summary": summary_table_dict,
            "num_pep_per_protein_dict": num_pep_per_protein_dict,
        }

        logger.info("Completed processing proteinGroups data")
        return result


# 1-1. PG:~Contaminants
def pg_contaminants(mq_data: pd.DataFrame, intensity_cols: List[str]) -> dict[Any, dict[str, Any]]:
    """
    Calculate the percentage of contaminants in each group.

    Args:
        mq_data: DataFrame containing protein groups data
        intensity_cols: List of intensity column names

    Returns:
        Dictionary mapping group names to contaminant percentages
    """
    with Timer(logger, "Calculating contaminant percentages"):
        if any(column not in mq_data.columns for column in intensity_cols):
            logger.warning("Some intensity columns not found in data")
            return None

        logger.debug(f"Calculating total intensity for {len(intensity_cols)} columns")
        df1 = mq_data[intensity_cols].sum().to_frame().reset_index()
        df1.columns = ["group", "total_intensity"]

        contaminant_count = mq_data[mq_data["potential contaminant"] == "+"].shape[0]
        logger.debug(f"Found {contaminant_count} contaminant entries")

        df2 = (
            mq_data[mq_data["potential contaminant"] == "+"][intensity_cols]
            .sum()
            .to_frame()
            .reset_index()
        )
        df2.columns = ["group", "contaminant_total_intensity"]

        result_df = pd.merge(df1, df2, on="group", how="inner")
        result_df["contaminant_percent"] = (
            result_df["contaminant_total_intensity"] / result_df["total_intensity"] * 100.00
        )

        result_dict = dict()
        for k, v in dict(zip(result_df["group"], result_df["contaminant_percent"])).items():
            result_dict[k] = {"Potential Contaminants": v}

        logger.info(f"Calculated contaminant percentages for {len(result_dict)} groups")
        return result_dict


# 1-2. PG: ~Intensity distribution
def pg_intensity_distr(mq_data, intensity_cols):

    if any(column not in mq_data.columns for column in intensity_cols):
        return None

    # Take the logarithm and remove zero values
    def box_fun(raw_df):
        log_df = raw_df.apply(lambda col: col.map(lambda x: 1 if (pd.isna(x) or x == 0) else x))
        log_df = np.log2(log_df).reset_index(drop=True)
        log_df_dict = log_df.to_dict(orient="list")
        log_df_dict = {
            key: [value for value in values if value != 0] for key, values in log_df_dict.items()
        }
        return log_df_dict

    if "potential contaminant" in mq_data.columns:
        raw_df = mq_data[intensity_cols]
        contaminant_df = mq_data[mq_data["potential contaminant"] == "+"][intensity_cols]

        # Box plot data
        boxplot_data = [box_fun(raw_df), box_fun(contaminant_df)]
    else:
        raw_df = mq_data[intensity_cols]

        # Box plot data
        boxplot_data = [box_fun(raw_df), None]

    result_dict = {"box": boxplot_data}

    return result_dict


# 1-3.proteinGroups.txt: PCA
def pg_pca(pg_data: pd.DataFrame, cols_name: List[str]):

    if any(column not in pg_data.columns for column in cols_name):
        return None

    if "potential contaminant" in pg_data.columns:
        pg_data = pg_data[pg_data["potential contaminant"] != "+"].copy()

    pca_df = pg_data[cols_name].copy().T

    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(pca_df)

    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_data)

    pca_result_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"])
    pca_result_df["Group"] = pca_df.index

    pca_dict = {}
    for raw_name, group in pca_result_df.groupby("Group"):
        pca_dict[raw_name] = {"x": group.iloc[0, 0], "y": group.iloc[0, 1]}

    return pca_dict


# 1-4.proteinGroups.txt: Summary Table
def pg_summary_table(pg_data, lfq_intensity, intensity):

    if "potential contaminant" in pg_data.columns:
        pg_data = pg_data[pg_data["potential contaminant"] != "+"].copy()

    num_proteins_identified = len(pg_data)

    intensitys = lfq_intensity + intensity
    num_proteins_quantified = len(pg_data[pg_data[intensitys].gt(0).any(axis=1)])

    result = {
        "num_proteins_identified": num_proteins_identified,
        "num_proteins_quantified": num_proteins_quantified,
    }

    return result


# 1-5.proteinGroups.txt: Number of Peptides identified Per Protein
def peptide_per_protein(pg_data):

    if "peptide counts (all)" not in pg_data.columns:
        return None

    if "potential contaminant" in pg_data.columns:
        pg_data = pg_data[pg_data["potential contaminant"] != "+"].copy()

    count = pg_data["peptide counts (all)"].apply(lambda x: int(x.split(";")[0]))
    count = count[count > 0]

    count_df = count.value_counts().sort_index().reset_index()
    count_df.columns = ["peptide_count", 'frequency']

    top_n = 49
    if len(count_df) > top_n:
        top_count_df = count_df.iloc[:top_n].copy()
        rest = count_df.iloc[top_n:]
        other_sum = rest["frequency"].sum()
        threshold_value = count_df.iloc[top_n]["peptide_count"]
        top_count_df.loc[len(top_count_df)] = [f"{threshold_value}+", other_sum]
        final_count_df = top_count_df.copy()
    else:
        final_count_df = count_df.copy()

    final_count_df["percentage"] = (final_count_df["frequency"] / final_count_df["frequency"].sum()) * 100

    result = [
        {
            str(k): {"Frequency": v} for k, v in zip(final_count_df["peptide_count"], final_count_df["frequency"])
        },
        {
            str(k): {"Percentage": v} for k, v in zip(final_count_df["peptide_count"], final_count_df["percentage"])
        },
    ]

    return result


# 2. summary.txt
def get_summary(file_path: Union[Path, str]):
    """
    Process summary.txt file and extract MS/MS identification percentages.

    Args:
        file_path: Path to the summary.txt file

    Returns:
        Dictionary mapping raw file names to MS/MS identification percentages
    """
    with Timer(logger, "Processing summary data"):
        mq_data = read(file_path, "summary")
        logger.info(f"Processing summary data with {len(mq_data)} entries")

        if any(column not in mq_data.columns for column in ["raw file", "ms/ms identified [%]"]):
            logger.warning("Required columns not found in summary file")
            return None

        if all(mq_data["ms/ms identified [%]"] == 0):
            logger.warning("All MS/MS identification percentages are zero")
            return None

        msms_identified = dict()
        for k, v in dict(zip(mq_data["raw file"], mq_data["ms/ms identified [%]"])).items():
            msms_identified[k] = {"Identified Rate": v}

        logger.info(
            f"Extracted MS/MS identification percentages for {len(msms_identified)} raw files"
        )
        return msms_identified


# 3. evidence.txt
def get_evidence(file_path):
    """
    Process evidence.txt file and extract various metrics.

    Args:
        file_path: Path to the evidence.txt file

    Returns:
        Dictionary containing various metrics extracted from the evidence file
    """
    mq_data = read(file_path, filter_type="Reverse")

    # Count peptides and MS/MS spectra
    total_entries = len(mq_data)

    unique_peptides = (
        mq_data["modified sequence"].nunique() if "modified sequence" in mq_data.columns else 0
    )
    msms_count = mq_data["ms/ms count"].sum() if "ms/ms count" in mq_data.columns else 0

    logger.info(
        f"Processing evidence data with {total_entries} entries, {unique_peptides} unique peptides, and {msms_count} MS/MS spectra"
    )

    # Summary table
    summary_table_stat = evidence_summary_table(mq_data)

    if not all(column in mq_data.columns for column in ["type"]):
        logger.error('Missing required columns (#Type) in "evidence.txt"!')
        raise ValueError('Missing required columns (#Type) in "evidence.txt"!')

    mq_data["is_transferred"] = mq_data["type"] == "MULTI-MATCH"

    evidence_df = mq_data[mq_data["type"] != "MULTI-MATCH"].copy()
    evidence_df_tf = mq_data[mq_data["type"] == "MULTI-MATCH"].copy()

    logger.info(
        f"Found {len(evidence_df)} direct matches and {len(evidence_df_tf)} transferred matches"
    )

    # Top Contaminants per Raw file
    if "potential contaminant" in evidence_df.columns:
        logger.debug("Processing top contaminants")
        top_cont_dict = evidence_top_contaminants(evidence_df, top_n=5)
    else:
        logger.debug("No contaminant column found")
        top_cont_dict = None

    # peptide intensity distribution
    logger.debug("Calculating peptide intensity distribution")
    peptide_intensity_dict = evidence_peptide_intensity(evidence_df)

    # Distribution of precursor charges
    logger.debug("Calculating charge distribution")
    charge_counts_dict = evidence_charge_distribution(evidence_df)

    # Modifications per Raw file
    logger.debug("Processing modifications")
    modified_dict = evidence_modified(evidence_df)

    # rt_count_dict
    logger.debug("Calculating retention time counts")
    rt_count_dict = evidence_rt_count(evidence_df)

    # Peak width over RT
    logger.debug("Calculating peak width over retention time")
    peak_rt_dict = evidence_peak_width_rt(evidence_df)

    # Oversampling
    logger.debug("Calculating oversampling")
    oversampling_dict = evidence_oversampling(evidence_df)

    # Uncalibrated Mass Error
    logger.debug("Calculating uncalibrated mass error")
    uncalibrated_mass_error = evidence_uncalibrated_mass_error(evidence_df)

    # Calibrated Mass Error
    logger.debug("Calculating calibrated mass error")
    calibrated_mass_error = evidence_calibrated_mass_error(evidence_df)

    # Peptide ID count
    logger.debug("Counting peptide IDs")
    peptide_id_count = evidence_peptide_count(evidence_df, evidence_df_tf)

    # ProteinGroups count
    logger.debug("Counting protein groups")
    protein_group_count = evidence_protein_count(evidence_df, evidence_df_tf)

    # Delta Mass [Da]
    maxquant_delta_mass_da = evidence_delta_mass_da(evidence_df)

    result = {
        "top_contaminants": top_cont_dict,
        "peptide_intensity": peptide_intensity_dict,
        "charge_counts": charge_counts_dict,
        "modified_percentage": modified_dict,
        "rt_counts": rt_count_dict,
        # 'all_evidence': mq_data,
        "evidence_df": evidence_df,
        # 'evidence_df_tf': evidence_df_tf,
        "peak_rt": peak_rt_dict,
        "oversampling": oversampling_dict,
        "uncalibrated_mass_error": uncalibrated_mass_error,
        "calibrated_mass_error": calibrated_mass_error,
        "peptide_id_count": peptide_id_count,
        "protein_group_count": protein_group_count,
        "summary_stat": summary_table_stat,
        "maxquant_delta_mass_da": maxquant_delta_mass_da,
    }

    logger.info("Completed processing evidence data")
    return result


# HeatMap
def calculate_heatmap(evidence_df, oversampling, msms_missed_cleavages):

    if any(x is None for x in (evidence_df, oversampling, msms_missed_cleavages)):
        return None
    
    if any(
        column not in evidence_df.columns
        for column in ["potential contaminant", "intensity", "raw file", "retention time", "charge"]
    ):
        return None
    
    if evidence_df[evidence_df["potential contaminant"] == "+"].empty:
        logger.info("The evidence.txt file does not contain any contaminants")

    evidence_data = evidence_df.copy()

    # 8. Pep Missing Values
    global_peps = evidence_df["modified sequence"].unique()
    global_peps_count = len(global_peps)

    heatmap_dict = dict()
    for raw_file, group in evidence_data[
        [
            "potential contaminant", "intensity", "retention time",
            "raw file", "modified sequence"
            ]
        ].groupby("raw file"):

        # 1. Contaminants
        contaminant = 1 - (group[group["potential contaminant"] == "+"]["intensity"].sum() 
                           / group["intensity"].sum())
        
        # 2. Peptide Intensity
        peptide_intensity = np.minimum(1.0, np.nanmedian(group["intensity"]) / (2**23))

        # 8. Pep Missing Values
        pep_missing_values = np.minimum(
            1.0,
            len(set(global_peps) & set(group["modified sequence"].unique())) / global_peps_count
        )

        heatmap_dict[raw_file] = {
            "Contaminants": contaminant,
            "Peptide Intensity": peptide_intensity,
            "ID rate over RT": qualUniform(group["retention time"]),      # 6. ID rate over RT
            "Pep Missing Values": pep_missing_values,
        }

    # 4. Missed Cleavages
    missed_cleavages = {key: value["0"] / 100 for key, value in msms_missed_cleavages.items()}

    # 5. Missed Cleavages Var
    mc_median = np.median(list(missed_cleavages.values()))
    missed_cleavages_var = dict(
        zip(
            missed_cleavages.keys(),
            list(map(lambda v: 1 - np.abs(v - mc_median), missed_cleavages.values())),
        )
    )
    for raw_file in missed_cleavages.keys():
        heatmap_dict[raw_file]["Missed Cleavages"] = missed_cleavages[raw_file]
        heatmap_dict[raw_file]["Missed Cleavages Var"] = missed_cleavages_var[raw_file]

    # 3. Charge
    charge = dict()
    for raw_file, group in evidence_data.loc[
        ~evidence_data["is_transferred"], ["charge", "raw file"]
        ].groupby("raw file"):
        charge[raw_file] = group["charge"].value_counts()[2] / len(group)
    charge_median = np.median(list(charge.values()))
    heatmap_charge = dict(
        zip(
            charge.keys(),
            list(map(lambda v: 1 - np.abs(v - charge_median), charge.values())),
        )
    )
    for raw_file in heatmap_charge.keys():
        heatmap_dict[raw_file]["Charge"] = heatmap_charge[raw_file]

    # 7. MS2 OverSampling
    for raw_file, value in oversampling.items():
        heatmap_dict[raw_file]["MS2 OverSampling"] = np.minimum(1.0, (value["1"] / 100))

    # Sort the xnames
    heatmap_xname_order = [
        "Contaminants",
        "Peptide Intensity",
        "Charge",
        "Missed Cleavages",
        "Missed Cleavages Var",
        "ID rate over RT",
        "MS2 OverSampling",
        "Pep Missing Values",
    ]

    for raw_file in heatmap_charge.keys():
        heatmap_dict[raw_file] = {
            key: heatmap_dict[raw_file][key]
            for key in heatmap_xname_order
            if key in heatmap_dict[raw_file].keys()
        }

    return heatmap_dict


# 3-1. evidence.txt: Top Contaminants per Raw file
def evidence_top_contaminants(evidence_df, top_n):
    if any(
        column not in evidence_df.columns
        for column in ["potential contaminant", "proteins", "intensity", "raw file"]
    ):
        return None

    if evidence_df[evidence_df["potential contaminant"] == "+"].empty:
        logger.info("The evidence.txt file does not contain any contaminants")
        return None
    
    evidence_data = evidence_df.copy()

    if "protein names" in evidence_data.columns:
        evidence_data["protein_name"] = evidence_data["protein names"].combine_first(
            evidence_data["proteins"]
        )
    else:
        evidence_data["protein_name"] = evidence_data["proteins"]

    sum_intensity = evidence_data["intensity"].sum()

    contaminant_df = evidence_data[evidence_data["potential contaminant"] == "+"]

    contaminant_count = (
        contaminant_df.groupby("protein_name")["intensity"].sum() / sum_intensity * 100
    )
    contaminant_count = contaminant_count.sort_values(ascending=False)
    top_contaminant = list(contaminant_count.head(top_n).index)

    contaminant_df.loc[~contaminant_df["protein_name"].isin(top_contaminant), "protein_name"] = (
        "Other"
    )

    intensity_per_file = evidence_data.groupby("raw file", as_index=False)["intensity"].sum()
    intensity_per_file.rename(columns={"intensity": "total_intensity"}, inplace=True)

    intensity_per_file_protein = contaminant_df.groupby(
        ["raw file", "protein_name"], as_index=False
    )["intensity"].sum()
    intensity_per_file_protein.rename(columns={"intensity": "contaminant_intensity"}, inplace=True)

    intensity_per_file_protein = pd.merge(
        intensity_per_file_protein, intensity_per_file, on="raw file"
    )

    intensity_per_file_protein["intensity_percent"] = (
        intensity_per_file_protein["contaminant_intensity"]
        / intensity_per_file_protein["total_intensity"]
        * 100
    )

    top_contaminant_dict = {}

    plot_dict = {}
    for raw_file, group in intensity_per_file_protein.groupby("raw file"):
        plot_dict[raw_file] = (
            group[["protein_name", "intensity_percent"]]
            .set_index("protein_name")["intensity_percent"]
            .to_dict()
        )

    top_contaminant_dict["plot_data"] = plot_dict
    top_contaminant_dict["cats"] = list(intensity_per_file_protein["protein_name"].unique())

    return top_contaminant_dict


# 3-2. evidence.txt: peptide intensity distribution
def evidence_peptide_intensity(evidence_df):
    if any(column not in evidence_df.columns for column in ["intensity", "raw file"]):
        return None

    # Take the logarithm and remove zero values
    def box_fun(intensity_df):
        intensity_df["intensity_processed"] = intensity_df["intensity"].apply(
            lambda x: 1 if (pd.isna(x) or x == 0) else x
        )
        intensity_df["log_intensity_processed"] = np.log2(intensity_df["intensity_processed"])
        box_dict = {}
        for raw_file, group in intensity_df.groupby("raw file"):
            log_intensity = group["log_intensity_processed"]
            box_dict[raw_file] = list(log_intensity[log_intensity != 0])
        return box_dict

    if "potential contaminant" in evidence_df.columns:
        contaminant_df = evidence_df[evidence_df["potential contaminant"] == "+"].copy()
        boxplot_data = [box_fun(evidence_df), box_fun(contaminant_df)]
    else:
        boxplot_data = [box_fun(evidence_df), None]

    result_dict = {"box": boxplot_data}

    return result_dict


# 3-3.evidence.txt: charge distribution
def evidence_charge_distribution(evidence_data):
    if any(column not in evidence_data.columns for column in ["charge", "raw file"]):
        return None

    if "potential contaminant" in evidence_data.columns:
        evidence_data = evidence_data[evidence_data["potential contaminant"] != "+"].copy()

    charge_counts = (
        evidence_data.groupby("raw file")["charge"].value_counts().reset_index(name="count")
    )

    plot_dict = {}
    for raw_file, group in charge_counts.groupby("raw file"):
        charge_counts_sorted = group.sort_values(by="charge")
        charge_counts_sorted["charge"] = charge_counts_sorted["charge"].astype(str)
        plot_dict[raw_file] = dict(
            zip(charge_counts_sorted["charge"], charge_counts_sorted["count"])
        )

    charge_dict = {}
    charge_dict["plot_data"] = plot_dict
    charge_dict["cats"] = list(map(str, sorted(charge_counts["charge"].unique())))

    return charge_dict


# 3-4.evidence.txt: Modifications per Raw file
def evidence_modified(evidence_data):
    if any(column not in evidence_data.columns for column in ["modifications", "raw file"]):
        return None

    if "potential contaminant" in evidence_data.columns:
        evidence_data = evidence_data[evidence_data["potential contaminant"] != "+"].copy()

    plot_dict = {}
    modified_cats = []

    for raw_file, group in evidence_data.groupby("raw file"):
        group_processed = mod_group_percentage(group)
        plot_dict[raw_file] = dict(
            zip(group_processed["modifications"], group_processed["percentage"])
        )
        modified_cats.extend(group_processed["modifications"])

    modified_dict = {}
    modified_dict["plot_data"] = plot_dict
    modified_dict["cats"] = list(sorted(modified_cats, key=lambda x: (x == "Modified (Total)", x)))

    return modified_dict


# 3-5.evidence.txt: IDs over RT
def evidence_rt_count(evidence_data):
    if any(column not in evidence_data.columns for column in ["retention time", "raw file"]):
        return None

    if "potential contaminant" in evidence_data.columns:
        evidence_data = evidence_data[evidence_data["potential contaminant"] != "+"].copy()

    rt_range = [evidence_data["retention time"].min(), evidence_data["retention time"].max()]

    def hist_compute(rt_list, rt_range):

        if rt_range[0] < 5:
            rt_range_min = 0
        else:
            rt_range_min = rt_range[0] - 5
            
        counts, bin_edges = np.histogram(
            rt_list, bins=np.arange(rt_range_min, rt_range[1] + 5, 1)
        )
        bin_mid = (bin_edges[:-1] + bin_edges[1:]) / 2
        rt_counts = pd.DataFrame({"retention_time": bin_mid, "counts": counts})

        return dict(zip(rt_counts["retention_time"], rt_counts["counts"]))

    rt_count_dict = {}
    for raw_file, group in evidence_data.groupby("raw file"):
        rt_count_dict[raw_file] = hist_compute(group["retention time"], rt_range)

    return rt_count_dict


# 3-6.evidence.txt: Peak width over RT
def evidence_peak_width_rt(evidence_data):
    if any(
        column not in evidence_data.columns
        for column in ["retention length", "retention time", "raw file"]
    ):
        return None

    if "potential contaminant" in evidence_data.columns:
        evidence_data = evidence_data[evidence_data["potential contaminant"] != "+"].copy()

    evidence_data = evidence_data[["retention length", "retention time", "raw file"]].copy()

    rt_range = [evidence_data["retention time"].min(), evidence_data["retention time"].max()]

    if rt_range[0] < 5:
        rt_range_min = 0
    else:
        rt_range_min = rt_range[0] - 5
    breaks = np.arange(rt_range_min, rt_range[1] + 5, 1)

    def rt_rl_compute(group_df):
        group_df["bin"] = np.digitize(group_df["retention time"], breaks, right=False)
        bin_group = group_df.groupby("bin")
        rt_rl_df = bin_group["retention length"].median().reset_index()
        rt_rl_df.rename(columns={"retention length": "median_RL"}, inplace=True)
        rt_rl_df["bin_RT"] = breaks[rt_rl_df["bin"] - 1]

        return rt_rl_df

    peak_width_rt_dict = {}
    for raw_file, group in evidence_data.groupby("raw file"):
        peak_width_rt_dict[raw_file] = dict(
            zip(rt_rl_compute(group)["bin_RT"], rt_rl_compute(group)["median_RL"])
        )

    return peak_width_rt_dict


# 3-7.evidence.txt: Oversampling (MS/MS counts per 3D-peak)
def evidence_oversampling(evidence_data):
    if any(column not in evidence_data.columns for column in ["ms/ms count", "raw file"]):
        return None

    if "potential contaminant" in evidence_data.columns:
        evidence_data = evidence_data[evidence_data["potential contaminant"] != "+"].copy()

    evidence_data = evidence_data[["ms/ms count", "raw file"]].copy()

    evidence_data["ms/ms count"] = evidence_data["ms/ms count"].apply(
        lambda x: ">=3" if x >= 3 else x
    )
    oversampling_df = evidence_data.groupby("raw file")["ms/ms count"].value_counts().reset_index()
    oversampling_df["ms/ms count"] = oversampling_df["ms/ms count"].astype(str)

    plot_dict = {}
    for raw_file, group in oversampling_df.groupby("raw file"):
        group["fraction"] = group["count"] / group["count"].sum() * 100
        plot_dict[raw_file] = dict(zip(group["ms/ms count"], group["fraction"]))

    oversampling = {}
    oversampling["plot_data"] = plot_dict
    oversampling["cats"] = list(oversampling_df["ms/ms count"].unique())

    return oversampling


# 3-8.evidence.txt: Uncalibrated mass error
def evidence_uncalibrated_mass_error(evidence_data):
    if any(
        column not in evidence_data.columns
        for column in ["uncalibrated mass error [ppm]", "raw file"]
    ):
        return None

    if "potential contaminant" in evidence_data.columns:
        evidence_data = evidence_data[evidence_data["potential contaminant"] != "+"].copy()

    uncalibrated_mass_error = {}
    for raw_file, group in evidence_data.groupby("raw file"):
        mass_error = list(
            group["uncalibrated mass error [ppm]"].map(lambda x: 0 if pd.isna(x) else x)
        )

        uncalibrated_mass_error[raw_file] = [value for value in mass_error if value != 0]

    return uncalibrated_mass_error


# 3-8.evidence.txt: Calibrated mass error
def evidence_calibrated_mass_error(evidence_data):
    if any(column not in evidence_data.columns for column in ["mass error [ppm]", "raw file"]):
        return None

    if "potential contaminant" in evidence_data.columns:
        evidence_data = evidence_data[evidence_data["potential contaminant"] != "+"].copy()

    evidence_data.dropna(subset=["mass error [ppm]"], inplace=True)

    count_bin = evidence_data["mass error [ppm]"].value_counts(sort=False, bins=1000)
    count_bin_data = dict()
    for index in count_bin.index:
        count_bin_data[float(index.mid)] = int(count_bin[index])

    frequency_bin = evidence_data["mass error [ppm]"].value_counts(sort=False, bins=1000, normalize=True)
    frequency_bin_data = dict()
    for index in frequency_bin.index:
        frequency_bin_data[float(index.mid)] = float(frequency_bin[index])

    result_dict = {
        "count": count_bin_data,
        "frequency": frequency_bin_data,
    }

    return result_dict


# 3-9.evidence.txt: Peptide ID count
def evidence_peptide_count(evidence_df, evidence_df_tf):
    if any(
        column not in evidence_df.columns
        for column in ["modified sequence", "is_transferred", "raw file"]
    ):
        return None

    if any(
        column not in evidence_df_tf.columns
        for column in ["modified sequence", "is_transferred", "raw file"]
    ):
        return None

    evidence_data = evidence_df.copy()
    evidence_data_tf = evidence_df_tf.copy()

    if "potential contaminant" in evidence_data.columns:
        evidence_data = evidence_data[evidence_data["potential contaminant"] != "+"].copy()

    if "potential contaminant" in evidence_data_tf.columns:
        evidence_data_tf = evidence_data_tf[
            evidence_data_tf["potential contaminant"] != "+"
        ].copy()

    required_cols = ["raw file", "is_transferred", "modified sequence"]
    evid_df = pd.concat(
        [evidence_data[required_cols], evidence_data_tf[required_cols]], axis=0, ignore_index=True
    )

    def get_peptide_counts(evd_df):

        peptide_counts = pd.DataFrame()
        for raw_file, group in evd_df.groupby("raw file"):
            pep_set_genuine_unique = group[~group["is_transferred"]]["modified sequence"].unique()
            pep_set_all_mb_runique = group[group["is_transferred"]]["modified sequence"].unique()
            pep_count_gen_and_mbr = len(
                set(pep_set_genuine_unique).intersection(pep_set_all_mb_runique)
            )
            pep_count_new_mbr = len(pep_set_all_mb_runique) - pep_count_gen_and_mbr
            pep_count_only_genuine = len(pep_set_genuine_unique) - pep_count_gen_and_mbr

            if any(evd_df["is_transferred"]):
                file_peptide_counts = pd.DataFrame(
                    {
                        "raw file": [raw_file, raw_file, raw_file],
                        "counts": [
                            pep_count_only_genuine,
                            pep_count_gen_and_mbr,
                            pep_count_new_mbr,
                        ],
                        "category": [
                            "Genuine (Exclusive)",
                            "Genuine + Transferred",
                            "Transferred (Exclusive)",
                        ],
                        "MBRgain": [
                            None,
                            None,
                            pep_count_new_mbr
                            / (pep_count_only_genuine + pep_count_gen_and_mbr)
                            * 100,
                        ],
                    }
                )
                categories = [
                    "Genuine (Exclusive)",
                    "Genuine + Transferred",
                    "Transferred (Exclusive)",
                ]
            else:
                file_peptide_counts = pd.DataFrame(
                    {
                        "raw file": [raw_file],
                        "counts": [pep_count_only_genuine],
                        "category": ["Genuine"],
                        "MBRgain": [None],
                    }
                )
                categories = ["Genuine"]

            peptide_counts = pd.concat(
                [peptide_counts, file_peptide_counts], axis=0, ignore_index=True
            )
        return peptide_counts, categories

    peptide_counts_df, cats = get_peptide_counts(evid_df)

    plot_data = {}
    for raw_file, group in peptide_counts_df.groupby("raw file"):
        plot_data[raw_file] = dict(zip(group["category"], group["counts"]))

    peptide_id_count = {}
    peptide_id_count["plot_data"] = plot_data
    peptide_id_count["cats"] = cats
    peptide_id_count["title_value"] = (
        "MBR gain: +{}%".format(round(peptide_counts_df["MBRgain"].mean(), 2))
        if any(evid_df["is_transferred"])
        else ""
    )

    return peptide_id_count


# 3-10.evidence.txt: ProteinGroups count
def evidence_protein_count(evidence_df, evidence_df_tf):
    if any(
        column not in evidence_df.columns
        for column in ["protein group ids", "is_transferred", "raw file"]
    ):
        return None

    if any(
        column not in evidence_df_tf.columns
        for column in ["protein group ids", "is_transferred", "raw file"]
    ):
        return None

    evidence_data = evidence_df.copy()
    evidence_data_tf = evidence_df_tf.copy()

    if "potential contaminant" in evidence_data.columns:
        evidence_data = evidence_data[evidence_data["potential contaminant"] != "+"].copy()
    if "potential contaminant" in evidence_data_tf.columns:
        evidence_data_tf = evidence_data_tf[
            evidence_data_tf["potential contaminant"] != "+"
        ].copy()

    required_cols = ["raw file", "is_transferred", "protein group ids"]
    evid_df = pd.concat(
        [evidence_data[required_cols], evidence_data_tf[required_cols]], axis=0, ignore_index=True
    )

    def get_protein_group_counts(evd_df):

        protein_group_counts = pd.DataFrame()
        for raw_file, group in evd_df.groupby("raw file"):

            group["protein_group_mtd"] = (
                group["protein group ids"] + "_" + group["is_transferred"].astype(str)
            )
            duplicated_protein_group = group[~group["protein_group_mtd"].duplicated()]

            protein_groups = duplicated_protein_group["protein group ids"].apply(
                lambda x: x.split(";")
            )
            protein_group_genuine_unique = (
                protein_groups[~duplicated_protein_group["is_transferred"]].explode().unique()
            )
            protein_group_mbr_unique = (
                protein_groups[duplicated_protein_group["is_transferred"]].explode().unique()
            )
            protein_group_gen_and_mbr = len(
                set(protein_group_genuine_unique).intersection(protein_group_mbr_unique)
            )
            protein_group_mbr_only = len(protein_group_mbr_unique) - protein_group_gen_and_mbr
            protein_group_genuine_only = (
                len(protein_group_genuine_unique) - protein_group_gen_and_mbr
            )

            if any(evd_df["is_transferred"]):
                file_protein_group_counts = pd.DataFrame(
                    {
                        "raw file": [raw_file, raw_file, raw_file],
                        "counts": [
                            protein_group_genuine_only,
                            protein_group_gen_and_mbr,
                            protein_group_mbr_only,
                        ],
                        "category": [
                            "Genuine (Exclusive)",
                            "Genuine + Transferred",
                            "Transferred (Exclusive)",
                        ],
                        "MBRgain": [
                            None,
                            None,
                            protein_group_mbr_only
                            / (protein_group_genuine_only + protein_group_gen_and_mbr)
                            * 100,
                        ],
                    }
                )
                categories = [
                    "Genuine (Exclusive)",
                    "Genuine + Transferred",
                    "Transferred (Exclusive)",
                ]
            else:
                file_protein_group_counts = pd.DataFrame(
                    {
                        "raw file": [raw_file],
                        "counts": [protein_group_genuine_only],
                        "category": ["Genuine"],
                        "MBRgain": [None],
                    }
                )
                categories = ["Genuine"]

            protein_group_counts = pd.concat(
                [protein_group_counts, file_protein_group_counts], axis=0, ignore_index=True
            )
        return protein_group_counts, categories

    protein_group_counts_df, cats = get_protein_group_counts(evid_df)

    plot_data = {}
    for raw_file, group in protein_group_counts_df.groupby("raw file"):
        plot_data[raw_file] = dict(zip(group["category"], group["counts"]))

    protein_group_count = {}
    protein_group_count["plot_data"] = plot_data
    protein_group_count["cats"] = cats
    protein_group_count["title_value"] = (
        "MBR gain: +{}%".format(round(protein_group_counts_df["MBRgain"].mean(), 2))
        if any(evid_df["is_transferred"])
        else ""
    )

    return protein_group_count


# 3-11.evidence.txt: Summary Table
def evidence_summary_table(evid_df):

    if "potential contaminant" in evid_df.columns:
        evid_df = evid_df[evid_df["potential contaminant"] != "+"].copy()

    unique_peptides = (
        evid_df["modified sequence"].nunique() if "modified sequence" in evid_df.columns else None
    )
    msms_count = evid_df["ms/ms count"].sum() if "ms/ms count" in evid_df.columns else None

    if unique_peptides is not None and msms_count is not None:
        
        return {
            "summary_identified_msms_count": int(msms_count),
            "summary_identified_peptides": int(unique_peptides),
        }

    return None


# 3-12.evidence.txt: Delta Mass [Da]
def evidence_delta_mass_da(evidence_data):

    if any(column not in evidence_data.columns for column in ["mass error [da]", "raw file"]):
        return None

    if "potential contaminant" in evidence_data.columns:
        evidence_data = evidence_data[evidence_data["potential contaminant"] != "+"].copy()

    evidence_data.dropna(subset=["mass error [da]"], inplace=True)

    count_bin = evidence_data["mass error [da]"].value_counts(sort=False, bins=1000)
    count_bin_data = dict()
    for index in count_bin.index:
        count_bin_data[float(index.mid)] = int(count_bin[index])

    frequency_bin = evidence_data["mass error [da]"].value_counts(sort=False, bins=1000, normalize=True)
    frequency_bin_data = dict()
    for index in frequency_bin.index:
        frequency_bin_data[float(index.mid)] = float(frequency_bin[index])

    delta_mass_da = {
        "count": count_bin_data,
        "frequency": frequency_bin_data,
    }

    return delta_mass_da


# 4.msms.txt
def get_msms(file_path: Union[Path, str], evidence_df: pd.DataFrame = None):
    """
    Process msms.txt file and extract various metrics.

    Args:
        file_path: Path to the msms.txt file
        evidence_df: Evidence DataFrame for cross-referencing

    Returns:
        Dictionary containing metrics extracted from the msms file
    """
    mq_data = read(file_path)

    # Count MS/MS spectra and peptides
    total_entries = len(mq_data)
    unique_peptides = mq_data["sequence"].nunique() if "sequence" in mq_data.columns else 0
    unique_proteins = mq_data["proteins"].nunique() if "proteins" in mq_data.columns else 0

    logger.info(
        f"Processing msms data with {total_entries} MS/MS spectra, {unique_peptides} unique peptides, and {unique_proteins} unique proteins"
    )

    if evidence_df is None:
        logger.warning("No evidence data provided, skipping missed cleavages calculation")
        return {"missed_cleavages": None}

    # Missed cleavages per Raw file
    logger.debug("Calculating missed cleavages")
    missed_cleavages = msms_missed_cleavages(mq_data, evidence_df)

    # MaxQuant: Search Engine Scores
    search_engine_scores_dict = search_engine_scores(mq_data)

    result = {
        # 'mq_data': mq_data,
        "missed_cleavages": missed_cleavages,
        "search_engine_scores": search_engine_scores_dict,
    }

    logger.info("Completed processing msms data")
    return result


# 4-1.msms.txt: Missed cleavages per Raw file
def msms_missed_cleavages(msms_df: pd.DataFrame, evidence_df: pd.DataFrame):
    if any(
        column not in msms_df.columns for column in ["evidence id", "missed cleavages", "raw file"]
    ):
        return None

    if any(column not in evidence_df.columns for column in ["id"]):
        return None

    # excludes contaminants & excludes Type == 'MULTI-MATCH'
    if "potential contaminant" in evidence_df.columns:
        evidence_not_contaminant_ids = evidence_df[evidence_df["potential contaminant"] != "+"][
            "id"
        ].copy()
    else:
        evidence_not_contaminant_ids = evidence_df["id"].copy()

    not_contaminant_index = msms_df["evidence id"][
        msms_df["evidence id"].isin(evidence_not_contaminant_ids)
    ].index

    msms_not_contaminant = msms_df.loc[
        not_contaminant_index, ["raw file", "missed cleavages", "evidence id"]
    ]
    msms_not_contaminant["missed cleavages"] = msms_not_contaminant["missed cleavages"].astype(
        "str"
    )

    plot_dict = {}
    for raw_file, group in msms_not_contaminant.groupby("raw file"):
        missed_cleavages_df = pd.DataFrame(group["missed cleavages"].value_counts().reset_index())
        missed_cleavages_df["percentage"] = (
            missed_cleavages_df["count"] / missed_cleavages_df["count"].sum() * 100
        )

        plot_dict[raw_file] = dict(
            zip(missed_cleavages_df["missed cleavages"], missed_cleavages_df["percentage"])
        )

    missed_cleavages_dict = {}
    missed_cleavages_dict["plot_data"] = plot_dict
    missed_cleavages_dict["cats"] = sorted(list(msms_not_contaminant["missed cleavages"].unique()))

    return missed_cleavages_dict


# 4-2.msms.txt: Search Engine Scores
def search_engine_scores(msms_df):

    if "score" not in msms_df.columns:
        return None

    bins_start = 0
    bins_end = 300
    bins_step = 6

    bins = list(range(bins_start, bins_end + 1, bins_step)) + [float("inf")]
    labels = [f"{i} ~ {i+bins_step}" for i in range(bins_start, bins_end, bins_step)] + [f"{bins_end} ~ inf"]

    plot_data = list()
    data_labels = list()
    for name, group in msms_df.groupby("raw file"):

        group["score_bin"] = pd.cut(group["score"], bins=bins, labels=labels, right=False)
        score_dist = group["score_bin"].value_counts().sort_index().reset_index()

        plot_data.append(
            {
                k: {"count": v} for k, v in zip(score_dist["score_bin"], score_dist["count"])
            }
        )
        data_labels.append(
            {
                "name": name,
                "ylab": "Counts"
            }
        )

    result = {
        "plot_data": plot_data,
        "data_labels": data_labels,
    }

    return result


# 5.msScans.txt
def get_msms_scans(file_path: Union[Path, str]):
    """
    Process msmsScans.txt or msScans.txt file and extract various metrics.

    Args:
        file_path: Path to the msmsScans.txt or msScans.txt file

    Returns:
        Dictionary containing metrics extracted from the msmsScans file
    """
    mq_data = read(file_path)

    # Count MS scans and get retention time range
    total_scans = len(mq_data)
    rt_min = mq_data["retention time"].min() if "retention time" in mq_data.columns else 0
    rt_max = mq_data["retention time"].max() if "retention time" in mq_data.columns else 0
    raw_files = mq_data["raw file"].nunique() if "raw file" in mq_data.columns else 0

    logger.info(
        f"Processing msmsScans data with {total_scans} scans across {raw_files} raw files (RT range: {rt_min:.2f}-{rt_max:.2f} min)"
    )

    # TODO check 'Scan event number'

    logger.debug("Rounding retention time values")
    mq_data["round_RT"] = mq_data["retention time"].apply(lambda x: round(x / 2) * 2)

    # Ion Injection Time over RT
    logger.debug("Calculating ion injection time over retention time")
    ion_injec_time_rt = msms_scans_ion_injec_time_rt(mq_data)

    # TopN over RT
    logger.debug("Calculating TopN over retention time")
    top_over_rt = msms_scans_top_over_rt(mq_data)

    # TopN
    logger.debug("Calculating TopN")
    top_n = msms_scans_top_n(mq_data)

    result = {
        "ion_injec_time_rt": ion_injec_time_rt,
        "top_n": top_n,
        "top_over_rt": top_over_rt,
        "summary_msms_spectra": total_scans,
    }

    logger.info("Completed processing msmsScans data")
    return result


# 5-1.msmsScans.txt: Ion Injection Time over RT
def msms_scans_ion_injec_time_rt(msms_scans_df):
    if any(column not in msms_scans_df.columns for column in ["raw file", "ion injection time"]):
        return None

    if msms_scans_df["ion injection time"].isna().all():
        return None

    median_ion_injec_time_df = (
        msms_scans_df.groupby(["raw file", "round_RT"])["ion injection time"]
        .median()
        .reset_index()
    )
    median_ion_injec_time_df = median_ion_injec_time_df.rename(
        columns={"ion injection time": "median_ion_injection_time"}
    )

    mean_ion_injec_time_df = (
        msms_scans_df.groupby(["raw file"])["ion injection time"].mean().reset_index()
    )
    mean_ion_injec_time_df = mean_ion_injec_time_df.rename(
        columns={"ion injection time": "mean_ion_injection_time"}
    )
    mean_ion_injec_time_df["int_mean_ion_injection_time"] = mean_ion_injec_time_df[
        "mean_ion_injection_time"
    ].apply(lambda x: int(x) if not pd.isna(x) else 0)
    mean_ion_injec_time_df["int_mean_ion_injection_time"] = mean_ion_injec_time_df[
        "int_mean_ion_injection_time"
    ].astype(str)
    mean_ion_injec_time_df["raw_file_mean_ion_time"] = (
        mean_ion_injec_time_df["raw file"]
        + " (~"
        + mean_ion_injec_time_df["int_mean_ion_injection_time"].astype(str)
        + "ms)"
    )

    result_df = pd.merge(
        median_ion_injec_time_df,
        mean_ion_injec_time_df[["raw file", "raw_file_mean_ion_time"]],
        on="raw file",
    )

    ion_injec_time_dict = {}
    for raw_file, group in result_df.groupby("raw_file_mean_ion_time"):
        ion_injec_time_dict[raw_file] = dict(
            zip(group["round_RT"], group["median_ion_injection_time"])
        )

    return ion_injec_time_dict


# 5-2.msmsScans.txt: TopN over RT
def msms_scans_top_over_rt(msms_scans_df):
    if any(
        column not in msms_scans_df.columns
        for column in ["raw file", "retention time", "scan event number"]
    ):
        return None

    # Scan event number:
    #   This number indicates which MS/MS scan this one is in the consecutive order of the MS/MS scans that are acquired after an MS scan.
    # Retention time:
    #   Time point along the elution profile at which the MS/MS data was recorded.
    msms_scans_data = msms_scans_df[
        ["raw file", "retention time", "round_RT", "scan event number"]
    ].copy()
    msms_scans_data = msms_scans_data.sort_values(by=["raw file", "retention time"])

    def find_local_maxima(list_data):
        local_max_indices = np.zeros(len(list_data), dtype=bool)
        for i in range(0, len(list_data) - 1):
            if np.isnan(list_data[i]):
                continue
            if list_data[i] > list_data[i + 1] and list_data[i] >= 0:
                local_max_indices[i] = True
        return local_max_indices

    local_max_se_number_df = pd.DataFrame()
    for raw_file, group in msms_scans_data.groupby("raw file"):
        local_max_se_number_df = pd.concat(
            [local_max_se_number_df, group[find_local_maxima(list(group["scan event number"]))]],
            axis=0,
            ignore_index=True,
        )

    local_max_se_number_df = local_max_se_number_df.rename(
        columns={"scan event number": "local_max_scan_event_number"}
    )
    median_se_number_df = (
        local_max_se_number_df.groupby(["raw file", "round_RT"])["local_max_scan_event_number"]
        .median()
        .reset_index()
    )

    scan_event_number_dict = {}
    for raw_file, group in median_se_number_df.groupby("raw file"):
        scan_event_number_dict[raw_file] = dict(
            zip(group["round_RT"], group["local_max_scan_event_number"])
        )

    return scan_event_number_dict


# 5-3.msmsScans.txt: TopN
def msms_scans_top_n(msms_scans_df):
    if any(column not in msms_scans_df.columns for column in ["raw file", "scan event number"]):
        return None

    msms_scans_data = msms_scans_df[["scan event number", "raw file"]].copy()

    while True:
        scan_event_index = 1 + np.where(np.diff(msms_scans_data["scan event number"]) > 1)[0]
        if len(scan_event_index) == 0:
            break
        msms_scans_data.loc[scan_event_index, "scan event number"] -= 1

    file_se_count_df = (
        msms_scans_data[["scan event number", "raw file"]].value_counts().reset_index(name="count")
    )
    max_se_number = max(file_se_count_df["scan event number"])

    def process_group(group_df, max_scan_event_number, raw_file):
        event_count = group_df["count"].values

        if not all(event_count[i] >= event_count[i + 1] for i in range(len(event_count) - 1)):
            raise ValueError("Scan event distribution is not monotonically increasing!")
        if max(group_df["scan event number"]) != len(group_df):
            raise ValueError("Scan event distribution has unexpected holes...!")

        event_pre = np.append(event_count[1:], 0)
        event_diff = event_count - event_pre

        se_number = group_df["scan event number"].values
        if max(se_number) < max_scan_event_number:
            event_diff = list(event_diff) + [0] * (max_scan_event_number - max(se_number))
            se_number = list(se_number) + list(
                range(max(se_number) + 1, max_scan_event_number + 1)
            )

        result_df = pd.DataFrame(
            {"raw file": raw_file, "scan event number": se_number, "count": event_diff}
        )
        return result_df

    file_se_count_ratio = pd.DataFrame()
    for raw_file, group in file_se_count_df.groupby("raw file"):
        file_se_count_ratio = pd.concat(
            [file_se_count_ratio, process_group(group, max_se_number, raw_file)],
            axis=0,
            ignore_index=True,
        )

    file_se_count_ratio["scan event number"] = file_se_count_ratio["scan event number"].astype(str)

    plot_dict = {}
    for raw_file, group in file_se_count_ratio.groupby("raw file"):
        plot_dict[raw_file] = dict(zip(group["scan event number"], group["count"]))

    se_count_dict = {}
    se_count_dict["plot_data"] = plot_dict
    se_count_dict["cats"] = [
        str(x) for x in reversed(file_se_count_ratio["scan event number"].unique())
    ]

    return se_count_dict


# 6.parameters.txt
def get_parameters(file_path: Path | ReadCsvBuffer[bytes] | ReadCsvBuffer[str] | str):
    """
    Process parameters.txt file and extract parameters table.

    Args:
        file_path: Path to the parameters.txt file

    Returns:
        Dictionary containing parameters table
    """
    mq_data = read(file_path)
    logger.info(f"Processing parameters data with {len(mq_data)} entries")

    logger.debug("Extracting parameters table")
    parameters_tb_dict = parameters_table(mq_data)

    result = {"parameters_tb_dict": parameters_tb_dict}

    logger.info("Completed processing parameters data")
    return result


# 6-1.parameters.txt: Parameters
def parameters_table(parameters_df):
    """
    Extract parameters table from parameters DataFrame.

    Args:
        parameters_df: DataFrame containing parameters data

    Returns:
        Dictionary containing parameters table
    """
    if any(column not in parameters_df.columns for column in ["parameter", "value"]):
        logger.warning("Required columns not found in parameters table")
        return None

    logger.debug("Filtering AIF parameters")
    parameters_data = parameters_df[~parameters_df["parameter"].str.startswith("AIF")]

    logger.debug("Processing FASTA file paths")
    fasta_files = parameters_data[parameters_data["parameter"] == "Fasta file"]["value"].values[0]
    fasta_files = fasta_files.split(";")
    logger.debug(f"Found {len(fasta_files)} FASTA files")

    def parse_location(location):
        if "\\" in location:
            location = location.replace("\\", "/")
        return os.path.basename(location)

    fasta_file_list = [parse_location(fasta_file) for fasta_file in fasta_files]
    fasta_file_list = ";".join(fasta_file_list)
    logger.debug(f"Processed FASTA file paths: {fasta_file_list}")

    logger.debug("Creating parameters table")
    table_data = parameters_data.drop_duplicates(subset="parameter", keep="first").reset_index(
        drop=True
    )
    table_data.loc[table_data["parameter"] == "Fasta file", "value"] = fasta_file_list

    parameters_dict = {}
    for index, row in table_data.iterrows():
        parameters_dict[index + 1] = row.to_dict()

    logger.debug(f"Created parameters table with {len(parameters_dict)} entries")
    return parameters_dict


def mod_group_percentage(group):

    if "Modifications" in group.columns:
        group.rename(columns={"Modifications": "modifications"}, inplace=True)

    counts = group["modifications"].str.split(",").explode().value_counts()
    percentage_df = (counts / len(group["modifications"]) * 100).reset_index()
    percentage_df.columns = ["modifications", "percentage"]

    # Modified (Total)
    percentage_df.loc[percentage_df["modifications"] == "Unmodified", "percentage"] = (
        100 - percentage_df.loc[percentage_df["modifications"] == "Unmodified", "percentage"]
    )
    percentage_df.loc[percentage_df["modifications"] == "Unmodified", "modifications"] = (
        "Modified (Total)"
    )

    return percentage_df
