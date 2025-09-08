import logging
import re
import pandas as pd
import numpy as np
import os

from ..common.file_utils import drop_empty_row
from statsmodels.nonparametric.smoothers_lowess import lowess


DEFAULT_BINS = 500

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def extract_condition_and_replicate(run_name):

    match = re.search(r"^(.*?)([A-Za-z]*)(\d+)$", run_name)

    if match:
        condition_base = match.group(1) + match.group(2)
        
        if condition_base.endswith("_"):
            condition_base = condition_base[: -1]

        replicate = int(match.group(3))

        return condition_base, replicate
    else:
        log.warning("Failed to identify condition groups in DIA report.tsv!")

def calculate_dia_intensity_std(df):

    df_sub = df.copy()
    df_sub[["run_condition", "run_replicate"]] = df_sub["Run"].apply(
        lambda x: pd.Series(extract_condition_and_replicate(x))
    )

    grouped_std = (
        df_sub
        .groupby(["run_condition", "Modified.Sequence"])["log_intensity"]
        .std()
        .reset_index(name="log_intensity_std")
    )

    plot_data = {
        condition: group["log_intensity_std"].dropna().tolist()
        for condition, group in grouped_std.groupby("run_condition")
    }

    return plot_data

def calculate_msms_count(df):
    count_df = df.groupby(
        [
            "Run",
            "Stripped.Sequence",
            "Precursor.Charge"
        ]
    )["MS2.Scan"].nunique().reset_index(name="msms_count")

    run_counts = count_df.groupby("Run")["msms_count"].value_counts().reset_index(name="msms_count_run")

    run_counts["msms_count_str"] = run_counts["msms_count"].apply(lambda x: ">=3" if x >= 3 else str(x))
    merged_df = run_counts.groupby(["Run", "msms_count_str"])["msms_count_run"].sum().reset_index()

    result_dict = dict()
    for run_name, group in merged_df.groupby("Run"):
        result_dict[str(run_name)] = dict(zip(group["msms_count_str"], group["msms_count_run"]))

    return {
         "plot_data": result_dict,
         "cats": list(merged_df["msms_count_str"].unique())
    }

# DIA-NN: Peptides Quantification Table
def create_peptides_table(report_df, sample_df, file_df):

    # Validation: remove rows with 0 or NA Precursor.Normalised values
    report_data = report_df[report_df["Precursor.Normalised"] > 0].copy()
    report_data = drop_empty_row(report_data, ["Protein.Names", "Stripped.Sequence"])
    
    report_data["BestSearchScore"] = 1 - report_data["Q.Value"]

    table_dict = dict()
    for sequence_protein, group in report_data.groupby(["Stripped.Sequence", "Protein.Names"]):

        table_dict[sequence_protein] = {
            "ProteinName": sequence_protein[1],
            "PeptideSequence": sequence_protein[0],
            "BestSearchScore": group["BestSearchScore"].min(),
            "Average Intensity": np.log10(
                group["Precursor.Normalised"].mean()
            ),
        }

    headers = {
        "ProteinName": {
            "title": "Protein Name",
            "description": "Name/Identifier(s) of the protein (group)",
            "minrange": "200",
        },
        "PeptideSequence": {"title": "Peptide Sequence"},
        "BestSearchScore": {"title": "Best Search Score", "format": "{:,.4f}"},
        "Average Intensity": {
            "title": "Average Intensity",
            "description": "Average intensity across all conditions",
            "format": "{:,.4f}",
        },
    }

    if not sample_df.empty and not file_df.empty:

        sample_cond_df = pd.merge(
            sample_df[["Sample", "MSstats_Condition"]],
            file_df[["Sample", "Spectra_Filepath"]],
            on="Sample"
        )
        sample_cond_df["Run"] = sample_cond_df["Spectra_Filepath"].apply(
            lambda x: os.path.splitext(x)[0]
        )

        cond_report_data = pd.merge(
            report_data[["Stripped.Sequence", "Protein.Names", "Precursor.Normalised", "Run"]],
            sample_cond_df[["Run", "MSstats_Condition"]].drop_duplicates(),
            on="Run"
        )

        for sequence_protein, group in cond_report_data.groupby(["Stripped.Sequence", "Protein.Names"]):
            
            condition_data = dict()
            for condition, sub_group in group.groupby("MSstats_Condition"):
                condition_data[str(condition)] = np.log10(
                    sub_group["Precursor.Normalised"].mean()
                )

            table_dict[sequence_protein].update(condition_data)

        for exp_condition in sample_df["MSstats_Condition"].drop_duplicates():

            headers[str(exp_condition)] = {
                "title": str(exp_condition),
                "description": "MSstats Condition",
                "format": "{:,.4f}"
            }

    result_dict = {i: v for i, (_, v) in enumerate(table_dict.items(), start=1)}

    return result_dict, headers

# DIA-NN: Protein Quantification Table
def create_protein_table(report_df, sample_df, file_df):

    # Validation: remove rows with 0 or NA Precursor.Normalised values
    report_data = report_df[report_df["Precursor.Normalised"] > 0].copy()
    report_data = drop_empty_row(report_data, ["Protein.Names", "Stripped.Sequence"])

    table_dict = dict()
    for protein_name, group in report_data.groupby("Protein.Names"):

        table_dict[protein_name] = {
            "ProteinName": protein_name,
            "Peptides_Number": group["Stripped.Sequence"].nunique(),
            "Average Intensity": np.log10(
                group["Precursor.Normalised"].mean()
            ),
        }

    headers = {
        "ProteinName": {
            "title": "Protein Name",
            "description": "Name/Identifier(s) of the protein (group)",
        },
        "Peptides_Number": {
            "title": "Number of Peptides",
            "description": "Number of peptides per proteins",
            "format": "{:,.0f}",
        },
        "Average Intensity": {
            "title": "Average Intensity",
            "description": "Average intensity across all conditions",
            "format": "{:,.4f}",
        },
    }

    if not sample_df.empty and not file_df.empty:

        sample_cond_df = pd.merge(
            sample_df[["Sample", "MSstats_Condition"]],
            file_df[["Sample", "Spectra_Filepath"]],
            on="Sample"
        )
        sample_cond_df["Run"] = sample_cond_df["Spectra_Filepath"].apply(
            lambda x: os.path.splitext(x)[0]
        )

        cond_report_data = pd.merge(
            report_data[["Stripped.Sequence", "Protein.Names", "Precursor.Normalised", "Run"]],
            sample_cond_df[["Run", "MSstats_Condition"]].drop_duplicates(),
            on="Run"
        )

        for protein_name, group in cond_report_data.groupby("Protein.Names"):
            
            condition_data = dict()
            for condition, sub_group in group.groupby("MSstats_Condition"):
                condition_data[str(condition)] = np.log10(
                    sub_group["Precursor.Normalised"].mean()
                )

            table_dict[protein_name].update(condition_data)

        for exp_condition in sample_df["MSstats_Condition"].drop_duplicates():

            headers[str(exp_condition)] = {
                "title": str(exp_condition),
                "description": "MSstats Condition",
                "format": "{:,.4f}"
            }

    result_dict = {i: v for i, (_, v) in enumerate(table_dict.items(), start=1)}

    return result_dict, headers

def cal_feature_avg_rt(report_data, col):

    # RT: the retention time (RT) of the PSM in minutes
    sub_df = report_data[[col, "Run", "RT"]].copy()

    # RT bin
    sub_df["RT_bin"] = pd.cut(sub_df["RT"], bins=DEFAULT_BINS)
    sub_df['RT_bin_mid'] = sub_df['RT_bin'].apply(lambda x: x.mid)
    result = sub_df.groupby(
        ["Run", "RT_bin_mid"],
        observed=False
    )[col].mean().reset_index()
    result[col] = result[col].fillna(0)

    plot_dict = {
        str(run): group.set_index("RT_bin_mid")[col].to_dict()
        for run, group in result.groupby("Run")
    }

    return plot_dict

# Lowess (Loess)
def cal_rt_irt_loess(report_df, frac=0.3, data_bins: int=DEFAULT_BINS):

    df = report_df.copy()
    
    # bin
    x_min, x_max = df["iRT"].min(), df["iRT"].max()
    bins = np.linspace(x_min, x_max, data_bins)
    
    plot_dict = dict()    
    for run, group in df.groupby("Run"):

        group_sorted = group.sort_values("iRT")
        x = group_sorted["iRT"].values
        y = group_sorted["RT"].values

        # lowess
        smoothed = lowess(y, x, frac=frac)
        smoothed_x = smoothed[:, 0]
        smoothed_y = smoothed[:, 1]

        bin_indices = np.digitize(smoothed_x, bins)
        binned_dict = dict()
        for i in range(1, len(bins)):
            mask = bin_indices == i
            if np.any(mask):
                x_bin_mean = float(smoothed_x[mask].mean())
                y_bin_mean = float(smoothed_y[mask].mean())
                binned_dict[x_bin_mean] = y_bin_mean

        plot_dict[run] = binned_dict
    
    return plot_dict