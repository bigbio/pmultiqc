import logging
import re
import pandas as pd
import numpy as np
import os

from ..common.file_utils import drop_empty_row

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

def get_best_intensity_column(df):
    """Determine the best intensity column to use, preferring non-normalized values"""
    if "Precursor.Quantity" in df.columns:
        return "Precursor.Quantity"
    elif "Ms1.Area" in df.columns:
        return "Ms1.Area"
    elif "Precursor.Normalised" in df.columns:
        return "Precursor.Normalised"
    else:
        return None

# DIA-NN: Peptides Quantification Table
def create_peptides_table(report_df, sample_df, file_df):

    # Get the best intensity column available
    intensity_col = get_best_intensity_column(report_df)
    if intensity_col is None:
        return {}, {}  # No intensity column available
    
    # Validation: remove rows with 0 or NA intensity values
    report_data = report_df[report_df[intensity_col] > 0].copy()
    report_data = drop_empty_row(report_data, ["Protein.Names", "Stripped.Sequence"])
    
    report_data["BestSearchScore"] = 1 - report_data["Q.Value"]

    table_dict = dict()
    for sequence_protein, group in report_data.groupby(["Stripped.Sequence", "Protein.Names"]):

        table_dict[sequence_protein] = {
            "ProteinName": sequence_protein[1],
            "PeptideSequence": sequence_protein[0],
            "BestSearchScore": group["BestSearchScore"].min(),
            "Average Intensity": np.log10(
                group[intensity_col].mean()
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
            report_data[["Stripped.Sequence", "Protein.Names", intensity_col, "Run"]],
            sample_cond_df[["Run", "MSstats_Condition"]].drop_duplicates(),
            on="Run"
        )

        for sequence_protein, group in cond_report_data.groupby(["Stripped.Sequence", "Protein.Names"]):
            
            condition_data = dict()
            for condition, sub_group in group.groupby("MSstats_Condition"):
                condition_data[str(condition)] = np.log10(
                    sub_group[intensity_col].mean()
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

    # Get the best intensity column available
    intensity_col = get_best_intensity_column(report_df)
    if intensity_col is None:
        return {}, {}  # No intensity column available

    # Validation: remove rows with 0 or NA intensity values
    report_data = report_df[report_df[intensity_col] > 0].copy()
    report_data = drop_empty_row(report_data, ["Protein.Names", "Stripped.Sequence"])

    table_dict = dict()
    for protein_name, group in report_data.groupby("Protein.Names"):

        table_dict[protein_name] = {
            "ProteinName": protein_name,
            "Peptides_Number": group["Stripped.Sequence"].nunique(),
            "Average Intensity": np.log10(
                group[intensity_col].mean()
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
            report_data[["Stripped.Sequence", "Protein.Names", intensity_col, "Run"]],
            sample_cond_df[["Run", "MSstats_Condition"]].drop_duplicates(),
            on="Run"
        )

        for protein_name, group in cond_report_data.groupby("Protein.Names"):
            
            condition_data = dict()
            for condition, sub_group in group.groupby("MSstats_Condition"):
                condition_data[str(condition)] = np.log10(
                    sub_group[intensity_col].mean()
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
