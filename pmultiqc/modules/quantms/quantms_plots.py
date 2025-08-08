import logging
from multiqc.plots import bargraph, box, table

from ..common.common_plots import (
    draw_ids_rt_count,
    remove_subtitle
)
from ..core.section_groups import add_sub_section
from ..maxquant.maxquant_utils import evidence_rt_count
from ..common.common_plots import draw_oversampling
from . import quantms_utils

import re
import numpy as np
import itertools


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

def draw_dia_intensitys(sub_section, report_df):

    df_sub = report_df[report_df["Precursor.Normalised"] > 0].copy()
    df_sub["log_intensity"] = np.log2(df_sub["Precursor.Normalised"])
    
    draw_dia_intensity_dis(sub_section, df_sub)

    if can_groupby_for_std(report_df, "Run"):
        draw_dia_intensity_std(sub_section, df_sub)

def draw_dia_ms2s(sub_section, df):

    draw_dia_ms2_charge(sub_section, df)

    if "MS2.Scan" in df.columns:
        msms_count_data = quantms_utils.calculate_msms_count(df)
        draw_oversampling(
            sub_section,
            msms_count_data["plot_data"],
            msms_count_data["cats"],
            False
        )

def draw_dia_time_mass(sub_section, df):

    df = df[["Run", "RT"]].copy()
    df.rename(
        columns={
            "Run": "raw file",
            "RT": "retention time"
        },
        inplace=True
    )

    ids_over_rt = evidence_rt_count(df)

    draw_ids_rt_count(
        sub_section,
        ids_over_rt,
        "dia"
    )

# Intensity Distribution
def draw_dia_intensity_dis(sub_section, df):

    box_data = {
        str(run): group["log_intensity"].dropna().tolist()
        for run, group in df.groupby("Run")
    }

    draw_config = {
        "id": "intensity_distribution_box",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Intensity Distribution",
        "tt_decimals": 5,
        "xlab": "log2(Precursor.Normalised)",
    }

    box_html = box.plot(
        list_of_data_by_sample=box_data,
        pconfig=draw_config
    )

    box_html = remove_subtitle(box_html)

    add_sub_section(
        sub_section=sub_section,
        plot=box_html,
        order=3,
        description="log2(Precursor.Normalised) for each Run.",
        helptext="""
            [DIA-NN: report.tsv] log2(Precursor.Normalised) for each Run.
            """
    )

# Charge-state of Per File
def draw_dia_ms2_charge(sub_section, df):

    df = df.copy()
    df["Precursor.Charge"] = df["Precursor.Charge"].astype("str")

    bar_data = (
        df.groupby("Run")["Precursor.Charge"]
        .value_counts()
        .sort_index()
        .unstack(fill_value=0)
        .to_dict(orient="index")
    )

    bar_data = {str(k): v for k, v in bar_data.items()}

    draw_config = {
        "id": "charge_state_of_per_file",
        "cpswitch": True,
        "title": "Charge-state of Per File",
        "tt_decimals": 0,
        "ylab": "Count",
    }

    bar_html = bargraph.plot(
        data=bar_data,
        pconfig=draw_config,
    )

    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=6,
        description="The distribution of the charge-state of the precursor ion.",
        helptext="""
            [DIA-NN: report.tsv] The distribution of the charge-state of the precursor ion (Precursor.Charge).
            """
    )

# DIA: Standard Deviation of Intensity
def can_groupby_for_std(df, col):
    unique_vals = df[col].drop_duplicates()
    
    regex = re.compile(r"^(.*?)([A-Za-z]*)(\d+)$")
    unmatched = [val for val in unique_vals if not regex.match(str(val))]
    
    if len(unmatched) > 0:
        return False
    else:
        return True

def draw_dia_intensity_std(sub_section, df):

    box_data = quantms_utils.calculate_dia_intensity_std(df)

    draw_box_config = {
        "id": "dia_std_intensity_box",
        "title": "Standard Deviation of Intensity",
        "cpswitch": False,
        "tt_decimals": 5,
        "xlab": "Standard Deviation of log2(Precursor.Normalised)",
    }
    
    box_html = box.plot(
        list_of_data_by_sample=box_data,
        pconfig=draw_box_config,
    )

    box_html = remove_subtitle(box_html)

    add_sub_section(
        sub_section=sub_section,
        plot=box_html,
        order=6,
        description="Standard deviation of intensity under different experimental conditions.",
        helptext="""
            [DIA-NN: report.tsv] First, identify the experimental conditions from the "Run" name. 
            Then, group the data by experimental condition and Modified.Sequence, and calculate 
            the standard deviation of log2(Precursor.Normalised).
            """
    )

# DIA-NN Quantification Table
def draw_diann_quant_table(sub_section, diann_report, sample_df, file_df):

    # Peptides Quantification Table
    peptides_table, peptides_headers = quantms_utils.create_peptides_table(
        diann_report,
        sample_df,
        file_df
    )
    draw_peptides_table(sub_section, peptides_table, peptides_headers)

    # Protein Quantification Table
    protein_table, protein_headers = quantms_utils.create_protein_table(
        diann_report,
        sample_df,
        file_df
    )
    draw_protein_table(sub_section, protein_table, protein_headers)

# Draw: Peptides Quantification Table (DIA-NN)
def draw_peptides_table(sub_section, table_data, headers):

    draw_config = {
        "id": "peptides_quantification_table",
        "title": "Peptides Quantification Table",
        "save_file": False,
        "sort_rows": False,
        "only_defined_headers": True,
        "col1_header": "PeptideID",
        "no_violin": True,
    }

    # only use the first 50 lines for the table
    display_rows = 50
    table_html = table.plot(
        dict(itertools.islice(table_data.items(), display_rows)),
        headers=headers,
        pconfig=draw_config,
    )

    add_sub_section(
        sub_section=sub_section,
        plot=table_html,
        order=1,
        description="""
            This plot shows the quantification information of peptides in the final result (DIA-NN report).
            """,
        helptext="""
            The quantification information of peptides is obtained from the DIA-NN output file. 
            The table shows the quantitative level and distribution of peptides in different study variables, 
            run and peptiforms. The distribution show all the intensity values in a bar plot above and below 
            the average intensity for all the fractions, runs and peptiforms.

            * BestSearchScore: It is equal to min(1 - Q.Value) for DIA-NN datasets.
            * Average Intensity: Average intensity of each peptide sequence across all conditions (0 or NA ignored).
            * Peptide intensity in each condition (Eg. `CT=Mixture;CN=UPS1;QY=0.1fmol`).
            """
    )

# Draw: Protein Quantification Table (DIA-NN)
def draw_protein_table(sub_sections, table_data, headers):

    draw_config = {
            "id": "protein_quant_result",
            "title": "Protein Quantification Table",
            "save_file": False,
            "sort_rows": False,
            "only_defined_headers": True,
            "col1_header": "ProteinID",
            "no_violin": True,
        }
    
    display_rows = 50
    table_html = table.plot(
        dict(itertools.islice(table_data.items(), display_rows)),
        headers=headers,
        pconfig=draw_config
    )

    add_sub_section(
        sub_section=sub_sections,
        plot=table_html,
        order=2,
        description="""
            This plot shows the quantification information of proteins in the final result (DIA-NN report).
            """,
        helptext="""
            The quantification information of proteins is obtained from the DIA-NN output file.
            The table shows the quantitative level and distribution of proteins in different study variables and run.

            * Peptides_Number: The number of peptides for each protein.
            * Average Intensity: Average intensity of each protein across all conditions (0 or NA ignored).
            * Protein intensity in each condition (Eg. `CT=Mixture;CN=UPS1;QY=0.1fmol`): Summarize intensity of peptides.
            """
    )
