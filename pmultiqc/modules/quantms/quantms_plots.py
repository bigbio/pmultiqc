import logging
from multiqc.plots import bargraph, box, table, linegraph

from ..common.common_plots import (
    draw_ids_rt_count,
    remove_subtitle
)
from ..core.section_groups import add_sub_section
from ..maxquant.maxquant_utils import evidence_rt_count
from . import quantms_utils, mzidentml_utils
from ..common.calc_utils import cal_delta_mass_dict

import re
import numpy as np
import itertools


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

def draw_dia_intensitys(sub_section, report_df):

    df_sub = report_df[report_df["Precursor.Quantity"] > 0].copy()
    df_sub["log_intensity"] = np.log2(df_sub["Precursor.Quantity"])
    
    draw_dia_intensity_dis(sub_section, df_sub)

    if can_groupby_for_std(report_df, "Run"):
        draw_dia_intensity_std(sub_section, df_sub)


def draw_dia_ms1(sub_section, df):

    # Ms1.Area: non-normalised MS1 peak area
    if "Ms1.Area" in df.columns:
        df_sub = df[df["Ms1.Area"] > 0][["Ms1.Area", "Run"]].copy()
        if len(df_sub) > 0:
            df_sub["log_ms1_area"] = np.log2(df_sub["Ms1.Area"])
            draw_dia_ms1_area(sub_section, df_sub)

def draw_dia_ms2s(sub_section, df):

    # Distribution of Precursor Charges
    if "Precursor.Charge" in df.columns:
        draw_dia_whole_exp_charge(sub_section, df)

    # Charge-state of Per File
    draw_dia_ms2_charge(sub_section, df)

    # MS2.Scan does not exist in the main report of DIA-NN versions 2.0 and later.

def draw_dia_time_mass(sub_section, df):

    # IDs over RT
    rt_df = df[["Run", "RT"]].copy()
    rt_df.rename(
        columns={
            "Run": "raw file",
            "RT": "retention time"
        },
        inplace=True
    )
    ids_over_rt = evidence_rt_count(rt_df)
    draw_ids_rt_count(
        sub_section,
        ids_over_rt,
        "dia"
    )

    # Ms1.Apex.Mz.Delta: difference between observed precursor m/z and the theoretical value
    if "Ms1.Apex.Mz.Delta" in df.columns:
        draw_dia_delta_mass(sub_section, df)

# DIA-NN: RT Quality Control
def draw_dia_rt_qc(sub_section, report_df):

    df = report_df.copy()

    # 1. Normalisation Factor over RT
    #   Normalisation.Factor: normalisation factor applied to the precursor in the specific run, 
    #   i.e. normalised quantity = normalisation factor X non-normalised quantity
    required_cols = ["Precursor.Normalised", "Precursor.Quantity"]
    if "Normalisation.Factor" not in df.columns and all(
        col in df.columns for col in required_cols
    ):
        df["Normalisation.Factor"] = df[required_cols[0]] / df[required_cols[1]]

    if "Normalisation.Factor" in df.columns:
        norm_factor_rt = quantms_utils.cal_feature_avg_rt(df, "Normalisation.Factor")
        draw_norm_factor_rt(sub_section, norm_factor_rt)

    # 2. FWHM over RT
    #   FWHM: estimated peak width at half-maximum
    if "FWHM" in df.columns:
        fwhm_rt = quantms_utils.cal_feature_avg_rt(df, "FWHM")
        draw_fwhm_rt(sub_section, fwhm_rt)
    
    # 3. Peak Width over RT
    #   RT.Start and RT.Stop peak boundaries
    if all(col in df.columns for col in ["RT.Start", "RT.Stop"]):
        df["peak_width"] = df["RT.Stop"] - df["RT.Start"]
        peak_width_rt = quantms_utils.cal_feature_avg_rt(df, "peak_width")
        draw_peak_width_rt(sub_section, peak_width_rt)

    # 4. Absolute RT Error over RT
    if all(col in df.columns for col in ["RT", "Predicted.RT"]):
        df["rt_error"] = abs(df["RT"] - df["Predicted.RT"])
        rt_error_rt = quantms_utils.cal_feature_avg_rt(df, "rt_error")
        draw_rt_error_rt(sub_section, rt_error_rt)

    # 5. loess(RT ~ iRT)
    if all(col in df.columns for col in ["RT", "iRT"]):
        rt_irt_loess = quantms_utils.cal_rt_irt_loess(df)
        draw_loess_rt_irt(sub_section, rt_irt_loess)


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
        "xlab": "log2(Precursor.Quantity)",
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
        description="log2(Precursor.Quantity) for each Run.",
        helptext="""
            [DIA-NN: report.tsv] log2(Precursor.Quantity) for each Run.
            """
    )

# Ms1.Area non-normalised MS1 peak area
def draw_dia_ms1_area(sub_section, df):

    box_data = {
        str(run): group["log_ms1_area"].dropna().tolist()
        for run, group in df.groupby("Run")
    }

    draw_config = {
        "id": "ms1_area_distribution_box",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Ms1 Area Distribution",
        "tt_decimals": 5,
        "xlab": "log2(Ms1.Area)",
    }

    box_html = box.plot(
        list_of_data_by_sample=box_data,
        pconfig=draw_config
    )

    box_html = remove_subtitle(box_html)

    add_sub_section(
        sub_section=sub_section,
        plot=box_html,
        order=5,
        description="log2(Ms1.Area) for each Run.",
        helptext="""
            [DIA-NN: report.tsv] log2(Ms1.Area) for each Run. Ms1.Area: non-normalised MS1 peak area.
            """
    )

# Distribution of Precursor Charges
def draw_dia_whole_exp_charge(sub_section, df):
    charge_df = df[["Precursor.Charge"]].copy()
    charge_df["Precursor.Charge"] = charge_df["Precursor.Charge"].astype("str")

    bar_data = {
        "Whole Experiment": charge_df[
            "Precursor.Charge"
        ].value_counts().sort_index().to_dict()
    }

    bar_data = {str(k): v for k, v in bar_data.items()}

    draw_config = {
        "id": "distribution_of_precursor_charges",
        "cpswitch": True,
        "title": "Distribution of Precursor Charges",
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
        order=5,
        description="""
            This is a bar chart representing the distribution of the precursor ion charges for a given whole experiment.
        """,
        helptext="""
            [DIA-NN: main report] distribution of the precursor ion charges for a given whole experiment.
            Precursor.Charge: the charge of the precursor.
            """
    )

# Charge-state of Per File
def draw_dia_ms2_charge(sub_section, df):

    df = df[["Precursor.Charge", "Run"]].copy()
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
            [DIA-NN: main report] The distribution of the charge-state of the precursor ion (Precursor.Charge).
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
        "xlab": "Standard Deviation of log2(Precursor.Quantity)",
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
            the standard deviation of log2(Precursor.Quantity).
            """
    )

# DIA-NN: Quantification Table
def draw_diann_quant_table(sub_section, diann_report, sample_df, file_df):

    # Peptides Quantification Table
    peptides_table, peptides_headers = quantms_utils.create_peptides_table(
        diann_report,
        sample_df,
        file_df
    )
    draw_peptides_table(
        sub_section,
        peptides_table,
        peptides_headers,
        "DIA-NN"
    )

    # Protein Quantification Table
    protein_table, protein_headers = quantms_utils.create_protein_table(
        diann_report,
        sample_df,
        file_df
    )
    draw_protein_table(
        sub_section,
        protein_table,
        protein_headers,
        "DIA-NN"
    )

# DIA-NN: Delta Mass
def draw_dia_delta_mass(sub_section, df):

    delta_df = df[df["Ms1.Apex.Mz.Delta"] != 0][["Ms1.Apex.Mz.Delta", "Run"]].copy()
    delta_mass = cal_delta_mass_dict(delta_df, "Ms1.Apex.Mz.Delta")

    x_values = list(delta_mass["count"].keys())
    range_step = (max(x_values) - min(x_values)) * 0.05

    data_label = [
        {
            "name": "Count (All Data)",
            "ylab": "Count",
            "tt_label": "{point.x} Mass delta counts: {point.y}",
            "xmax": max(abs(x) for x in x_values) + range_step,
            "xmin": -(max(abs(x) for x in x_values) + range_step),
            "ymin": 0,
        },
        {
            "name": "Relative Frequency (All Data)",
            "ylab": "Relative Frequency",
            "tt_label": "{point.x} Mass delta relative frequency: {point.y}",
            "xmax": max(abs(x) for x in x_values) + range_step,
            "xmin": -(max(abs(x) for x in x_values) + range_step),
            "ymin": 0,
        },
    ]

    pconfig = {
        "id": "delta_mass",
        "title": "Delta Mass",
        "colors": {"count": "#b2df8a", "relative_frequency": "#b2df8a"},
        "xlab": "Ms1.Apex.Mz.Delta",
        "data_labels": data_label,
        "style": "lines",
    }

    line_html = linegraph.plot(
        [
            {"count": delta_mass["count"]},
            {"relative_frequency": delta_mass["frequency"]}
        ],
        pconfig
    )

    add_sub_section(
        sub_section=sub_section,
        plot=line_html,
        order=3,
        description="""
            This plot is based on the "Ms1.Apex.Mz.Delta" column from the DIA-NN main report.
        """,
        helptext="""
            [DIA-NN: main report] 
            Ms1.Apex.Mz.Delta: difference between observed precursor m/z and the theoretical value.
        """
    )

# DIA-NN: Normalisation Factor over RT
# Normalisation.Factor average ~ RT bin
def draw_norm_factor_rt(sub_section, plot_data):

    draw_config = {
        "id": "normalisation_factor_over_rt",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Normalisation Factor over RT",
        "ymin": 0,
        "tt_decimals": 3,
        "ylab": "Normalisation Factor",
        "xlab": "Retention time [min]",
        "showlegend": True,
    }

    linegraph_html = linegraph.plot(
        data=plot_data,
        pconfig=draw_config
    )

    linegraph_html = remove_subtitle(linegraph_html)

    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=1,
        description="""
            Distribution of Normalisation.Factor with retention time, derived from the main report.
        """,
        helptext="""
            [DIA-NN: main report] Distribution of Normalisation.Factor with retention time (RT) for each run.
            RT: the retention time (RT) of the PSM in minutes. Normalisation.Factor: normalisation factor 
            applied to the precursor in the specific run, 
            i.e. normalised quantity = normalisation factor X non-normalised quantity
        """,
    )

# DIA-NN: FWHM over RT
def draw_fwhm_rt(sub_section, plot_data):

    draw_config = {
        "id": "fwhm_over_rt",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "FWHM over RT",
        "ymin": 0,
        "tt_decimals": 3,
        "ylab": "FWHM",
        "xlab": "Retention time [min]",
        "showlegend": True,
    }

    linegraph_html = linegraph.plot(
        data=plot_data,
        pconfig=draw_config
    )

    linegraph_html = remove_subtitle(linegraph_html)

    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=2,
        description="""
            Distribution of FWHM with retention time, derived from the main report. 
            FWHM: estimated peak width at half-maximum.
        """,
        helptext="""
            [DIA-NN: main report] Distribution of FWHM with retention time (RT) for each run.
            RT: the retention time (RT) of the PSM in minutes. FWHM: estimated peak width at half-maximum; 
            note that the accuracy of such estimates sometimes strongly depends on the DIA cycle time and 
            sample injection amount, i.e. they can only be used to evaluate chromatographic performance in 
            direct comparisons with similar settings, including the scan window; another caveat is that 
            FWHM does not reflect any peak tailing.
        """,
    )

# DIA-NN: Peak Width over RT
def draw_peak_width_rt(sub_section, plot_data):

    draw_config = {
        "id": "peak_width_over_rt",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Peak Width over RT",
        "ymin": 0,
        "tt_decimals": 3,
        "ylab": "Peak Width",
        "xlab": "Retention time [min]",
        "showlegend": True,
    }

    linegraph_html = linegraph.plot(
        data=plot_data,
        pconfig=draw_config
    )

    linegraph_html = remove_subtitle(linegraph_html)

    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=3,
        description="""
            Distribution of peak width with retention time, derived from the main report. 
            Peak Width = RT.Stop - RT.Start.
        """,
        helptext="""
            [DIA-NN: main report] Distribution of peak width with retention time (RT) for each run.
            RT: the retention time (RT) of the PSM in minutes. RT.Start and RT.Stop: peak boundaries.
        """,
    )

# DIA-NN: Absolute RT Error over RT
def draw_rt_error_rt(sub_section, plot_data):

    draw_config = {
        "id": "rt_error_over_rt",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Absolute RT Error over RT",
        "ymin": 0,
        "tt_decimals": 3,
        "ylab": "|RT - Predicted.RT|",
        "xlab": "Retention time [min]",
        "showlegend": True,
    }

    linegraph_html = linegraph.plot(
        data=plot_data,
        pconfig=draw_config
    )

    linegraph_html = remove_subtitle(linegraph_html)

    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=4,
        description="""
            Distribution of rt error with retention time, derived from the main report. 
        """,
        helptext="""
            [DIA-NN: main report] Distribution of peak width with retention time (RT) for each run.
            RT: the retention time (RT) of the PSM in minutes. 
            Predicted.RT: predicted RT based on the iRT.
        """,
    )

# DIA-NN: LOESS RT ~ iRT
def draw_loess_rt_irt(sub_section, plot_data):

    draw_config = {
        "id": "loess_rt_irt",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "LOESS RT ~ iRT",
        "ymin": 0,
        "tt_decimals": 3,
        "ylab": "RT",
        "xlab": "iRT",
        "showlegend": True,
    }

    linegraph_html = linegraph.plot(
        data=plot_data,
        pconfig=draw_config
    )

    linegraph_html = remove_subtitle(linegraph_html)

    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=4,
        description="""
            Distribution of LOESS RT ~ iRT for each run, derived from the main report.
        """,
        helptext="""
            [DIA-NN: main report] Distribution of LOESS RT ~ iRT for each run.
            RT: the retention time (RT) of the PSM in minutes. 
            iRT: reference RT as recorded in the spectral library.
        """,
    )

# mzIdentML: Quantification Table
def draw_mzid_quant_table(sub_section, mzid_mzml_df):

    # Peptides Quantification Table
    peptides_table, peptides_headers = mzidentml_utils.create_peptides_table(
        mzid_mzml_df
    )
    draw_peptides_table(
        sub_section,
        peptides_table,
        peptides_headers,
        "mzIdentML"
    )

    # Protein Quantification Table
    protein_table, protein_headers = mzidentml_utils.create_protein_table(
        mzid_mzml_df
    )
    draw_protein_table(
        sub_section,
        protein_table,
        protein_headers,
        "mzIdentML"
    )


# Draw: Peptides Quantification Table
def draw_peptides_table(sub_section, table_data, headers, report_type):

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

    if report_type == "DIA-NN":
        description_text = """
            This plot shows the quantification information of peptides in the final result (DIA-NN report).
            """
        helptext_text = """
            The quantification information of peptides is obtained from the DIA-NN output file. 
            The table shows the quantitative level and distribution of peptides in different study variables, 
            run and peptiforms. The distribution show all the intensity values in a bar plot above and below 
            the average intensity for all the fractions, runs and peptiforms.

            * BestSearchScore: It is equal to min(1 - Q.Value) for DIA-NN datasets.
            * Average Intensity: Average intensity of each peptide sequence across all conditions (0 or NA ignored).
            * Peptide intensity in each condition (Eg. `CT=Mixture;CN=UPS1;QY=0.1fmol`).
            """
    elif report_type == "mzIdentML":
        description_text = """
            This plot shows the quantification information of peptides in the final result (mzIdentML).
            """
        helptext_text = """
            The quantification information of peptides is obtained from the mzIdentML. 
            The table shows the quantitative level and distribution of peptides in different study variables, 
            run and peptiforms. The distribution show all the intensity values in a bar plot above and below 
            the average intensity for all the fractions, runs and peptiforms.

            * BestSearchScore: It is equal to max(search_engine_score) for mzIdentML datasets.
            * Average Intensity: Average intensity of each peptide sequence (0 or NA ignored).
            """
    else:
        description_text = ""
        helptext_text = ""

    add_sub_section(
        sub_section=sub_section,
        plot=table_html,
        order=1,
        description=description_text,
        helptext=helptext_text
    )

# Draw: Protein Quantification Table
def draw_protein_table(sub_sections, table_data, headers, report_type):

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

    if report_type == "DIA-NN":
        description_text = """
            This plot shows the quantification information of proteins in the final result (DIA-NN report).
            """
        helptext_text = """
            The quantification information of proteins is obtained from the DIA-NN output file.
            The table shows the quantitative level and distribution of proteins in different study variables and run.

            * Peptides_Number: The number of peptides for each protein.
            * Average Intensity: Average intensity of each protein across all conditions (0 or NA ignored).
            * Protein intensity in each condition (Eg. `CT=Mixture;CN=UPS1;QY=0.1fmol`): Summarize intensity of peptides.
            """
    elif report_type == "mzIdentML":
        description_text = """
            This plot shows the quantification information of proteins in the final result (mzIdentML).
            """
        helptext_text = """
            The quantification information of proteins is obtained from the mzIdentML.
            The table shows the quantitative level and distribution of proteins in different study variables and run.

            * Peptides_Number: The number of peptides for each protein.
            * Average Intensity: Average intensity of each protein(0 or NA ignored).
            """
    else:
        description_text = ""
        helptext_text = ""

    add_sub_section(
        sub_section=sub_sections,
        plot=table_html,
        order=2,
        description=description_text,
        helptext=helptext_text
    )