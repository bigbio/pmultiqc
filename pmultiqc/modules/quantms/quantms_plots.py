import itertools
import logging

import numpy as np
from multiqc.plots import table

from . import quantms_utils, dia_plots
from ..common.common_plots import draw_ids_rt_count
from ..core.section_groups import add_sub_section
from ..maxquant.maxquant_utils import evidence_rt_count

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def draw_dia_heatmap(sub_section, report_df, heatmap_color):

    log.info("Compute the Heatmap.")
    heatmap_data = quantms_utils.cal_dia_heatmap(report_df)
    dia_plots.draw_heatmap(sub_section, heatmap_color, heatmap_data)
    log.info("Heatmap calculation is done.")


def draw_dia_intensitys(sub_section, report_df):

    df_sub = report_df[report_df["Precursor.Quantity"] > 0].copy()
    df_sub["log_intensity"] = np.log2(df_sub["Precursor.Quantity"])

    dia_plots.draw_dia_intensity_dis(sub_section, df_sub)

    if dia_plots.can_groupby_for_std(report_df, "Run"):
        dia_plots.draw_dia_intensity_std(sub_section, df_sub)


def draw_dia_ms1(sub_section, df):

    # Ms1.Area: non-normalised MS1 peak area
    if "Ms1.Area" in df.columns:
        df_sub = df[df["Ms1.Area"] > 0][["Ms1.Area", "Run"]].copy()
        if len(df_sub) > 0:
            df_sub["log_ms1_area"] = np.log2(df_sub["Ms1.Area"])
            dia_plots.draw_dia_ms1_area(sub_section, df_sub)


def draw_dia_ms2s(sub_section, df):

    # Distribution of Precursor Charges
    if "Precursor.Charge" in df.columns:
        dia_plots.draw_dia_whole_exp_charge(sub_section, df)

    # Charge-state of Per File
    dia_plots.draw_dia_ms2_charge(sub_section, df)


def draw_dia_mass_error(sub_section, df):

    # Ms1.Apex.Mz.Delta: difference between observed precursor m/z and the theoretical value
    if "Ms1.Apex.Mz.Delta" in df.columns:
        dia_plots.draw_dia_delta_mass(sub_section, df)


# DIA-NN: RT Quality Control
def draw_dia_rt_qc(sub_section, report_df):

    df = report_df.copy()

    # IDs over RT
    draw_dia_ids_rt(sub_section, df)

    # 1. Normalisation Factor over RT
    #   Normalisation.Factor: normalisation factor applied to the precursor in the specific run,
    #   i.e. normalised quantity = normalisation factor X non-normalised quantity
    if "Normalisation.Factor" in df.columns:
        log.info("Draw[rt_qc]: norm_factor_rt")
        norm_factor_rt = quantms_utils.cal_feature_avg_rt(df, "Normalisation.Factor")
        dia_plots.draw_norm_factor_rt(sub_section, norm_factor_rt)

    # 2. FWHM over RT
    #   FWHM: estimated peak width at half-maximum
    if "FWHM" in df.columns:
        log.info("Draw[rt_qc]: draw_fwhm_rt")
        fwhm_rt = quantms_utils.cal_feature_avg_rt(df, "FWHM")
        dia_plots.draw_fwhm_rt(sub_section, fwhm_rt)

    # 3. Peak Width over RT
    #   RT.Start and RT.Stop peak boundaries
    if all(col in df.columns for col in ["RT.Start", "RT.Stop"]):
        log.info("Draw[rt_qc]: draw_peak_width_rt")
        df["peak_width"] = df["RT.Stop"] - df["RT.Start"]
        peak_width_rt = quantms_utils.cal_feature_avg_rt(df, "peak_width")
        dia_plots.draw_peak_width_rt(sub_section, peak_width_rt)

    # 4. Absolute RT Error over RT
    if all(col in df.columns for col in ["RT", "Predicted.RT"]):
        log.info("Draw[rt_qc]: draw_rt_error_rt")
        df["rt_error"] = abs(df["RT"] - df["Predicted.RT"])
        rt_error_rt = quantms_utils.cal_feature_avg_rt(df, "rt_error")
        dia_plots.draw_rt_error_rt(sub_section, rt_error_rt)

    # 5. loess(RT ~ iRT)
    if all(col in df.columns for col in ["RT", "iRT"]):
        log.info("Draw[rt_qc]: draw_loess_rt_irt")
        rt_irt_loess = quantms_utils.cal_rt_irt_loess(df)
        if rt_irt_loess is not None:
            dia_plots.draw_loess_rt_irt(sub_section, rt_irt_loess)


# DIA-NN: IDs over RT
def draw_dia_ids_rt(sub_section, report_df):

    rt_df = report_df[["Run", "RT"]].copy()
    rt_df.rename(columns={"Run": "raw file", "RT": "retention time"}, inplace=True)
    ids_over_rt = evidence_rt_count(rt_df)
    draw_ids_rt_count(sub_section, ids_over_rt, "dia")


# DIA-NN: Quantification Table
def draw_diann_quant_table(sub_section, diann_report, sample_df, file_df):

    # Peptides Quantification Table
    peptides_table, peptides_headers = quantms_utils.create_peptides_table(
        diann_report, sample_df, file_df
    )
    draw_peptides_table(sub_section, peptides_table, peptides_headers, "DIA-NN")

    # Protein Quantification Table
    protein_table, protein_headers = quantms_utils.create_protein_table(
        diann_report, sample_df, file_df
    )
    draw_protein_table(sub_section, protein_table, protein_headers, "DIA-NN")


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
        helptext=helptext_text,
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
        pconfig=draw_config,
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
        helptext=helptext_text,
    )