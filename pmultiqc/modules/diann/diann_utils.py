import logging
import os
import pandas as pd

from pmultiqc.modules.common.stats import qual_uniform

import numpy as np
from pmultiqc.modules.diann import diann_plots
from pmultiqc.modules.common.plots import (
    draw_peptides_table,
    draw_protein_table,
    draw_ids_rt_count
)
from pmultiqc.modules.common.file_utils import drop_empty_row
from pmultiqc.modules.quantms.quantms_utils import (
    create_peptides_table,
    create_protein_table
)
from pmultiqc.modules.common.utils import evidence_rt_count
from statsmodels.nonparametric.smoothers_lowess import lowess


DEFAULT_BINS = 500

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def get_diann_path(find_log_files):
    """
    Find and return the DIA-NN report file path.
    
    Args:
        find_log_files: Function to search for log files
        
    Returns:
        tuple: (diann_report_path, enable_dia) where:
            - diann_report_path: Full path to the DIA-NN report file
            - enable_dia: Boolean indicating if DIA mode is enabled
            
    Raises:
        ValueError: If no DIA-NN report file is found
    """
    # DIA-NN report file types to search for (in order of preference)
    file_types = ["pmultiqc/diann_report_parquet", "pmultiqc/diann_report_tsv"]
    
    for file_type in file_types:
        for file_info in find_log_files(file_type, filecontents=False):
            diann_report_path = os.path.join(file_info["root"], file_info["fn"])
            return diann_report_path, True
    
    # If no report file found, raise an error
    raise ValueError("DIA-NN report not found. Please check your data!")

def draw_dia_intensitys(sub_section, report_df):

    df_sub = report_df[report_df["Precursor.Quantity"] > 0].copy()
    df_sub["log_intensity"] = np.log2(df_sub["Precursor.Quantity"])
    
    diann_plots.draw_dia_intensity_dis(sub_section, df_sub)

    if diann_plots.can_groupby_for_std(report_df, "Run"):
        diann_plots.draw_dia_intensity_std(sub_section, df_sub)


def draw_dia_heatmap(sub_section, report_df, heatmap_color):

    log.info("Compute the Heatmap.")
    heatmap_data = cal_dia_heatmap(report_df)
    diann_plots.draw_heatmap(sub_section, heatmap_color, heatmap_data)
    log.info("Heatmap calculation is done.")


def draw_dia_ms1(sub_section, df):

    # Ms1.Area: non-normalised MS1 peak area
    if "Ms1.Area" in df.columns:
        df_sub = df[df["Ms1.Area"] > 0][["Ms1.Area", "Run"]].copy()
        if len(df_sub) > 0:
            df_sub["log_ms1_area"] = np.log2(df_sub["Ms1.Area"])
            diann_plots.draw_dia_ms1_area(sub_section, df_sub)

def draw_dia_ms2s(sub_section, df):

    # Distribution of Precursor Charges
    if "Precursor.Charge" in df.columns:
        diann_plots.draw_dia_whole_exp_charge(sub_section, df)

    # Charge-state of Per File
    diann_plots.draw_dia_ms2_charge(sub_section, df)


def draw_dia_mass_error(sub_section, df):

    # Ms1.Apex.Mz.Delta: difference between observed precursor m/z and the theoretical value
    if "Ms1.Apex.Mz.Delta" in df.columns:
        diann_plots.draw_dia_delta_mass(sub_section, df)

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
        norm_factor_rt = cal_feature_avg_rt(df, "Normalisation.Factor")
        diann_plots.draw_norm_factor_rt(sub_section, norm_factor_rt)

    # 2. FWHM over RT
    #   FWHM: estimated peak width at half-maximum
    if "FWHM" in df.columns:
        log.info("Draw[rt_qc]: draw_fwhm_rt")
        fwhm_rt = cal_feature_avg_rt(df, "FWHM")
        diann_plots.draw_fwhm_rt(sub_section, fwhm_rt)
    
    # 3. Peak Width over RT
    #   RT.Start and RT.Stop peak boundaries
    if all(col in df.columns for col in ["RT.Start", "RT.Stop"]):
        log.info("Draw[rt_qc]: draw_peak_width_rt")
        df["peak_width"] = df["RT.Stop"] - df["RT.Start"]
        peak_width_rt = cal_feature_avg_rt(df, "peak_width")
        diann_plots.draw_peak_width_rt(sub_section, peak_width_rt)

    # 4. Absolute RT Error over RT
    if all(col in df.columns for col in ["RT", "Predicted.RT"]):
        log.info("Draw[rt_qc]: draw_rt_error_rt")
        df["rt_error"] = abs(df["RT"] - df["Predicted.RT"])
        rt_error_rt = cal_feature_avg_rt(df, "rt_error")
        diann_plots.draw_rt_error_rt(sub_section, rt_error_rt)

    # 5. loess(RT ~ iRT)
    if all(col in df.columns for col in ["RT", "iRT"]):
        log.info("Draw[rt_qc]: draw_loess_rt_irt")
        rt_irt_loess = cal_rt_irt_loess(df)
        if rt_irt_loess is not None:
            diann_plots.draw_loess_rt_irt(sub_section, rt_irt_loess)

# DIA-NN: Quantification Table
def draw_diann_quant_table(sub_section, diann_report, sample_df, file_df):

    # Peptides Quantification Table
    peptides_table, peptides_headers = create_peptides_table(
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
    protein_table, protein_headers = create_protein_table(
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


# DIA-NN: HeatMap
def cal_dia_heatmap(report_df):

    # "Contaminants" & "Peptide Intensity"
    pep_intensity = heatmap_cont_pep_intensity(report_df)

    # missed tryptic cleavages: there is no available data

    return pep_intensity


def heatmap_cont_pep_intensity(report_df):

    df = report_df[
        [
            "Run", "Protein.Names", "Precursor.Quantity", "RT", "Predicted.RT", 
            "Precursor.Charge", "Normalisation.Factor", "RT.Stop", "RT.Start"
        ]
    ].copy()

    # TODO "CON"?
    df["is_contaminant"] = df["Protein.Names"].str.startswith("CON", na=False)

    # 3. "Charge"
    charge = dict()
    for raw_file, group in df[["Run", "Precursor.Charge"]].groupby("Run"):
        charge[raw_file] = group["Precursor.Charge"].value_counts()[2] / len(group)
    charge_median = np.median(list(charge.values()))
    heatmap_charge = dict(
        zip(
            charge.keys(),
            list(map(lambda v: float(1 - np.abs(v - charge_median)), charge.values())),
        )
    )

    heatmap_dict = {}
    for run, group in df.groupby("Run"):
        
        # 1. "Contaminants"
        cont_intensity_sum = group[group["is_contaminant"]]["Precursor.Quantity"].sum()
        if np.isnan(cont_intensity_sum) or cont_intensity_sum == 0:
            contaminant = 1
        else:
            intensity_sum = group["Precursor.Quantity"].sum()
            contaminant = 1 - cont_intensity_sum / intensity_sum

        # 2. "Peptide Intensity"
        pep_median = np.nanmedian(group["Precursor.Quantity"].to_numpy())
        pep_intensity = float(np.fmin(1.0, pep_median / (2**23)))


        # 4. "RT Alignment"
        rt_alignment = max(
            0.0,
            1 - float(np.mean(np.abs(group["RT"] - group["Predicted.RT"])))
        )

        # 5. "ID rate over RT"
        ids_rate_over_rt = qual_uniform(group["RT"])

        # 6. Normalization Factor MAD
        def mean_abs_dev(x):
            mean_x = x.mean()
            return float(1- np.mean(np.abs(x - mean_x)))

        norm_factor_mad = mean_abs_dev(group["Normalisation.Factor"])

        # 7. Peak Width = RT.Stop - RT.Start
        peak_width = max(
            0.0,
            1 - float(np.mean(group["RT.Stop"] - group["RT.Start"]))
        )

        # All Dict
        heatmap_dict[run] = {
            "Contaminants": contaminant,
            "Peptide Intensity": pep_intensity,
            "Charge": heatmap_charge.get(run, 0),
            "RT Alignment": rt_alignment,
            "ID rate over RT": ids_rate_over_rt,
            "Norm Factor": norm_factor_mad,
            "Peak Width": peak_width,
        }

    return heatmap_dict

# DIA-NN: IDs over RT
def draw_dia_ids_rt(sub_section, report_df):

    rt_df = report_df[["Run", "RT"]].copy()
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

# DIA-NN: Lowess (Loess)
def cal_rt_irt_loess(report_df, frac=0.3, data_bins: int=DEFAULT_BINS):

    if len(report_df) > 1000000:
        log.warning(f"Dataset too large ({len(report_df)} rows). Skipping LOWESS computation.")
        return None

    log.info("Start compute loess...")
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
