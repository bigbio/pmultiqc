import re
import pandas as pd

from multiqc.plots import heatmap, box, bargraph, linegraph

from pmultiqc.modules.common.plots.general import (
    plot_html_check,
    plot_data_check
)
from pmultiqc.modules.common.stats import cal_delta_mass_dict
from pmultiqc.modules.core.section_groups import add_sub_section
from pmultiqc.modules.common.plots import dia as dia_plots
from pmultiqc.modules.common.common_utils import group_charge
from pmultiqc.modules.common.logging import get_logger

log = get_logger("pmultiqc.modules.common.plots.dia")

# DIA-NN: HeatMap
def draw_heatmap(sub_section, hm_colors, heatmap_data):

    pconfig = {
        "id": "heatmap",
        "title": "HeatMap",
        "min": 0,
        "max": 1,
        "xlab": "Metrics",
        "ylab": "RawName",
        "zlab": "Score",
        "tt_decimals": 4,
        "square": False,
        "colstops": hm_colors,
        "cluster_rows": False,
        "cluster_cols": False,
        "save_data_file": False,
    }

    hm_html = heatmap.plot(data=heatmap_data, pconfig=pconfig)

    hm_html = plot_html_check(hm_html)

    add_sub_section(
        sub_section=sub_section,
        plot=hm_html,
        order=2,
        description="""
            This heatmap provides an overview of the performance of the quantms DIA (DIA-NN) results.
        """,
        helptext="""
            This plot shows the pipeline performance overview. Some metrics are calculated.
            *Heatmap score[Contaminants]: as fraction of summed intensity with 0 = sample full of contaminants; 
                1 = no contaminants
            *Heatmap score[Pep Intensity (>23.0)]: Linear scale of the median intensity reaching the threshold, 
                i.e. reaching 2^21 of 2^23 gives score 0.25.
            *Heatmap score[Charge]: Deviation of the charge 2 proportion from a representative 
                Raw file (median). For typtic digests, peptides of charge 2 (one N-terminal and one at 
                tryptic C-terminal R or K residue) should be dominant. Ionization issues (voltage?), 
                in-source fragmentation, missed cleavages and buffer irregularities can cause a shift 
                (see Bittremieux 2017, DOI: 10.1002/mas.21544)
            *Heatmap score[RT Alignment]: Compute 1 minus the mean absolute difference between 'RT' and 'Predicted.RT', 
                and take the maximum of this value and 0. 1: |RT - Predicted.RT| = 0
            *Heatmap score [ID rate over RT]: Judge column occupancy over retention time. 
                Ideally, the LC gradient is chosen such that the number of identifications 
                (here, after FDR filtering) is uniform over time, to ensure consistent instrument duty cycles. 
                Sharp peaks and uneven distribution of identifications over time indicate potential for LC gradient 
                optimization.Scored using 'Uniform' scoring function. i.e. constant receives good score, extreme shapes are bad
            *Heatmap score [Norm Factor]: Computes the mean absolute deviation (MAD) of 'Normalisation.Factor' from its mean.
                0 = high variability in normalization factors; 1 = perfectly consistent normalization factors
            *Heatmap score [Peak Width]: Average peak width (RT.Stop - RT.Start). 1 = peak width equals 0; 
                0 = peak width equals 1 or greater
        """,
    )


# Intensity Distribution
def draw_dia_intensity_dis(sub_section, df, sdrf_file_df):

    df_sub = df[["Run", "Modified.Sequence", "Protein.Group", "log_intensity"]].copy()

    if not sdrf_file_df.empty:

        df_sub = df_sub.merge(
            sdrf_file_df[["Sample", "Run"]].drop_duplicates(),
            on="Run"
        )

        df_sub["Sample"] = df_sub["Sample"].astype(int)

        box_data = [
            {
                (
                    f"Sample {str(run)}"
                    if data_type == "Sample"
                    else str(run)
                ): group["log_intensity"].dropna().tolist()
                for run, group in df_sub.groupby(data_type, sort=True)
            }
            for data_type in ["Run", "Sample"]
        ]

        plot_label = ["by Run", "by Sample"]

    else:
        box_data = [
            {
                str(run): group["log_intensity"].dropna().tolist()
                for run, group in df.groupby("Run")
            }
        ]
        plot_label = ["by Run"]

    draw_config = {
        "id": "intensity_distribution_box",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Intensity Distribution",
        "tt_decimals": 5,
        "xlab": "log2(Precursor.Quantity)",
        "data_labels": plot_label,
        "sort_samples": False,
        "save_data_file": False,
    }

    box_html = box.plot(list_of_data_by_sample=box_data, pconfig=draw_config)

    # box_html.flat
    box_html = plot_data_check(
        plot_data=box_data,
        plot_html=box_html,
        log_text="pmultiqc.modules.common.plots.dia",
        function_name="draw_dia_intensity_dis"
    )
    box_html = plot_html_check(box_html)

    add_sub_section(
        sub_section=sub_section,
        plot=box_html,
        order=3,
        description="log2(Precursor.Quantity) for each Run (or Sample).",
        helptext="""
            [DIA-NN: main report] log2(Precursor.Quantity) for each Run (or Sample).
            """,
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
        "save_data_file": False,
    }

    box_html = box.plot(list_of_data_by_sample=box_data, pconfig=draw_config)

    # box_html.flat
    box_html = plot_data_check(
        plot_data=box_data,
        plot_html=box_html,
        log_text="pmultiqc.modules.common.plots.dia",
        function_name="draw_dia_ms1_area"
    )
    box_html = plot_html_check(box_html)

    add_sub_section(
        sub_section=sub_section,
        plot=box_html,
        order=5,
        description="log2(Ms1.Area) for each Run.",
        helptext="""
            [DIA-NN: report.tsv] log2(Ms1.Area) for each Run. Ms1.Area: non-normalised MS1 peak area.
            """,
    )

# Distribution of Precursor Charges
def draw_dia_whole_exp_charge(sub_section, df):
    charge_df = df[["Precursor.Charge"]].copy()
    charge_df["Precursor.Charge"] = charge_df["Precursor.Charge"].astype("str")

    bar_data = {
        "Whole Experiment": charge_df["Precursor.Charge"].value_counts().sort_index().to_dict()
    }

    bar_data = {str(k): v for k, v in bar_data.items()}

    draw_config = {
        "id": "distribution_of_precursor_charges",
        "cpswitch": True,
        "title": "Distribution of Precursor Charges",
        "tt_decimals": 0,
        "ylab": "Count",
        "save_data_file": False,
    }

    bar_html = bargraph.plot(
        data=bar_data,
        pconfig=draw_config,
    )

    bar_html = plot_html_check(bar_html)

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
            """,
    )


# Charge-state of Per File
def draw_dia_ms2_charge(sub_section, df, sdrf_file_df):

    df = df[["Precursor.Charge", "Run"]].copy()
    df["Precursor.Charge"] = df["Precursor.Charge"].astype("str")

    stat_data_by_run = group_charge(df, "Run", "Precursor.Charge")

    if sdrf_file_df.empty:

        bar_data = stat_data_by_run.to_dict(orient="index")
        plot_label = ["by Run"]

    else:
        df = df.merge(
            right=sdrf_file_df[["Sample", "Run"]].drop_duplicates(),
            on="Run"
        )

        stat_data_by_sample = group_charge(df, "Sample", "Precursor.Charge")

        bar_data = [
            stat_data_by_run.to_dict(orient="index"),
            stat_data_by_sample.to_dict(orient="index")
        ]
        plot_label = ["by Run", "by Sample"]

    draw_config = {
        "id": "charge_state_of_per_file",
        "cpswitch": True,
        "title": "Charge-state of Per File",
        "tt_decimals": 0,
        "ylab": "Count",
        "data_labels": plot_label,
        "save_data_file": False,
    }

    bar_html = bargraph.plot(
        data=bar_data,
        pconfig=draw_config,
    )

    bar_html = plot_html_check(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=6,
        description="The distribution of the charge-state of the precursor ion.",
        helptext="""
            [DIA-NN: main report] The distribution of the charge-state of the precursor ion (Precursor.Charge).
            """,
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


def draw_dia_intensity_std(sub_section, df, sdrf_file_df):

    box_data = calculate_dia_intensity_std(df, sdrf_file_df)

    if not box_data:
        return

    draw_box_config = {
        "id": "dia_std_intensity_box",
        "title": "Standard Deviation of Intensity",
        "cpswitch": False,
        "tt_decimals": 5,
        "xlab": "Standard Deviation of log2(Precursor.Quantity)",
        "save_data_file": False,
    }

    box_html = box.plot(
        list_of_data_by_sample=box_data,
        pconfig=draw_box_config,
    )

    # box_html.flat
    box_html = plot_data_check(
        plot_data=box_data,
        plot_html=box_html,
        log_text="pmultiqc.modules.common.plots.dia",
        function_name="draw_dia_intensity_std"
    )
    box_html = plot_html_check(box_html)

    add_sub_section(
        sub_section=sub_section,
        plot=box_html,
        order=6,
        description="Standard deviation of intensity by sample (experimental conditions).",
        helptext="""
            [DIA-NN: report.tsv] Sample grouping is derived from the SDRF when available; 
            otherwise, it is parsed from "Run" names. 
            First, identify the experimental conditions from the "Run" name. 
            Then, group the data by experimental condition and Modified.Sequence, and calculate 
            the standard deviation of log2(Precursor.Quantity).
            """,
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
        "save_data_file": False,
    }

    line_html = linegraph.plot(
        [{"count": delta_mass["count"]}, {"relative_frequency": delta_mass["frequency"]}], pconfig
    )

    add_sub_section(
        sub_section=sub_section,
        plot=line_html,
        order=1,
        description="""
            This plot is based on the "Ms1.Apex.Mz.Delta" column from the DIA-NN main report.
        """,
        helptext="""
            [DIA-NN: main report] 
            Ms1.Apex.Mz.Delta: difference between observed precursor m/z and the theoretical value.
        """,
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
        "save_data_file": False,
    }

    linegraph_html = linegraph.plot(data=plot_data, pconfig=draw_config)

    linegraph_html = plot_html_check(linegraph_html)

    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=2,
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
        "save_data_file": False,
    }

    linegraph_html = linegraph.plot(data=plot_data, pconfig=draw_config)

    linegraph_html = plot_html_check(linegraph_html)

    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=3,
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
        "save_data_file": False,
    }

    linegraph_html = linegraph.plot(data=plot_data, pconfig=draw_config)

    linegraph_html = plot_html_check(linegraph_html)

    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=4,
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
        "save_data_file": False,
    }

    linegraph_html = linegraph.plot(data=plot_data, pconfig=draw_config)

    linegraph_html = plot_html_check(linegraph_html)

    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=5,
        description="""
            Distribution of rt error with retention time, derived from the main report. 
        """,
        helptext="""
            [DIA-NN: main report] Distribution of absolute RT error (|RT - Predicted.RT|) with retention time (RT) for each run.
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
        "save_data_file": False,
    }

    linegraph_html = linegraph.plot(data=plot_data, pconfig=draw_config)

    linegraph_html = plot_html_check(linegraph_html)

    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=6,
        description="""
            Distribution of LOESS RT ~ iRT for each run, derived from the main report.
        """,
        helptext="""
            [DIA-NN: main report] Distribution of LOESS RT ~ iRT for each run.
            RT: the retention time (RT) of the PSM in minutes. 
            iRT: reference RT as recorded in the spectral library.
        """,
    )

def calculate_dia_intensity_std(df, sdrf_file_df):

    df_sub = df[["Run", "Modified.Sequence", "Protein.Group", "log_intensity"]].copy()

    if not sdrf_file_df.empty:

        df_sub = df_sub.merge(
            sdrf_file_df[["Sample", "Run"]].drop_duplicates(),
            on="Run"
        )

        grouped_std = (
            df_sub.groupby(["Sample", "Modified.Sequence"])["log_intensity"]
            .std()
            .reset_index(name="log_intensity_std")
        )

        plot_data = {
            f"Sample {str(sample)}": group["log_intensity_std"].dropna().tolist()
            for sample, group in grouped_std.groupby("Sample")
        }

        return plot_data

    if dia_plots.can_groupby_for_std(df_sub, "Run"):

        log.info("No SDRF available; experimental grouping was parsed using Run names.")

        df_sub[["run_condition", "run_replicate"]] = df_sub["Run"].apply(
            lambda x: pd.Series(extract_condition_and_replicate(x))
        )

        grouped_std = (
            df_sub.groupby(["run_condition", "Modified.Sequence"])["log_intensity"]
            .std()
            .reset_index(name="log_intensity_std")
        )

        plot_data = {
            condition: group["log_intensity_std"].dropna().tolist()
            for condition, group in grouped_std.groupby("run_condition")
        }

        return plot_data
    else:
        log.warning("No SDRF available; failed to parse experimental groups; SD Intensity not generated.")

def extract_condition_and_replicate(run_name):

    match = re.search(r"^(.*?)([A-Za-z]*)(\d+)$", run_name)

    if match:
        condition_base = (match.group(1) + match.group(2)).rstrip("_")
        replicate = int(match.group(3))
        return condition_base, replicate
    else:
        log.warning("Failed to parse condition/replicate from Run='%s' in DIA report.tsv", run_name)
        # Fallback: keep full run name as condition, unknown replicate
        return run_name, None

# re-export by moving file; contents will be identical after move.
