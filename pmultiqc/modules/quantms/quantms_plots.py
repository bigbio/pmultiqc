import logging
from multiqc.plots import bargraph, box, table

from ..common.common_plots import (
    draw_ids_rt_count,
    remove_subtitle
)
from ..core.section_groups import add_sub_section
from ..maxquant.maxquant_utils import evidence_rt_count
from ..common.common_plots import draw_oversampling
from . import quantms_utils, mzidentml_utils

import re
import numpy as np
import pandas as pd
import itertools


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

def draw_dia_intensitys(sub_section, report_df):

    # Use non-normalized Precursor.Quantity instead of Precursor.Normalised
    # Also try Ms1.Area if available
    intensity_col = None
    if "Precursor.Quantity" in report_df.columns:
        intensity_col = "Precursor.Quantity"
    elif "Ms1.Area" in report_df.columns:
        intensity_col = "Ms1.Area"
    elif "Precursor.Normalised" in report_df.columns:
        # Fallback to normalized if non-normalized not available
        intensity_col = "Precursor.Normalised"
    
    if intensity_col is None:
        return  # No intensity column available
    
    # Filter out zero values before log-transform
    df_sub = report_df[report_df[intensity_col] > 0].copy()
    df_sub["log_intensity"] = np.log2(df_sub[intensity_col])
    df_sub["intensity_column"] = intensity_col  # Store which column was used
    
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
    
    # Add missed cleavages distribution if the peptide sequence information is available
    draw_dia_missed_cleavages(sub_section, df)

def draw_dia_time_mass(sub_section, df):

    # First draw RT distribution (existing functionality)
    df_rt = df[["Run", "RT"]].copy()
    df_rt.rename(
        columns={
            "Run": "raw file",
            "RT": "retention time"
        },
        inplace=True
    )

    ids_over_rt = evidence_rt_count(df_rt)

    draw_ids_rt_count(
        sub_section,
        ids_over_rt,
        "dia"
    )

    # Add mass error trends if available (DIA-NN 2.1+)
    # MS1 and MS2 mass errors should exclude zero values
    mass_error_cols = []
    if "Mass.Error..ppm." in df.columns:
        mass_error_cols.append("Mass.Error..ppm.")
    if "MS1.Mass.Error..ppm." in df.columns:
        mass_error_cols.append("MS1.Mass.Error..ppm.")
    if "MS2.Mass.Error..ppm." in df.columns:
        mass_error_cols.append("MS2.Mass.Error..ppm.")
    
    # Also check for older naming conventions
    if "Mass.Error.ppm" in df.columns:
        mass_error_cols.append("Mass.Error.ppm")
    if "MS1.Mass.Error.ppm" in df.columns:
        mass_error_cols.append("MS1.Mass.Error.ppm")
    if "MS2.Mass.Error.ppm" in df.columns:
        mass_error_cols.append("MS2.Mass.Error.ppm")

    for mass_error_col in mass_error_cols:
        draw_dia_mass_error_trends(sub_section, df, mass_error_col)
    
    # Add new RT quality control plots if columns are available
    draw_rt_quality_plots(sub_section, df)

def draw_dia_mass_error_trends(sub_section, df, mass_error_col):
    """Draw mass error trends vs RT or m/z bins"""
    
    # Filter out zero values as they may mean 'no data'
    df_filtered = df[df[mass_error_col] != 0].copy()
    
    if df_filtered.empty:
        return  # No data to plot
    
    # Bin by RT (use 10-minute bins)
    rt_bins = np.arange(0, df_filtered["RT"].max() + 10, 10)
    df_filtered["RT_bin"] = pd.cut(df_filtered["RT"], bins=rt_bins, labels=rt_bins[:-1])
    
    # Calculate average mass error per RT bin for each run
    mass_error_vs_rt = {}
    for run, group in df_filtered.groupby("Run"):
        rt_binned = group.groupby("RT_bin")[mass_error_col].mean().dropna()
        if len(rt_binned) > 0:
            mass_error_vs_rt[str(run)] = rt_binned.to_dict()
    
    if not mass_error_vs_rt:
        return  # No data to plot
        
    # Convert to format expected by linegraph
    plot_data = {}
    for run, data in mass_error_vs_rt.items():
        plot_data[run] = {float(rt_bin): mass_error for rt_bin, mass_error in data.items()}
    
    draw_config = {
        "id": f"mass_error_vs_rt_{mass_error_col.replace('.', '_').replace(' ', '_')}",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": f"{mass_error_col.replace('.', ' ')} vs RT",
        "ylab": f"{mass_error_col.replace('.', ' ')} (ppm)",
        "xlab": "Retention Time (min)",
        "tt_decimals": 3,
        "showlegend": True,
        "style": "lines+markers",
    }

    from multiqc.plots import linegraph
    
    line_html = linegraph.plot(
        data=plot_data,
        pconfig=draw_config
    )

    line_html = remove_subtitle(line_html)

    # Determine error type for description
    error_type = "MS1" if "MS1" in mass_error_col else "MS2" if "MS2" in mass_error_col else "Mass"
    
    add_sub_section(
        sub_section=sub_section,
        plot=line_html,
        order=7 if "MS1" in mass_error_col else 8,
        description=f"{error_type} mass error trends over retention time.",
        helptext=f"""
            [DIA-NN: report.tsv] {error_type} mass error in ppm averaged across retention time bins. 
            Zero values are excluded as they may indicate missing data. Available in DIA-NN 2.1+.
            """
    )

def draw_rt_quality_plots(sub_section, df):
    """Draw RT quality control plots as requested in the issue"""
    
    # 1. Normalisation Factor average ~ RT bin (DIA-NN 2.0+)
    if "Normalisation.Factor" in df.columns:
        draw_rt_binned_average(sub_section, df, "Normalisation.Factor", "RT", 
                               "Normalisation Factor vs RT", 9)
    
    # Alternative calculation for older versions: Normalisation.Factor = Precursor.Normalised/Precursor.Quantity
    elif "Precursor.Normalised" in df.columns and "Precursor.Quantity" in df.columns:
        df_norm = df[(df["Precursor.Normalised"] > 0) & (df["Precursor.Quantity"] > 0)].copy()
        df_norm["Normalisation.Factor"] = df_norm["Precursor.Normalised"] / df_norm["Precursor.Quantity"]
        draw_rt_binned_average(sub_section, df_norm, "Normalisation.Factor", "RT", 
                               "Normalisation Factor vs RT (calculated)", 9)
    
    # 2. FWHM average ~ RT bin (DIA-NN 2.0+)
    if "FWHM" in df.columns:
        draw_rt_binned_average(sub_section, df, "FWHM", "RT", "FWHM vs RT", 10)
    
    # 3. (Stop.RT - Start.RT) average ~ RT bin (DIA-NN 2.0+)
    if "Stop.RT" in df.columns and "Start.RT" in df.columns:
        df_rt_width = df[(df["Stop.RT"] > 0) & (df["Start.RT"] > 0)].copy()
        df_rt_width["RT.Width"] = df_rt_width["Stop.RT"] - df_rt_width["Start.RT"]
        draw_rt_binned_average(sub_section, df_rt_width, "RT.Width", "RT", "RT Width vs RT", 11)
    
    # 4. abs(RT - Predicted.RT) average vs RT bin
    if "Predicted.RT" in df.columns:
        df_rt_pred = df[df["Predicted.RT"] > 0].copy()
        df_rt_pred["RT.Deviation"] = abs(df_rt_pred["RT"] - df_rt_pred["Predicted.RT"])
        draw_rt_binned_average(sub_section, df_rt_pred, "RT.Deviation", "RT", 
                               "RT Deviation vs RT", 12)

def draw_rt_binned_average(sub_section, df, y_col, x_col, title, order):
    """Helper function to draw RT binned averages"""
    
    # Create RT bins (5-minute bins for better resolution)
    rt_bins = np.arange(0, df[x_col].max() + 5, 5)
    df_binned = df.copy()
    df_binned[f"{x_col}_bin"] = pd.cut(df_binned[x_col], bins=rt_bins, labels=rt_bins[:-1])
    
    # Calculate average per RT bin for each run
    plot_data = {}
    for run, group in df_binned.groupby("Run"):
        rt_binned = group.groupby(f"{x_col}_bin")[y_col].mean().dropna()
        if len(rt_binned) > 0:
            plot_data[str(run)] = {float(rt_bin): avg_val for rt_bin, avg_val in rt_binned.items()}
    
    if not plot_data:
        return  # No data to plot
    
    draw_config = {
        "id": f"rt_quality_{y_col.replace('.', '_').lower()}",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": title,
        "ylab": f"Average {y_col}",
        "xlab": "Retention Time (min)",
        "tt_decimals": 3,
        "showlegend": True,
        "style": "lines+markers",
    }

    from multiqc.plots import linegraph
    
    line_html = linegraph.plot(
        data=plot_data,
        pconfig=draw_config
    )

    line_html = remove_subtitle(line_html)

    # Generate description based on the metric
    if "Normalisation" in y_col:
        description = "Normalisation factor trends over retention time bins."
        helptext = """
            [DIA-NN: report.tsv] Average normalisation factor across RT bins for each run.
            Shows how normalisation varies across the gradient. Available in DIA-NN 2.0+.
            """
    elif "FWHM" in y_col:
        description = "Full Width at Half Maximum (FWHM) trends over retention time."
        helptext = """
            [DIA-NN: report.tsv] Average peak width (FWHM) across RT bins for each run.
            Indicates chromatographic peak quality across the gradient. Available in DIA-NN 2.0+.
            """
    elif "Width" in y_col:
        description = "RT peak width (Stop.RT - Start.RT) trends over retention time."
        helptext = """
            [DIA-NN: report.tsv] Average RT peak width across RT bins for each run.
            Shows chromatographic peak shape consistency. Available in DIA-NN 2.0+.
            """
    elif "Deviation" in y_col:
        description = "RT prediction accuracy across retention time."
        helptext = """
            [DIA-NN: report.tsv] Average absolute deviation between observed and predicted RT.
            Lower values indicate better RT prediction accuracy.
            """
    else:
        description = f"{y_col} trends over retention time."
        helptext = f"[DIA-NN: report.tsv] Average {y_col} across RT bins for each run."

    add_sub_section(
        sub_section=sub_section,
        plot=line_html,
        order=order,
        description=description,
        helptext=helptext
    )

# Intensity Distribution
def draw_dia_intensity_dis(sub_section, df):

    box_data = {
        str(run): group["log_intensity"].dropna().tolist()
        for run, group in df.groupby("Run")
    }

    # Determine which intensity column was used
    intensity_col = df["intensity_column"].iloc[0] if "intensity_column" in df.columns else "Precursor.Normalised"
    
    draw_config = {
        "id": "intensity_distribution_box",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Intensity Distribution",
        "tt_decimals": 5,
        "xlab": f"log2({intensity_col})",
    }

    box_html = box.plot(
        list_of_data_by_sample=box_data,
        pconfig=draw_config
    )

    box_html = remove_subtitle(box_html)

    # Dynamic description based on column used
    if intensity_col == "Precursor.Quantity":
        description_text = "log2(Precursor.Quantity) for each Run (non-normalized intensity)."
        helptext_text = """
            [DIA-NN: report.tsv] log2(Precursor.Quantity) for each Run. Shows the non-normalized 
            intensity values as recommended. Zero values are filtered out before log-transform.
            """
    elif intensity_col == "Ms1.Area":
        description_text = "log2(Ms1.Area) for each Run."
        helptext_text = """
            [DIA-NN: report.tsv] log2(Ms1.Area) for each Run. MS1 area values are used when 
            Precursor.Quantity is not available. Zero values are filtered out before log-transform.
            """
    else:
        description_text = "log2(Precursor.Normalised) for each Run (fallback)."
        helptext_text = """
            [DIA-NN: report.tsv] log2(Precursor.Normalised) for each Run. Using normalized values 
            as fallback when non-normalized quantities are not available.
            """

    add_sub_section(
        sub_section=sub_section,
        plot=box_html,
        order=3,
        description=description_text,
        helptext=helptext_text
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

def draw_dia_missed_cleavages(sub_section, df):
    """Draw missed tryptic cleavages distribution per run"""
    
    # We need peptide sequence information to calculate missed cleavages
    sequence_col = None
    if "Stripped.Sequence" in df.columns:
        sequence_col = "Stripped.Sequence"
    elif "Sequence" in df.columns:
        sequence_col = "Sequence"
    elif "Modified.Sequence" in df.columns:
        # Strip modifications from Modified.Sequence
        df = df.copy()
        pattern = re.compile(r"\(.*?\)")
        df["clean_sequence"] = df["Modified.Sequence"].apply(lambda x: re.sub(pattern, "", x))
        sequence_col = "clean_sequence"
    
    if sequence_col is None:
        return  # No sequence information available
    
    def count_missed_cleavages(sequence):
        """Count missed cleavages in a peptide sequence (tryptic cleavages)"""
        if pd.isna(sequence) or len(sequence) < 2:
            return 0
        
        # Trypsin cleaves after K and R, except when followed by P
        missed = 0
        for i in range(len(sequence) - 1):
            if sequence[i] in ['K', 'R'] and sequence[i + 1] != 'P':
                missed += 1
        
        return missed
    
    # Calculate missed cleavages for each peptide
    df_mc = df.copy()
    df_mc["missed_cleavages"] = df_mc[sequence_col].apply(count_missed_cleavages)
    
    # Create bar chart data similar to charge distribution
    bar_data = (
        df_mc.groupby("Run")["missed_cleavages"]
        .value_counts()
        .sort_index()
        .unstack(fill_value=0)
        .to_dict(orient="index")
    )
    
    bar_data = {str(k): v for k, v in bar_data.items()}
    
    if not bar_data:
        return  # No data to plot
    
    draw_config = {
        "id": "missed_cleavages_distribution",
        "cpswitch": True,
        "title": "Missed Tryptic Cleavages Distribution",
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
        order=7,
        description="Distribution of missed tryptic cleavages per run.",
        helptext="""
            [DIA-NN: report.tsv] Distribution of missed tryptic cleavages in identified peptides.
            Trypsin cleaves after K and R residues (except when followed by P). High numbers of 
            missed cleavages may indicate incomplete digestion or non-specific cleavage.
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

    # Determine which intensity column was used
    intensity_col = df["intensity_column"].iloc[0] if "intensity_column" in df.columns else "Precursor.Normalised"

    draw_box_config = {
        "id": "dia_std_intensity_box",
        "title": "Standard Deviation of Intensity",
        "cpswitch": False,
        "tt_decimals": 5,
        "xlab": f"Standard Deviation of log2({intensity_col})",
    }
    
    box_html = box.plot(
        list_of_data_by_sample=box_data,
        pconfig=draw_box_config,
    )

    box_html = remove_subtitle(box_html)

    # Dynamic description based on column used
    if intensity_col == "Precursor.Quantity":
        description_text = "Standard deviation of non-normalized intensity under different experimental conditions."
        helptext_text = """
            [DIA-NN: report.tsv] First, identify the experimental conditions from the "Run" name. 
            Then, group the data by experimental condition and Modified.Sequence, and calculate 
            the standard deviation of log2(Precursor.Quantity).
            """
    elif intensity_col == "Ms1.Area":
        description_text = "Standard deviation of MS1 area under different experimental conditions."
        helptext_text = """
            [DIA-NN: report.tsv] First, identify the experimental conditions from the "Run" name. 
            Then, group the data by experimental condition and Modified.Sequence, and calculate 
            the standard deviation of log2(Ms1.Area).
            """
    else:
        description_text = "Standard deviation of intensity under different experimental conditions."
        helptext_text = """
            [DIA-NN: report.tsv] First, identify the experimental conditions from the "Run" name. 
            Then, group the data by experimental condition and Modified.Sequence, and calculate 
            the standard deviation of log2(Precursor.Normalised).
            """

    add_sub_section(
        sub_section=sub_section,
        plot=box_html,
        order=6,
        description=description_text,
        helptext=helptext_text
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