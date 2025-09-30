from matplotlib.colors import LinearSegmentedColormap, to_hex 
from multiqc.plots import (
    bargraph,
    linegraph,
    heatmap,
    table
)
import itertools
import numpy as np
from multiqc.plots.table_object import InputRow
from multiqc.types import SampleGroup, SampleName
from typing import Dict, List
import re
from collections import OrderedDict
from multiqc import config

from pmultiqc.modules.core.section_groups import add_sub_section
from pmultiqc.modules.common.common_utils import read_openms_design, condition_split


# HeatMap color list
color_map = LinearSegmentedColormap.from_list("red_green", ["#ff0000", "#00ff00"])
HEATMAP_COLOR_LIST = [[s, to_hex(color_map(s))] for s in [round(i * 0.1, 1) for i in range(11)]]


# MS/MS Identified Per Raw File
def draw_ms_ms_identified(sub_section, msms_identified_percent):

    draw_config = {
        "id": "msms_identified_per_raw_file",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "MS/MS Identified Per Raw File",
        "tt_decimals": 2,
        "ylab": "MS/MS Identified [%]",
    }

    bar_html = bargraph.plot(data=msms_identified_percent, pconfig=draw_config)

    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=7,
        description="MS/MS identification rate per raw file.",
        helptext="MS/MS identification rate per raw file (quantms data from mzTab and mzML files; MaxQuant data from summary.txt)",
    )


# Potential Contaminants per Group
def draw_potential_contaminants(sub_section, contaminant_percent, is_maxquant):

    draw_config = {
        "id": "potential_contaminants_per_group",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Potential Contaminants Per File",
        "tt_decimals": 2,
        "ylab": "Percent [%]",
    }

    bar_html = bargraph.plot(data=contaminant_percent, pconfig=draw_config)

    bar_html = remove_subtitle(bar_html)

    if is_maxquant:
        description_text = "Potential contaminants per group from proteinGroups.txt."
        help_text = """
            External protein contamination should be controlled for, therefore MaxQuant ships with a 
            comprehensive, yet customizable protein contamination database, which is searched by MaxQuant 
            by default. 
            
            A contamination plot derived from the proteinGroups (PG) table, showing the 
            fraction of total protein intensity attributable to contaminants. 
            
            Note that this plot is based on experimental groups, and therefore may not correspond 1:1 to Raw files.
            """
    else:
        description_text = "Potential contaminants per file from mzTab."
        help_text = """
            A contamination plot derived from the Peptide section of the mzTab file, showing the 
            fraction of total intensity attributable to contaminants.
            """

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=2,
        description=description_text,
        helptext=help_text,
    )


# Top5 Contaminants Per Raw File
def draw_top_n_contaminants(sub_section, top_contaminants_data):

    draw_config = {
        "id": "top_contaminants_per_raw_file",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Top5 Contaminants Per Raw File",
        "tt_decimals": 2,
        "ylab": "Identified [%]",
    }

    bar_html = bargraph.plot(
        top_contaminants_data["plot_data"],
        cats=top_contaminants_data["cats"],
        pconfig=draw_config,
    )

    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=1,
        description="The five most abundant external protein contaminants by Raw file",
        helptext="""
            pmultiqc will explicitly show the five most abundant external protein contaminants 
            (as detected via MaxQuant's contaminants FASTA file) by Raw file, and summarize the 
            remaining contaminants as 'other'. This allows to track down which proteins exactly 
            contaminate your sample. Low contamination is obviously better.
            
            If you see less than 5 contaminants, it either means there are actually less, or that 
            one (or more) of the shortened contaminant names subsume multiple of the top5 
            contaminants (since they start with the same prefix).
            """,
    )


# Charge-state of Per File
def draw_charge_state(sub_section, charge_data, report_type):

    draw_config = {
        "id": "charge_state_of_per_raw_file",
        "cpswitch": True,
        "title": "Charge-state of Per File",
        "tt_decimals": 0,
        "ylab": "Count",
    }

    bar_html = bargraph.plot(
        data=charge_data["plot_data"], cats=charge_data["cats"], pconfig=draw_config
    )

    bar_html = remove_subtitle(bar_html)

    help_text = """Charge distribution per Raw file. For typtic digests, peptides of charge 2
    (one N-terminal and one at tryptic C-terminal R or K residue) should be dominant.
    Ionization issues (voltage?), in-source fragmentation, missed cleavages and buffer irregularities can
    cause a shift (see Bittremieux 2017, DOI: 10.1002/mas.21544).
    The charge distribution should be similar across Raw files.
    Consistent charge distribution is paramount for comparable 3D-peak intensities across samples."""

    if report_type == "MaxQuant":
        description_text = "The distribution of the charge-state of the precursor ion, excluding potential contaminants."
        help_text += "<p>This plot ignores charge states of contaminants.<p>"
    elif report_type == "mzIdentML":
        description_text = "The distribution of the charge-state of the precursor ion."
    else:
        description_text = "The distribution of precursor ion charge states (based on mzTab data)."
        help_text += "<p>Precursor ion charge states are based on mzTab data.<p>"

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=6,
        description=description_text,
        helptext=help_text,
    )


# Modifications
def draw_modifications(sub_section, modified_data):

    draw_config = {
        "id": "modifications_per_raw_file",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Modifications Per Raw File",
        "stacking": "group",
        "tt_decimals": 3,
        "ylab": "Occurence [%]",
    }

    bar_html = bargraph.plot(
        data=modified_data["plot_data"], cats=modified_data["cats"], pconfig=draw_config
    )

    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=6,
        description="""
            Compute an occurrence table of modifications (e.g. Oxidation (M)) for all peptides, including the unmodified (but without contaminants).
            """,
        helptext="""
Post-translational modifications contained within the identified peptide sequence.<br>

<p>The plot will show percentages, i.e. is normalized by the total number of peptide sequences (where different charge state counts as a separate peptide) per Raw file.
The sum of frequencies may exceed 100% per Raw file, since a peptide can have multiple modifications.</p>

E.g. given three peptides in a single Raw file                 <br>
1. _M(Oxidation (M))LVLDEADEM(Oxidation (M))LNK_               <br>
2. _(Acetyl (Protein N-term))M(Oxidation (M))YGLLLENLSEYIK_    <br>
3. DPFIANGER                                                   <br>

<p>, the following frequencies arise:</p>

* 33% of 'Acetyl (Protein N-term)' <br>
* 33% of 'Oxidation (M)'           <br>
* 33% of '2 Oxidation (M)'         <br>
* 33% of 'Unmodified'              <br>

<p>Thus, 33% of sequences are unmodified, implying 66% are modified at least once. 
If a modification, e.g. Oxidation(M), occurs multiple times in a single peptide it's listed as a separate modification (e.g. '2 Oxidation (M)' for double oxidation of a single peptide).</p>
""",
    )


def draw_oversampling(sub_section, oversampling, oversampling_plot, is_maxquant):

    if is_maxquant:
        draw_config = {
            "id": "oversampling_distribution",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "MS/MS Counts Per 3D-peak",
            "tt_decimals": 2,
            "ylab": "MS/MS Counts Per 3D-peak [%]",
        }

        bar_html = bargraph.plot(
            data=oversampling["plot_data"], cats=oversampling["cats"], pconfig=draw_config
        )

    else:
        draw_config = {
            "id": "oversampling_distribution",
            "cpswitch": True,
            "cpswitch_c_active": False,
            "title": "MS/MS Counts Per 3D-peak",
            "ylab": "MS/MS Counts Per 3D-peak [%]",
            "tt_decimals": 0,
        }

        bar_html = bargraph.plot(data=oversampling, cats=oversampling_plot, pconfig=draw_config)

    bar_html = remove_subtitle(bar_html)

    helptext = """
                For high complexity samples, oversampling of individual 3D-peaks automatically leads to 
                undersampling or even omission of other 3D-peaks, reducing the number of identified peptides. 
                Oversampling occurs in low-complexity samples or long LC gradients, as well as undersized dynamic 
                exclusion windows for data independent acquisitions.
                """
    if is_maxquant:
        helptext += "<p>If DIA-Data: this metric is skipped.</p>"

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=7,
        description="""
            An oversampled 3D-peak is defined as a peak whose peptide ion 
            (same sequence and same charge state) was identified by at least two distinct MS2 spectra 
            in the same Raw file.
            """,
        helptext=helptext,
    )


# Missed Cleavages Per Raw File
def draw_msms_missed_cleavages(sub_section, missed_cleavages_data, is_maxquant):

    draw_config = {
        "id": "missed_cleavages_per_raw_file",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Missed Cleavages Per Raw File",
        "tt_decimals": 2,
        "ylab": "Missed Cleavages [%]",
    }

    bar_html = bargraph.plot(
        data=missed_cleavages_data["plot_data"],
        cats=missed_cleavages_data["cats"],
        pconfig=draw_config,
    )

    bar_html = remove_subtitle(bar_html)

    helptext = """
                Under optimal digestion conditions (high enzyme grade etc.), only few missed cleavages (MC) are expected. In 
                general, increased MC counts also increase the number of peptide signals, thus cluttering the available 
                space and potentially provoking overlapping peptide signals, biasing peptide quantification.
                Thus, low MC counts should be favored. Interestingly, it has been shown recently that 
                incorporation of peptides with missed cleavages does not negatively influence protein quantification (see 
                [Chiva, C., Ortega, M., and Sabido, E. Influence of the Digestion Technique, Protease, and Missed 
                Cleavage Peptides in Protein Quantitation. J. Proteome Res. 2014, 13, 3979-86](https://doi.org/10.1021/pr500294d) ). 
                However this is true only if all samples show the same degree of digestion. High missed cleavage values 
                can indicate for example, either a) failed digestion, b) a high (post-digestion) protein contamination, or 
                c) a sample with high amounts of unspecifically degraded peptides which are not digested by trypsin. 

                If MC>=1 is high (>20%) you should re-analyse with increased missed cleavages parameters and compare the number of peptides.
                Usually high MC correlates with bad identification rates, since many spectra cannot be matched to the forward database.
                """

    if is_maxquant:
        description_text = "[Excludes Contaminants] Missed Cleavages per raw file."
        helptext += "<p>In the rare case that 'no enzyme' was specified in MaxQuant, neither scores nor plots are shown.</p>"
    else:
        description_text = "Missed Cleavages per raw file."

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=5,
        description=description_text,
        helptext=helptext,
    )


# IDs over RT
def draw_ids_rt_count(sub_section, rt_count_data, report_type):

    draw_config = {
        "id": "IDs_over_RT",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "IDs over RT",
        "ymin": 0,
        "tt_decimals": 3,
        "ylab": "Count",
        "xlab": "Retention time [min]",
        "showlegend": True,
    }

    linegraph_html = linegraph.plot(data=rt_count_data, pconfig=draw_config)

    linegraph_html = remove_subtitle(linegraph_html)

    if report_type == "maxquant":
        description_text = "Distribution of retention time, derived from the evidence table."
        help_text = """
            The uncalibrated retention time in minutes in the elution profile of the precursor ion, 
            and does not include potential contaminants.
            """
    elif report_type == "dia":
        description_text = "Distribution of retention time, derived from the main report."
        help_text = """
            [DIA-NN: main report] Distribution of retention time (RT) for each run.
            """
    elif report_type == "mzIdentML":
        description_text = "Distribution of retention time, derived from the mzIdentML (or mzML)."
        help_text = """
            The uncalibrated retention time in minutes in the elution profile of the precursor ion.
            """
    else:
        description_text = "Distribution of retention time, derived from the mzTab."
        help_text = """
            The uncalibrated retention time in minutes in the elution profile of the precursor ion.
            """

    help_text += """
            <p>This plot allows to judge column occupancy over retention time.
            Ideally, the LC gradient is chosen such that the number of identifications (here, after FDR filtering) is
            uniform over time, to ensure consistent instrument duty cycles. Sharp peaks and uneven distribution of
            identifications over time indicate potential for LC gradient optimization.
            See [Moruz 2014, DOI: 10.1002/pmic.201400036](https://pubmed.ncbi.nlm.nih.gov/24700534/) for details.</p>
"""

    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=1,
        description=description_text,
        helptext=help_text,
    )


# MaxQuant: Delta Mass [Da], Delta Mass [ppm]
# quantms: Delta Mass [ppm]
def draw_delta_mass_da_ppm(sub_section, delta_mass, delta_mass_type):

    # MaxQuant: Delta Mass [Da]
    if delta_mass_type == "Mass Error [Da]":
        plot_id = "delta_mass_da"
        plot_title = "Delta Mass [Da]"
        plot_xlab = "Experimental m/z - Theoretical m/z"
        sub_section_order = 1
        description_text = """
            This plot is based on the "Mass Error [Da]" column from the evidence.txt generated by MaxQuant.
            """
        help_text = """
            Mass error of the recalibrated mass-over-charge value of the precursor ion in comparison to the
            predicted monoisotopic mass of the identified peptide sequence in milli-Dalton.
            """

    # MaxQuant: Delta Mass [ppm]
    elif delta_mass_type == "Mass Error [ppm]":
        plot_id = "delta_mass_ppm"
        plot_title = "Delta Mass [ppm]"
        plot_xlab = "Delta Mass [ppm]"
        sub_section_order = 2
        description_text = """
            This plot is based on the "Mass Error [ppm]" column from the evidence.txt generated by MaxQuant.
            """
        help_text = """
            Mass error of the recalibrated mass-over-charge value of the precursor ion in comparison to the 
            predicted monoisotopic mass of the identified peptide sequence in parts per million.
            
            Ppm errors should be centered on zero and their spread is expected to be significantly smaller than before calibration.
            """

    # quantms: Delta Mass [ppm]
    elif delta_mass_type == "quantms_ppm":
        plot_id = "delta_mass_ppm"
        plot_title = "Delta Mass [ppm]"
        plot_xlab = "Delta Mass [ppm]"
        sub_section_order = 2
        description_text = """
            Delta Mass [ppm] calculated from mzTab.
            """
        help_text = """
            Delta Mass [ppm] calculated from mzTab: ((𝑒𝑥𝑝𝑒𝑟𝑖𝑚𝑒𝑛𝑡𝑎𝑙 𝑚/𝑧 - 𝑡ℎ𝑒𝑜𝑟𝑒𝑡𝑖𝑐𝑎𝑙 𝑚/𝑧) / (𝑡ℎ𝑒𝑜𝑟𝑒𝑡𝑖𝑐𝑎𝑙 𝑚/𝑧)) × 10^6.
            """

    x_values = list(delta_mass["count"].keys())

    range_threshold = 10
    if max(abs(x) for x in x_values) > range_threshold:
        range_abs = range_threshold
    else:
        range_abs = 1
    range_step = (max(x_values) - min(x_values)) * 0.05

    if max(abs(x) for x in x_values) > range_abs:

        delta_mass_range = {k: v for k, v in delta_mass["count"].items() if abs(k) <= range_abs}
        delta_mass_percent_range = {
            k: v for k, v in delta_mass["frequency"].items() if abs(k) <= range_abs
        }

        x_values_adj = list(delta_mass_range.keys())
        range_step_adj = (max(x_values_adj) - min(x_values_adj)) * 0.05

        data_label = [
            {
                "name": f"Count (range: -{range_abs} to {range_abs})",
                "ylab": "Count",
                "tt_label": "{point.x} Mass delta counts: {point.y}",
                "xmax": max(abs(x) for x in x_values_adj) + range_step_adj,
                "xmin": -(max(abs(x) for x in x_values_adj) + range_step_adj),
                "ymin": 0,
            },
            {
                "name": f"Relative Frequency (range: -{range_abs} to {range_abs})",
                "ylab": "Relative Frequency",
                "tt_label": "{point.x} Mass delta relative frequency: {point.y}",
                "xmax": max(abs(x) for x in x_values_adj) + range_step_adj,
                "xmin": -(max(abs(x) for x in x_values_adj) + range_step_adj),
                "ymin": 0,
            },
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
            "id": plot_id,
            "title": plot_title,
            "colors": {"count": "#b2df8a", "relative_frequency": "#b2df8a"},
            "xlab": plot_xlab,
            "data_labels": data_label,
            "style": "lines+markers",
        }

        line_html = linegraph.plot(
            [
                {"count": delta_mass_range},
                {"relative_frequency": delta_mass_percent_range},
                {"count": delta_mass["count"]},
                {"relative_frequency": delta_mass["frequency"]},
            ],
            pconfig,
        )

    else:

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
            "id": plot_id,
            "title": plot_title,
            "colors": {"count": "#b2df8a", "relative_frequency": "#b2df8a"},
            "xlab": plot_xlab,
            "data_labels": data_label,
            "style": "lines+markers",
        }

        line_html = linegraph.plot(
            [{"count": delta_mass["count"]}, {"relative_frequency": delta_mass["frequency"]}],
            pconfig,
        )

    add_sub_section(
        sub_section=sub_section,
        plot=line_html,
        order=sub_section_order,
        description=description_text,
        helptext=help_text,
    )


def draw_heatmap(
    sub_sections, hm_colors, heatmap_data, heatmap_xnames, heatmap_ynames, is_maxquant
):

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
    }

    if is_maxquant:

        hm_html = heatmap.plot(data=heatmap_data, pconfig=pconfig)

        description_text = "This heatmap provides an overview of the performance of MaxQuant."

    else:

        hm_html = heatmap.plot(heatmap_data, heatmap_xnames, heatmap_ynames, pconfig)

        description_text = "This heatmap provides an overview of the performance of quantms."

    hm_html = remove_subtitle(hm_html)

    add_sub_section(
        sub_section=sub_sections,
        plot=hm_html,
        order=2,
        description=description_text,
        helptext="""
            This plot shows the pipeline performance overview. Some metrics are calculated           
            *Heatmap score[Contaminants]: as fraction of summed intensity with 0 = sample full of contaminants; 
                1 = no contaminants
            *Heatmap score[Pep Intensity (>23.0)]: Linear scale of the median intensity reaching the threshold, 
                i.e. reaching 2^21 of 2^23 gives score 0.25
            *Heatmap score[Charge]: Deviation of the charge 2 proportion from a representative 
                Raw file (median). For typtic digests, peptides of charge 2 (one N-terminal and one at 
                tryptic C-terminal R or K residue) should be dominant. Ionization issues (voltage?), 
                in-source fragmentation, missed cleavages and buffer irregularities can cause a shift 
                (see Bittremieux 2017, DOI: 10.1002/mas.21544)
            *Heatmap score [Missed Cleavages]: the fraction (0% - 100%) of fully cleaved peptides per Raw file
            *Heatmap score [Missed Cleavages Var]: each Raw file is scored for its deviation from the 'average' digestion 
                state of the current study
            *Heatmap score [ID rate over RT]: Judge column occupancy over retention time. 
                Ideally, the LC gradient is chosen such that the number of identifications 
                (here, after FDR filtering) is uniform over time, to ensure consistent instrument duty cycles. 
                Sharp peaks and uneven distribution of identifications over time indicate potential for LC gradient 
                optimization.Scored using 'Uniform' scoring function. i.e. constant receives good score, extreme shapes are bad
            *Heatmap score [MS2 Oversampling]: The percentage of non-oversampled 3D-peaks. An oversampled 
                3D-peak is defined as a peak whose peptide ion (same sequence and same charge state) was 
                identified by at least two distinct MS2 spectra in the same Raw file. For high complexity samples, 
                oversampling of individual 3D-peaks automatically leads to undersampling or even omission of other 3D-peaks, 
                reducing the number of identified peptides
            *Heatmap score [Pep Missing Values]: Linear scale of the fraction of missing peptides
        """,
    )


def remove_subtitle(plot_html):

    for dataset in plot_html.datasets:

        if "subtitle" in dataset.dconfig:
            dataset.dconfig["subtitle"] = ""

        title_text = ""
        if dataset.layout and "title" in dataset.layout:
            title_text = dataset.layout["title"].get("text", "")
            if "<br><sup>" in title_text:
                dataset.layout["title"]["text"] = title_text.split("<br><sup>")[0]

    return plot_html

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

# self.sub_sections["experiment"]
def draw_exp_design(sub_sections, exp_design):
    # Currently this only supports the OpenMS two-table format (default in quantms pipeline)
    # One table format would actually be even easier. You can just use pandas.read_tsv
    sample_df, file_df = read_openms_design(exp_design)

    exp_design_runs = np.unique(file_df["Run"].tolist())

    is_bruker = False
    if file_df["Spectra_Filepath"][0].endswith((".d", ".d.tar")):
        is_bruker = True

    # Create table plot
    pattern = r'^(\w+=[^=;]+)(;\w+=[^=;]+)*$'
    is_multi_conditions = all(sample_df["MSstats_Condition"].apply(lambda x: bool(re.match(pattern, str(x)))))

    rows_by_group: Dict[SampleGroup, List[InputRow]] = {}

    if is_multi_conditions:
        for sample in sorted(sample_df["Sample"].tolist(), key=lambda x: int(x)):
            file_df_sample = file_df[file_df["Sample"] == sample].copy()
            sample_df_slice = sample_df[sample_df["Sample"] == sample].copy()
            row_data: List[InputRow] = []

            sample_data = {}
            for k, v in condition_split(sample_df_slice["MSstats_Condition"].iloc[0]).items():
                sample_data["MSstats_Condition_" + str(k)] = v
            
            sample_data["MSstats_BioReplicate"] = sample_df_slice["MSstats_BioReplicate"].iloc[0]
            sample_data["Fraction_Group"] = ""
            sample_data["Fraction"] = ""
            sample_data["Label"] = ""

            row_data.append(
                InputRow(
                    sample=SampleName(sample),
                    data=sample_data,
                )
            )

            for _, row in file_df_sample.iterrows():

                sample_data = {}
                for k, _ in condition_split(sample_df_slice["MSstats_Condition"].iloc[0]).items():
                    sample_data["MSstats_Condition_" + str(k)] = ""

                sample_data["MSstats_BioReplicate"] = ""
                sample_data["Fraction_Group"] = row["Fraction_Group"]
                sample_data["Fraction"] = row["Fraction"]
                sample_data["Label"] = row["Label"]

                row_data.append(
                    InputRow(
                        sample=SampleName(row["Run"]),
                        data=sample_data,
                    )
                )
            group_name: SampleGroup = SampleGroup(sample)
            rows_by_group[group_name] = row_data
        headers = {}
        headers["Sample"] = {
            "title": "Sample [Spectra File]",
            "description": "",
            "scale": False,
        }
        for k, _ in condition_split(sample_df_slice["MSstats_Condition"].iloc[0]).items():
            headers["MSstats_Condition_" + str(k)] ={
                "title": "MSstats Condition: " + str(k),
                "description": "",
                "scale": False,
            }
        headers["MSstats_BioReplicate"] = {
            "title": "MSstats BioReplicate",
            "description": "",
            "scale": False,
        }
        headers["Fraction_Group"] = {
            "title": "Fraction Group",
            "description": "",
            "scale": False,
        }
        headers["Fraction"] = {
            "title": "Fraction",
            "description": "Fraction Identifier",
            "scale": False,
        }
        headers["Label"] = {
            "title": "Label",
            "description": "",
            "scale": False,
        }
    else:
        for sample in sorted(sample_df["Sample"].tolist(), key=lambda x: int(x)):
            file_df_sample = file_df[file_df["Sample"] == sample].copy()
            sample_df_slice = sample_df[sample_df["Sample"] == sample].copy()
            row_data: List[InputRow] = []
            row_data.append(
                InputRow(
                    sample=SampleName(sample),
                    data={
                        "MSstats_Condition": sample_df_slice["MSstats_Condition"].iloc[0],
                        "MSstats_BioReplicate": sample_df_slice["MSstats_BioReplicate"].iloc[0],
                        "Fraction_Group": "",
                        "Fraction": "",
                        "Label": "",
                    },
                )
            )
            for _, row in file_df_sample.iterrows():
                row_data.append(
                    InputRow(
                        sample=SampleName(row["Run"]),
                        data={
                            "MSstats_Condition": "",
                            "MSstats_BioReplicate": "",
                            "Fraction_Group": row["Fraction_Group"],
                            "Fraction": row["Fraction"],
                            "Label": row["Label"],
                        },
                    )
                )
            group_name: SampleGroup = SampleGroup(sample)
            rows_by_group[group_name] = row_data

        headers = {
            "Sample": {
                "title": "Sample [Spectra File]",
                "description": "",
                "scale": False,
            },
            "MSstats_Condition": {
                "title": "MSstats Condition",
                "description": "MSstats Condition",
                "scale": False,
            },
            "MSstats_BioReplicate": {
                "title": "MSstats BioReplicate",
                "description": "MSstats BioReplicate",
                "scale": False,
            },
            "Fraction_Group": {
                "title": "Fraction Group",
                "description": "Fraction Group",
                "scale": False,
            },
            "Fraction": {
                "title": "Fraction",
                "description": "Fraction Identifier",
                "scale": False,
            },
            "Label": {
                "title": "Label",
                "description": "Label",
                "scale": False,
            },
        }

    pconfig = {
        "id": "experimental_design",
        "title": "Experimental Design",
        "save_file": False,
        "raw_data_fn": "multiqc_Experimental_Design_table",
        "no_violin": True,
    }
    table_html = table.plot(rows_by_group, headers, pconfig)
    add_sub_section(
        sub_section=sub_sections,
        plot=table_html,
        order=1,
        description="""
            This table shows the design of the experiment. I.e., which files and channels correspond to which sample/condition/fraction.
            """,
        helptext="""
            You can see details about it in 
            https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/classOpenMS_1_1ExperimentalDesign.html
            """
    )

    return sample_df, file_df, exp_design_runs, is_bruker, is_multi_conditions

# self.sub_sections["summary"]
def draw_summary_protein_ident_table(
        sub_sections,
        enable_dia: bool = False,
        total_peptide_count: int = 0,
        total_protein_quantified: int = 0,
        total_ms2_spectra: int = 0,
        total_ms2_spectral_identified: int = 0,
        total_protein_identified: int = 0,
    ):
    headers = OrderedDict()
    if enable_dia:
        summary_table = {
            total_peptide_count: {"#Proteins Quantified": total_protein_quantified}
        }
        col_header = "#Peptides Quantified"
    else:
        summary_table = {
            total_ms2_spectra: {
                "#Identified MS2 Spectra": total_ms2_spectral_identified
            }
        }
        coverage = total_ms2_spectral_identified / total_ms2_spectra * 100
        summary_table[total_ms2_spectra]["%Identified MS2 Spectra"] = coverage
        summary_table[total_ms2_spectra][
            "#Peptides Identified"
        ] = total_peptide_count
        summary_table[total_ms2_spectra][
            "#Proteins Identified"
        ] = total_protein_identified

        if not config.kwargs.get("mzid_plugin", False):
            summary_table[total_ms2_spectra][
                "#Proteins Quantified"
            ] = total_protein_quantified

        headers["#Identified MS2 Spectra"] = {
            "description": "Total number of MS/MS spectra identified",
        }
        headers["%Identified MS2 Spectra"] = {
            "description": "Percentage of Identified MS/MS Spectra",
            "format": "{:,.2f}",
            "suffix": "%",
        }
        col_header = "#MS2 Spectra"

    # Create table plot
    pconfig = {
        "id": "identification_summary_table",  # ID used for the table
        "title": "Summary Table",  # Title of the table. Used in the column config modal
        "save_file": False,  # Whether to save the table data to a file
        "raw_data_fn": "multiqc_summary_table_table",  # File basename to use for raw data file
        "sort_rows": False,  # Whether to sort rows alphabetically
        "only_defined_headers": False,  # Only show columns that are defined in the headers config
        "col1_header": col_header,
        # 'format': '{:,.0f}',  # The header used for the first column
        "scale": "Set1",
    }

    table_html = table.plot(summary_table, headers, pconfig)

    if config.kwargs.get("mzid_plugin", False) or config.kwargs.get("parse_diann", False):
        description_str = "This plot shows the summary statistics of the submitted data."
        # TODO: add description here @Yasset
        helptext_str = """
            This plot shows the summary statistics of the submitted data.
            """
    else:
        description_str = "This table shows the quantms pipeline summary statistics."
        # TODO: add description here @Yasset
        helptext_str = """
            This table shows the quantms pipeline summary statistics.
            """

    add_sub_section(
        sub_section=sub_sections,
        plot=table_html,
        order=1,
        description=description_str,
        helptext=helptext_str
    )

# self.sub_sections["summary"]
def draw_quantms_identi_num(
        sub_sections,
        enable_exp,
        enable_sdrf,
        is_multi_conditions,
        sample_df,
        file_df,
        cal_num_table_data
    ):

    # Create table data
    rows_by_group: Dict[SampleGroup, List[InputRow]] = {}

    if enable_exp or enable_sdrf:
            
        if is_multi_conditions:
            for sample in sorted(sample_df["Sample"].tolist(), key=lambda x: int(x)):
                file_df_sample = file_df[file_df["Sample"] == sample].copy()
                sample_df_slice = sample_df[sample_df["Sample"] == sample].copy()
                row_data: List[InputRow] = []

                sample_data = {}
                for k, v in condition_split(sample_df_slice["MSstats_Condition"].iloc[0]).items():
                    sample_data["MSstats_Condition_" + str(k)] = v
                
                sample_data["MSstats_BioReplicate"] = sample_df_slice["MSstats_BioReplicate"].iloc[0]
                sample_data["Fraction"] = ""
                sample_data["Peptide_Num"] = ""
                sample_data["Unique_Peptide_Num"] = ""
                sample_data["Modified_Peptide_Num"] = ""
                sample_data["Protein_Num"] = ""

                row_data.append(
                    InputRow(
                        sample=SampleName(sample),
                        data=sample_data,
                    )
                )

                for _, row in file_df_sample.iterrows():

                    sample_data = {}
                    for k, _ in condition_split(sample_df_slice["MSstats_Condition"].iloc[0]).items():
                        sample_data["MSstats_Condition_" + str(k)] = ""

                    sample_data["Fraction"] = row["Fraction"]
                    sample_data["Peptide_Num"] = cal_num_table_data[row["Run"]]["peptide_num"]
                    sample_data["Unique_Peptide_Num"] = cal_num_table_data[row["Run"]]["unique_peptide_num"]
                    sample_data["Modified_Peptide_Num"] = cal_num_table_data[row["Run"]]["modified_peptide_num"]
                    sample_data["Protein_Num"] = cal_num_table_data[row["Run"]]["protein_num"]

                    row_data.append(
                        InputRow(
                            sample=SampleName(row["Run"]),
                            data=sample_data,
                        )
                    ) 
                group_name: SampleGroup = SampleGroup(sample)
                rows_by_group[group_name] = row_data
                
            headers = {}
            for k, _ in condition_split(sample_df_slice["MSstats_Condition"].iloc[0]).items():
                headers["MSstats_Condition_" + str(k)] ={
                    "title": "MSstats Condition: " + str(k),
                    "description": "",
                    "scale": False,
                }
            headers["Fraction"] = {
                "title": "Fraction",
                "description": "Fraction Identifier",
                "scale": False,
            }
            headers["Peptide_Num"] = {
                "title": "#Peptide IDs",
                "description": "The number of identified PSMs in the pipeline",
            }
            headers["Unique_Peptide_Num"] = {
                "title": "#Unambiguous Peptide IDs",
                "description": "The number of unique peptides in the pipeline. Those that match only one protein in the provided database",
            }
            headers["Modified_Peptide_Num"] = {
                "title": "#Modified Peptide IDs",
                "description": "Number of modified identified peptides in the pipeline",
            }
            headers["Protein_Num"] = {
                "title": "#Protein (group) IDs",
                "description": "The number of identified protein(group)s in the pipeline",
            }
        else:
            for sample in sorted(sample_df["Sample"].tolist(), key=lambda x: int(x)):
                file_df_sample = file_df[file_df["Sample"] == sample].copy()
                sample_df_slice = sample_df[sample_df["Sample"] == sample].copy()
                row_data: List[InputRow] = []
                row_data.append(
                    InputRow(
                        sample=SampleName(sample),
                        data={
                            "MSstats_Condition": sample_df_slice["MSstats_Condition"].iloc[0],
                            "Fraction": "",
                            "Peptide_Num": "",
                            "Unique_Peptide_Num": "",
                            "Modified_Peptide_Num": "",
                            "Protein_Num": "",
                        },
                    )
                )
                for _, row in file_df_sample.iterrows():
                    row_data.append(
                        InputRow(
                            sample=SampleName(row["Run"]),
                            data={
                                "MSstats_Condition": "",
                                "Fraction": row["Fraction"],
                                "Peptide_Num": cal_num_table_data[row["Run"]]["peptide_num"],
                                "Unique_Peptide_Num": cal_num_table_data[row["Run"]]["unique_peptide_num"],
                                "Modified_Peptide_Num": cal_num_table_data[row["Run"]]["modified_peptide_num"],
                                "Protein_Num": cal_num_table_data[row["Run"]]["protein_num"],
                                },
                            )
                        )
                group_name: SampleGroup = SampleGroup(sample)
                rows_by_group[group_name] = row_data

            headers = {
                "MSstats_Condition": {
                    "title": "MSstats_Condition",
                    "description": "MSstats Condition",
                    "scale": False,
                },
                "Fraction": {
                    "title": "Fraction",
                    "description": "Fraction Identifier",
                    "scale": False,
                },
                "Peptide_Num": {
                    "title": "#Peptide IDs",
                    "description": "The number of identified PSMs in the pipeline",
                },
                "Unique_Peptide_Num": {
                    "title": "#Unambiguous Peptide IDs",
                    "description": "The number of unique peptides in the pipeline. Those that match only one protein in the provided database",
                },
                "Modified_Peptide_Num": {
                    "title": "#Modified Peptide IDs",
                    "description": "Number of modified identified peptides in the pipeline",
                },
                "Protein_Num": {
                    "title": "#Protein (group) IDs",
                    "description": "The number of identified protein(group)s in the pipeline",
                },
            }

    else:
        rows_by_group = dict()
        for sample, value in cal_num_table_data.items():
            rows_by_group[sample] = {
                "Peptide_Num": value["peptide_num"],
                "Unique_Peptide_Num": value["unique_peptide_num"],
                "Modified_Peptide_Num": value["modified_peptide_num"],
                "Protein_Num": value["protein_num"],
            }
        
        headers = {
            "Peptide_Num": {
                "title": "#Peptide IDs",
                "description": "The number of identified PSMs in the pipeline",
            },
            "Unique_Peptide_Num": {
                "title": "#Unambiguous Peptide IDs",
                "description": "The number of unique peptides in the pipeline. Those that match only one protein in the provided database",
            },
            "Modified_Peptide_Num": {
                "title": "#Modified Peptide IDs",
                "description": "Number of modified identified peptides in the pipeline",
            },
            "Protein_Num": {
                "title": "#Protein (group) IDs",
                "description": "The number of identified protein(group)s in the pipeline",
            },
        }

    # Create table plot
    pconfig = {
        "id": "pipeline_result_statistics",
        "title": "Pipeline Result Statistics",
        "save_file": False,
        "raw_data_fn": "multiqc_pipeline_result_statistics",
        "no_violin": True,
    }
    table_html = table.plot(rows_by_group, headers, pconfig)
    add_sub_section(
        sub_section=sub_sections,
        plot=table_html,
        order=3,
        description="This plot shows the quantms pipeline final result.",
        helptext="""
            Including Sample Name, Possible Study Variables, identified the number of peptide in the pipeline,
            and identified the number of modified peptide in the pipeline, eg. All data in this table are obtained 
            from the out_msstats file. You can also remove the decoy with the `remove_decoy` parameter.
            """
    )

# self.sub_sections["identification"]
def draw_num_pep_per_protein(
        sub_sections,
        pep_plot
    ):

    if any([len(i) >= 100 for i in pep_plot.dict["data"].values()]):
        data_labels = ["Frequency", "Percentage"]
    else:
        data_labels = [
                {
                    "name": "Frequency",
                    "ylab": "Frequency",
                    "tt_suffix": "",
                    "tt_decimals": 0,
                },
                {
                    "name": "Percentage",
                    "ylab": "Percentage [%]",
                    "tt_suffix": "%",
                    "tt_decimals": 2,
                },
            ]

    pconfig = {
        "id": "number_of_peptides_per_proteins",
        "cpswitch": False,
        "title": "Number of Peptides identified Per Protein",
        "data_labels": data_labels,
    }

    bar_html = bargraph.plot(
        [pep_plot.dict["data"]["frequency"], pep_plot.dict["data"]["percentage"]],
        ["Frequency", "Percentage"],
        pconfig,
    )
    bar_html = remove_subtitle(bar_html)

    if config.kwargs.get("mzid_plugin", False):
        description_str = (
            "This plot shows the number of peptides per protein in the submitted data"
        )
        helptext_str = """
                    Proteins supported by more peptide identifications can constitute more confident results.
                """
    else:
        description_str = "This plot shows the number of peptides per protein in quantms pipeline final result"
        helptext_str = """
                    This statistic is extracted from the out_msstats file. Proteins supported by more peptide 
                    identifications can constitute more confident results.
                """
    
    add_sub_section(
        sub_section=sub_sections,
        plot=bar_html,
        order=1,
        description=description_str,
        helptext=helptext_str
    )

# self.sub_sections["ms2"]
def draw_peaks_per_ms2(
        sub_sections,
        peaks_ms2_plot,
        ms_info
    ):

    pconfig = {
        "id": "peaks_per_ms2",
        "cpswitch": False,
        "title": "Number of Peaks per MS/MS spectrum",
        "stacking": "group",
        "logswitch": True,
        "logswitch_active": True,
        "logswitch_label": "Log10 Scale",
        "ylab": "Count",
        "tt_decimals": 0,
    }
    # if config.kwargs.get("mzid_plugin", False) and self.mgf_paths:
    #     cats = self.mgf_peaks_ms2_plot.dict["cats"]
    # else:
    #     cats = self.mzml_peaks_ms2_plot.dict["cats"]

    cats = peaks_ms2_plot.dict["cats"]

    bar_html = bargraph.plot(ms_info["peaks_per_ms2"], cats, pconfig)
    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_sections,
        plot=bar_html,
        order=1,
        description="""
            This chart represents a histogram containing the number of peaks per MS/MS spectrum 
            in a given experiment.
            """,
        helptext="""
            This chart assumes centroid data. Too few peaks can identify poor fragmentation 
                or a detector fault, as opposed to a large number of peaks representing very noisy spectra. 
                This chart is extensively dependent on the pre-processing steps performed to the spectra 
                (centroiding, deconvolution, peak picking approach, etc).
                """,
        )

# self.sub_sections["ms2"]
def draw_peak_intensity_distribution(
        sub_sections,
        peak_distribution_plot,
        ms_info
    ):

    pconfig = {
        "id": "peak_intensity_distribution",
        "title": "Peak Intensity Distribution",
        "cpswitch": False,
        "stacking": "group",
        "logswitch": True,
        "logswitch_active": True,
        "logswitch_label": "Log10 Scale",
        "ylab": "Count",
        "tt_decimals": 0,
    }
    # if config.kwargs.get("mzid_plugin", False) and self.mgf_paths:
    #     cats = self.mgf_peak_distribution_plot.dict["cats"]
    # else:
    #     cats = self.mzml_peak_distribution_plot.dict["cats"]
    cats = peak_distribution_plot.dict["cats"]

    bar_html = bargraph.plot(ms_info["peak_distribution"], cats, pconfig)
    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_sections,
        plot=bar_html,
        order=2,
        description="""
            This is a histogram representing the ion intensity vs.
            the frequency for all MS2 spectra in a whole given experiment.
            It is possible to filter the information for all, identified and unidentified spectra.
            This plot can give a general estimation of the noise level of the spectra.
            """,
        helptext="""
            Generally, one should expect to have a high number of low intensity noise peaks with a low number 
            of high intensity signal peaks. 
            A disproportionate number of high signal peaks may indicate heavy spectrum pre-filtering or 
            potential experimental problems. In the case of data reuse this plot can be useful in 
            identifying the requirement for pre-processing of the spectra prior to any downstream analysis. 
            The quality of the identifications is not linked to this data as most search engines perform internal 
            spectrum pre-processing before matching the spectra. Thus, the spectra reported are not 
            necessarily pre-processed since the search engine may have applied the pre-processing step 
            internally. This pre-processing is not necessarily reported in the experimental metadata.
            """,
    )

# self.sub_sections["ms1"]
def draw_ms_information(
        sub_sections,
        ms1_tic=None,
        ms1_bpc=None,
        ms1_peaks=None,
        ms1_general_stats=None
    ):

    if ms1_tic:
        ms1_tic_config = {
            "id": "ms1_tic",
            "tt_label": "<b>{point.x} Ion Count:</b> {point.y}",
            "title": "Total Ion Chromatograms",
            "ylab": "Ion Count",
            "xlab": "Retention Time (seconds)",
            "ymin": 0,
            "showlegend": True,
        }
        ms1_tic_html = linegraph.plot(ms1_tic, ms1_tic_config)

        ms1_tic_html = remove_subtitle(ms1_tic_html)

        add_sub_section(
            sub_section=sub_sections,
            plot=ms1_tic_html,
            order=1,
            description="MS1 quality control information extracted from the spectrum files.",
            helptext="""
                This plot displays Total Ion Chromatograms (TICs) derived from MS1 scans across all analyzed samples. 
                The x-axis represents retention time, and the y-axis shows the total ion intensity at each time point. 
                Each colored trace corresponds to a different sample. The TIC provides a global view of the ion signal
                throughout the LC-MS/MS run, reflecting when compounds elute from the chromatography column. 
                Key aspects to assess include:
                
                * Overall intensity pattern: A consistent baseline and similar peak profiles across samples indicate good reproducibility.
                * Major peak alignment: Prominent peaks appearing at similar retention times suggest stable chromatographic performance.
                * Signal-to-noise ratio: High peaks relative to baseline noise reflect better sensitivity.
                * Chromatographic resolution: Sharp, well-separated peaks indicate effective separation.
                * Signal drift: A gradual decline in signal intensity across the run may point to source contamination or chromatography issues.
                
                Deviations such as shifted retention times, missing peaks, or inconsistent intensities may signal problems 
                in sample preparation, LC conditions, or mass spectrometer performance that require further investigation.
                """,
        )

    if ms1_bpc:

        ms1_bpc_config = {
            "id": "ms1_bpc",
            "tt_label": "<b>{point.x} Ion Count:</b> {point.y}",
            "title": "MS1 Base Peak Chromatograms",
            "ylab": "Ion Count",
            "xlab": "Retention Time (seconds)",
            "ymin": 0,
            "showlegend": True,
        }

        ms1_bpc_html = linegraph.plot(ms1_bpc, ms1_bpc_config)
        ms1_bpc_html = remove_subtitle(ms1_bpc_html)

        add_sub_section(
            sub_section=sub_sections,
            plot=ms1_bpc_html,
            order=2,
            description="MS1 base peak chromatograms extracted from the spectrum files.",
            helptext="""
                The Base Peak Chromatogram (BPC) displays the intensity of the most abundant ion at each retention
                time point across your LC-MS run. Unlike the Total Ion Chromatogram (TIC) which shows the summed
                intensity of all ions, the BPC highlights the strongest signals, providing better visualization
                of compounds with high abundance while reducing baseline noise. This makes it particularly useful
                for identifying major components in complex samples, monitoring dominant species, and providing
                clearer peak visualization when signal-to-noise ratio is a concern. Comparing BPC patterns across
                samples allows you to evaluate consistency in the detection of high-abundance compounds
                and can reveal significant variations in sample composition or instrument performance.
                """,
        )

    if ms1_peaks:

        ms1_peaks_config = {
            "id": "ms1_peaks",
            "tt_label": "<b>{point.x} Peak Count:</b> {point.y}",
            "title": "MS1 Peaks",
            "ylab": "Peak Count",
            "xlab": "Retention Time (seconds)",
            "ymin": 0,
            "showlegend": True,
        }

        ms1_peaks_html = linegraph.plot(ms1_peaks, ms1_peaks_config)
        ms1_peaks_html = remove_subtitle(ms1_peaks_html)

        add_sub_section(
            sub_section=sub_sections,
            plot=ms1_peaks_html,
            order=3,
            description="MS1 Peaks from the spectrum files",
            helptext="""
                This plot shows the number of peaks detected in MS1 scans over the course of each sample run.
                The x-axis represents retention time (in minutes), while the y-axis displays the number of
                distinct ion signals (peaks) identified in each MS1 scan. The MS1 peak count reflects spectral
                complexity and provides insight into instrument performance during the LC-MS analysis.
                Key aspects to consider include:
                
                * Overall pattern: Peak counts typically increase during the elution of complex mixtures and
                decrease during column washing or re-equilibration phases.
                * Peak density: Higher counts suggest more complex spectra, potentially indicating a greater
                number of compounds present at that time point."
                * Peak Consistency across samples: Similar profiles among replicates or related samples
                indicate good analytical reproducibility.
                * Sudden drops: Abrupt decreases in peak count may point to transient ionization issues,
                spray instability, or chromatographic disruptions.
                * Baseline values: The minimum peak count observed reflects the level of background noise
                or instrument sensitivity in the absence of eluting compounds.
                
                Monitoring MS1 peak counts complements total ion chromatogram (TIC) and base peak chromatogram
                (BPC) data, offering an additional layer of quality control related to signal complexity,
                instrument stability, and sample composition.
                """,
        )

    if ms1_general_stats:
        tconfig = {
            "id": "ms_general_stats",
            "title": "General stats for MS1 information",
            "only_defined_headers": True,
            "col1_header": "File",
        }
        headers = {
            "File": {
                "title": "File",
            },
            "AcquisitionDateTime": {
                "title": "Acquisition Date Time",
            },
            "log10(TotalCurrent)": {
                "title": "log10(Total Current)",
                "format": "{:,.4f}",
            },
            "log10(ScanCurrent)": {
                "title": "log10(Scan Current)",
                "format": "{:,.4f}",
            },
        }
        table_html = table.plot(ms1_general_stats, headers=headers, pconfig=tconfig)

        add_sub_section(
            sub_section=sub_sections,
            plot=table_html,
            order=4,
            description="General stats for MS1 information extracted from the spectrum files.",
            helptext="""
                This table presents general statistics for MS1 information extracted from mass spectrometry data files."
                It displays MS runs with their acquisition dates and times.
                For each file, the table shows two key metrics: TotalCurrent (the sum of all MS1 ion intensities throughout the run) and ScanCurrent (the sum of MS2 ion intensities).
                These values provide a quick overview of the total ion signals detected during both survey scans (MS1) and fragmentation scans (MS2),
                allowing for comparison of overall signal intensity across samples. Consistent TotalCurrent and ScanCurrent values across similar samples typically
                indicate good reproducibility in the mass spectrometry analysis, while significant variations may suggest issues with sample preparation,
                instrument performance, or ionization efficiency. The blue shading helps visualize the relative intensity differences between samples.
                """,
        )

# self.sub_sections["identification"]
def draw_quantms_identification(
        sub_sections,
        cal_num_table_data=None,
        mzml_table=None,
        quantms_missed_cleavages=None,
        quantms_modified=None,
        identified_msms_spectra=None
    ):

    # 1.ProteinGroups Count
    draw_config = {
        "id": "protein_group_count",
        "cpswitch": False,
        "title": "ProteinGroups Count",
        "tt_decimals": 0,
        "ylab": "Count",
    }

    if cal_num_table_data:
        protein_count = {
            sample: {"Count": info["protein_num"]}
            for sample, info in cal_num_table_data.items()
        }
        peptide_count = {
            sample: {"Count": info["peptide_num"]}
            for sample, info in cal_num_table_data.items()
        }
    else:
        return

    bar_html = bargraph.plot(
        protein_count,
        pconfig=draw_config,
    )
    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_sections,
        plot=bar_html,
        order=3,
        description="Number of protein groups per raw file.",
        helptext="""
            Based on statistics calculated from mzTab, mzIdentML (mzid), or DIA-NN report files.
            """,
    )

    # 2.Peptide ID Count
    draw_config = {
        "id": "peptide_id_count",
        "cpswitch": False,
        "title": "Peptide ID Count",
        "tt_decimals": 0,
        "ylab": "Count",
    }
    bar_html = bargraph.plot(
        peptide_count,
        pconfig=draw_config,
    )
    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_sections,
        plot=bar_html,
        order=4,
        description="""
            Number of unique (i.e. not counted twice) peptide sequences including modifications per Raw file.
            """,
        helptext="""
            Based on statistics calculated from mzTab, mzIdentML (mzid), or DIA-NN report files.
            """,
    )

    # 3.Missed Cleavages Per Raw File
    if quantms_missed_cleavages:
        mc_group = {}
        for sample, counts in quantms_missed_cleavages.items():
            mc_group[sample] = {"0": 0, "1": 0, ">=2": 0}
            for mc, count in counts.items():
                if mc == 0:
                    mc_group[sample]["0"] += count
                elif mc == 1:
                    mc_group[sample]["1"] += count
                else:
                    mc_group[sample][">=2"] += count

        mc_group_ratio = {}
        for sample, counts in mc_group.items():
            total = sum(counts.values())
            if total == 0:
                mc_group_ratio[sample] = {"0": 0, "1": 0, ">=2": 0}
                continue
            mc_group_ratio[sample] = {
                group: count / total * 100 for group, count in counts.items()
            }

        mc_data = {"plot_data": mc_group_ratio, "cats": ["0", "1", ">=2"]}
        draw_msms_missed_cleavages(sub_sections, mc_data, False)

    # 4.Modifications Per Raw File
    if quantms_modified:
        draw_modifications(sub_sections, quantms_modified)

    # 5.MS/MS Identified per Raw File
    if identified_msms_spectra and mzml_table:

        msms_identified_rate = dict()

        for m in identified_msms_spectra.keys():
            identified_ms2 = identified_msms_spectra[m].get("Identified", 0)
            all_ms2 = mzml_table.get(m, {}).get("MS2_Num", 0)

            if all_ms2 > 0:
                msms_identified_rate[m] = {"Identified Rate": (identified_ms2 / all_ms2) * 100}

        draw_ms_ms_identified(sub_sections, msms_identified_rate)

# self.sub_sections["ms2"]
def draw_precursor_charge_distribution(
        sub_sections,
        charge_plot=None,
        ms_info=None
    ):
    pconfig = {
        "id": "distribution_of_precursor_charges",
        "title": "Distribution of Precursor Charges",
        "cpswitch": True,
        "tt_decimals": 0,
        "ylab": "Count",
    }
    # if config.kwargs.get("mzid_plugin", False) and self.mgf_paths:
    #     cats = self.mgf_charge_plot.dict["cats"]
    # else:
    #     cats = self.mzml_charge_plot.dict["cats"]
    cats = charge_plot.dict["cats"]

    bar_html = bargraph.plot(ms_info["charge_distribution"], cats, pconfig)
    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_sections,
        plot=bar_html,
        order=5,
        description="""
            This is a bar chart representing the distribution of the precursor ion charges for a given whole experiment.
            """,
        helptext="""
            This information can be used to identify potential ionization problems 
            including many 1+ charges from an ESI ionization source or an unexpected 
            distribution of charges. MALDI experiments are expected to contain almost exclusively 1+ 
            charged ions. An unexpected charge distribution may furthermore be caused by specific search 
            engine parameter settings such as limiting the search to specific ion charges.
            """,
    )
