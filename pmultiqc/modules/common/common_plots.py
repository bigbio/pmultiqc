from multiqc.plots import (
    bargraph,
    linegraph,
    heatmap
)
from ..core.section_groups import add_sub_section


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

    bar_html = bargraph.plot(
        data=msms_identified_percent,
        pconfig=draw_config
    )

    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=7,
        description="MS/MS identification rate per Raw file.",
        helptext="MS/MS identification rate per raw file (quantms data from mzTab and mzML files; MaxQuant data from summary.txt)"
    )

# Potential Contaminants per Group
def draw_potential_contaminants(
        sub_section,
        contaminant_percent,
        is_maxquant
    ):

    draw_config = {
        "id": "potential_contaminants_per_group",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Potential Contaminants Per File",
        "tt_decimals": 2,
        "ylab": "Percent [%]",
    }

    bar_html = bargraph.plot(
        data=contaminant_percent,
        pconfig=draw_config
    )

    bar_html = remove_subtitle(bar_html)

    if is_maxquant:
        description_text = "Potential contaminants per group from proteinGroups.txt."
        help_text = """
            External protein contamination should be controlled for, therefore MaxQuant ships with a 
            comprehensive, yet customizable protein contamination database, which is searched by MaxQuant 
            by default. A contamination plot derived from the proteinGroups (PG) table, showing the 
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
        helptext=help_text
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
            """
    )

# Charge-state of Per File
def draw_charge_state(sub_section, charge_data, is_maxquant):

    draw_config = {
        "id": "charge_state_of_per_raw_file",
        "cpswitch": True,
        "title": "Charge-state of Per File",
        "tt_decimals": 0,
        "ylab": "Count",
    }

    bar_html = bargraph.plot(
        data=charge_data["plot_data"],
        cats=charge_data["cats"],
        pconfig=draw_config
    )

    bar_html = remove_subtitle(bar_html)

    if is_maxquant:
        description_text = "The distribution of the charge-state of the precursor ion, excluding potential contaminants."
        help_text = "The distribution of the charge-state of the precursor ion, excluding potential contaminants."
    else:
        description_text = "The distribution of precursor ion charge states (based on mzTab data)."
        help_text = "The distribution of precursor ion charge states (based on mzTab data)."

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=6,
        description=description_text,
        helptext=help_text
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
        data=modified_data["plot_data"],
        cats=modified_data["cats"],
        pconfig=draw_config
    )

    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=6,
        description="""
            Compute an occurence table of modifications (e.g. Oxidation (M)) for all peptides, including the unmodified.
            """,
        helptext="Post-translational modifications contained within the identified peptide sequence."
    )

def draw_oversampling(
        sub_section,
        oversampling,
        oversampling_plot,
        is_maxquant
    ):

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
            data=oversampling["plot_data"],
            cats=oversampling["cats"],
            pconfig=draw_config
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

        bar_html = bargraph.plot(
            data=oversampling,
            cats=oversampling_plot,
            pconfig=draw_config
        )

    bar_html = remove_subtitle(bar_html)
    
    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=7,
        description="""
            An oversampled 3D-peak is defined as a peak whose peptide ion 
            (same sequence and same charge state) was identified by at least two distinct MS2 spectra 
            in the same Raw file.
            """,
        helptext="""
            For high complexity samples, oversampling of individual 3D-peaks automatically leads to 
            undersampling or even omission of other 3D-peaks, reducing the number of identified peptides. 
            Oversampling occurs in low-complexity samples or long LC gradients, as well as undersized dynamic 
            exclusion windows for data independent acquisitions.
            """
    )

# Missed Cleavages Per Raw File
def draw_msms_missed_cleavages(
        sub_section,
        missed_cleavages_data,
        is_maxquant
    ):

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

    if is_maxquant:
        description_text = "[Excludes Contaminants] Missed Cleavages per raw file."
    else:
        description_text = "Missed Cleavages per raw file."

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=5,
        description=description_text,
        helptext="""
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

            If MC>=1 is high (>20%) you should increase the missed cleavages settings in MaxQuant and compare the number of peptides.
            Usually high MC correlates with bad identification rates, since many spectra cannot be matched to the forward database.

            In the rare case that 'no enzyme' was specified in MaxQuant, neither scores nor plots are shown.
            """
    )

# IDs over RT
def draw_ids_rt_count(sub_section, rt_count_data, data_type):

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

    linegraph_html = linegraph.plot(
        data=rt_count_data,
        pconfig=draw_config
    )

    linegraph_html = remove_subtitle(linegraph_html)

    if data_type == "maxquant":
        description_text = "Distribution of retention time, derived from the evidence table."
        help_text = """
            The uncalibrated retention time in minutes in the elution profile of the precursor ion, 
            and does not include potential contaminants.
            """
    elif data_type == "dia":
        description_text = "Distribution of retention time, derived from the report.tsv."
        help_text = """
            [DIA-NN: report.tsv] Distribution of retention time (RT) for each run.
            """
    else:
        description_text = "Distribution of retention time, derived from the mzTab."
        help_text = """
            The uncalibrated retention time in minutes in the elution profile of the precursor ion.
            """
        
    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=1,
        description=description_text,
        helptext=help_text
    )


# MaxQuant: Delta Mass [Da], Delta Mass [ppm]
# quantms: Delta Mass [ppm]
def draw_delta_mass_da_ppm(sub_section, delta_mass, delta_mass_type):

    # MaxQuant: Delta Mass [Da]
    if delta_mass_type == "Mass Error [Da]":
        plot_id = "delta_mass_da"
        plot_title = "Delta Mass [Da]"
        plot_xlab = "Experimental m/z - Theoretical m/z"
        sub_section_order = 3
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
        sub_section_order = 4
        description_text = """
            This plot is based on the "Mass Error [ppm]" column from the evidence.txt generated by MaxQuant.
            """
        help_text = """
            Mass error of the recalibrated mass-over-charge value of the precursor ion in comparison to the 
            predicted monoisotopic mass of the identified peptide sequence in parts per million.
            """
    
    # quantms: Delta Mass [ppm]
    elif delta_mass_type == "quantms_ppm":
        plot_id = "delta_mass_ppm"
        plot_title = "Delta Mass [ppm]"
        plot_xlab = "Delta Mass [ppm]"
        sub_section_order = 4
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
        delta_mass_percent_range = {k: v for k, v in delta_mass["frequency"].items() if abs(k) <= range_abs}

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
                {"relative_frequency": delta_mass["frequency"]}
            ],
            pconfig
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
            [
                {"count": delta_mass["count"]},
                {"relative_frequency": delta_mass["frequency"]}
            ],
            pconfig
        )

    add_sub_section(
        sub_section=sub_section,
        plot=line_html,
        order=sub_section_order,
        description=description_text,
        helptext=help_text
    )


def draw_heatmap(
        sub_sections,
        hm_colors,
        heatmap_data,
        heatmap_xnames,
        heatmap_ynames,
        is_maxquant
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

        hm_html = heatmap.plot(
            data=heatmap_data,
            pconfig=pconfig
        )
        
        description_text = "This heatmap provides an overview of the performance of the MaxQuant results."

    else:

        hm_html = heatmap.plot(
            heatmap_data,
            heatmap_xnames,
            heatmap_ynames,
            pconfig
        )

        description_text = "This heatmap provides an overview of the performance of the quantms."

    hm_html = remove_subtitle(hm_html)

    add_sub_section(
        sub_section=sub_sections,
        plot=hm_html,
        order=2,
        description=description_text,
        helptext="""
            This plot shows the pipeline performance overview. Some metrics are calculated.
            
            * Heatmap score[Contaminants]: as fraction of summed intensity with 0 = sample full of contaminants; 
                1 = no contaminants
            * Heatmap score[Pep Intensity (>23.0)]: Linear scale of the median intensity reaching the threshold, 
                i.e. reaching 2^21 of 2^23 gives score 0.25.
            * Heatmap score[Charge]: Deviation of the charge 2 proportion from a representative 
                Raw file (median). For typtic digests, peptides of charge 2 (one N-terminal and one at 
                tryptic C-terminal R or K residue) should be dominant. Ionization issues (voltage?), 
                in-source fragmentation, missed cleavages and buffer irregularities can cause a shift 
                (see Bittremieux 2017, DOI: 10.1002/mas.21544).
            * Heatmap score [Missed Cleavages]: the fraction (0% - 100%) of fully cleaved peptides per Raw file
            * Heatmap score [Missed Cleavages Var]: each Raw file is scored for its deviation from the ‘average’ digestion 
                state of the current study.
            * Heatmap score [ID rate over RT]: Judge column occupancy over retention time. 
                Ideally, the LC gradient is chosen such that the number of identifications 
                (here, after FDR filtering) is uniform over time, to ensure consistent instrument duty cycles. 
                Sharp peaks and uneven distribution of identifications over time indicate potential for LC gradient 
                optimization.Scored using ‘Uniform’ scoring function. i.e. constant receives good score, extreme shapes are bad.
            * Heatmap score [MS2 Oversampling]: The percentage of non-oversampled 3D-peaks. An oversampled 
                3D-peak is defined as a peak whose peptide ion (same sequence and same charge state) was 
                identified by at least two distinct MS2 spectra in the same Raw file. For high complexity samples, 
                oversampling of individual 3D-peaks automatically leads to undersampling or even omission of other 3D-peaks, 
                reducing the number of identified peptides.
            * Heatmap score [Pep Missing Values]: Linear scale of the fraction of missing peptides.
            
            """
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
