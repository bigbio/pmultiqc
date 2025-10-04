"""
General plotting functions and utilities for PMultiQC

This module contains general plotting functions like heatmaps and utility functions.
"""

from matplotlib.colors import LinearSegmentedColormap, to_hex
from multiqc.plots import heatmap
from pmultiqc.modules.core.section_groups import add_sub_section
from .contaminants import remove_subtitle

# HeatMap color list
color_map = LinearSegmentedColormap.from_list("red_green", ["#ff0000", "#00ff00"])
HEATMAP_COLOR_LIST = [[s, to_hex(color_map(s))] for s in [round(i * 0.1, 1) for i in range(11)]]


def draw_heatmap(sub_sections, hm_colors, heatmap_data, heatmap_xnames, heatmap_ynames, is_maxquant):
    """
    Draw heatmap plot.

    Args:
        sub_sections: MultiQC sub-sections for adding plots
        hm_colors: Color configuration for heatmap
        heatmap_data: Dictionary containing heatmap data
        heatmap_xnames: List of x-axis names
        heatmap_ynames: List of y-axis names
        is_maxquant: Boolean indicating if data is from MaxQuant
    """
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
