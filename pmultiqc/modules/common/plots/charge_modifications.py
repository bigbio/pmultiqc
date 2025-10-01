"""
Charge state and modification-related plotting functions for PMultiQC

This module contains functions for plotting charge state distributions and modifications.
"""

from multiqc.plots import bargraph
from pmultiqc.modules.core.section_groups import add_sub_section
from .contaminants import remove_subtitle


def draw_charge_state(sub_section, charge_data, report_type):
    """
    Draw charge state distribution plot per file.

    Args:
        sub_section: MultiQC sub-section for adding plots
        charge_data: Dictionary containing charge distribution data
        report_type: String indicating the type of report (MaxQuant, mzIdentML, etc.)
    """
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


def draw_modifications(sub_section, modified_data):
    """
    Draw modifications plot per raw file.

    Args:
        sub_section: MultiQC sub-section for adding plots
        modified_data: Dictionary containing modification data
    """
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


def draw_precursor_charge_distribution(sub_sections, charge_plot=None, ms_info=None):
    """
    Draw precursor charge distribution plot.

    Args:
        sub_sections: MultiQC sub-sections for adding plots
        charge_plot: Plot data object containing charge statistics
        ms_info: Dictionary containing MS information
    """
    pconfig = {
        "id": "distribution_of_precursor_charges",
        "title": "Distribution of Precursor Charges",
        "cpswitch": True,
        "tt_decimals": 0,
        "ylab": "Count",
    }
    
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
