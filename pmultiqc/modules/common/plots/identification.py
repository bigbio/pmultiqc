"""
Identification-related plotting functions for PMultiQC

This module contains functions for plotting identification statistics and results.
"""

from multiqc.plots import bargraph
from pmultiqc.modules.core.section_groups import add_sub_section
from .contaminants import remove_subtitle


def draw_ms_ms_identified(sub_section, msms_identified_percent):
    """
    Draw MS/MS identified plot per raw file.

    Args:
        sub_section: MultiQC sub-section for adding plots
        msms_identified_percent: Dictionary with MS/MS identification percentage data
    """
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


def draw_num_pep_per_protein(sub_sections, pep_plot):
    """
    Draw number of peptides per protein plot.

    Args:
        sub_sections: MultiQC sub-sections for adding plots
        pep_plot: Plot data object containing peptide statistics
    """
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

    from multiqc import config
    
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


def draw_quantms_identification(sub_sections, cal_num_table_data=None, mzml_table=None, 
                               quantms_missed_cleavages=None, quantms_modified=None, 
                               identified_msms_spectra=None):
    """
    Draw quantms identification plots including protein counts, peptide counts, missed cleavages, etc.

    Args:
        sub_sections: MultiQC sub-sections for adding plots
        cal_num_table_data: Dictionary with calculated number table data
        mzml_table: Dictionary with mzML table data
        quantms_missed_cleavages: Dictionary with missed cleavages data
        quantms_modified: Dictionary with modification data
        identified_msms_spectra: Dictionary with identified MS/MS spectra data
    """
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

    bar_html = bargraph.plot(protein_count, pconfig=draw_config)
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
    bar_html = bargraph.plot(peptide_count, pconfig=draw_config)
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
        from .mass_analysis import draw_msms_missed_cleavages
        
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
        from .charge_modifications import draw_modifications
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
