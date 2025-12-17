from typing import Dict, List
from collections import OrderedDict

from multiqc.plots import bargraph, linegraph, table
from multiqc.types import SampleGroup, SampleName
from multiqc.plots.table_object import InputRow

from pmultiqc.modules.core.section_groups import add_sub_section
from pmultiqc.modules.common.common_utils import condition_split
from pmultiqc.modules.common.plots.general import remove_subtitle
from pmultiqc.modules.common.tooltip_config import apply_tooltip_config


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

    help_text = (
        "Charge distribution per Raw file. For typtic digests, peptides of charge 2 "
        "(one N-terminal and one at tryptic C-terminal R or K residue) should be dominant. "
        "Ionization issues (voltage?), in-source fragmentation, missed cleavages and buffer irregularities can "
        "cause a shift (see Bittremieux 2017, DOI: 10.1002/mas.21544). The charge distribution should be similar "
        "across Raw files. Consistent charge distribution is paramount for comparable 3D-peak intensities across samples."
    )

    if report_type == "MaxQuant":
        description_text = (
            "The distribution of the charge-state of the precursor ion, excluding potential contaminants."
        )
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


def draw_delta_mass_da_ppm(sub_section, delta_mass, delta_mass_type):
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

    elif delta_mass_type == "quantms_ppm":
        plot_id = "delta_mass_ppm"
        plot_title = "Delta Mass [ppm]"
        plot_xlab = "Delta Mass [ppm]"
        sub_section_order = 2
        description_text = """
            Delta Mass [ppm] calculated from mzTab.
            """
        help_text = """
            Delta Mass [ppm] calculated from mzTab: ((experimental m/z - theoretical m/z) / (theoretical m/z)) x 10^6.
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


def draw_quantms_identification(
        sub_sections,
        cal_num_table_data=None,
        mzml_table=None,
        quantms_missed_cleavages=None,
        quantms_modified=None,
        identified_msms_spectra=None
):
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

    if quantms_modified:
        draw_modifications(sub_sections, quantms_modified)

    if identified_msms_spectra and mzml_table:

        msms_identified_rate = dict()
        for m in identified_msms_spectra.keys():
            identified_ms2 = identified_msms_spectra[m].get("Identified", 0)
            all_ms2 = mzml_table.get(m, {}).get("MS2_Num", 0)
            if all_ms2 > 0:
                msms_identified_rate[m] = {"Identified Rate": (identified_ms2 / all_ms2) * 100}

        draw_ms_ms_identified(sub_sections, msms_identified_rate)


def draw_summary_protein_ident_table(
        sub_sections,
        enable_dia: bool = False,
        total_peptide_count: int = 0,
        total_protein_quantified: int = 0,
        total_ms2_spectra: int = 0,
        total_ms2_spectra_identified: int = 0,
        total_protein_identified: int = 0,
        enable_mzid: bool = False
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
                "#Identified MS2 Spectra": total_ms2_spectra_identified
            }
        }
        coverage = (total_ms2_spectra_identified / total_ms2_spectra * 100) if total_ms2_spectra > 0 else 0.0
        summary_table[total_ms2_spectra]["%Identified MS2 Spectra"] = coverage
        summary_table[total_ms2_spectra]["#Peptides Identified"] = total_peptide_count
        summary_table[total_ms2_spectra]["#Proteins Identified"] = total_protein_identified

        if not enable_mzid:
            summary_table[total_ms2_spectra]["#Proteins Quantified"] = total_protein_quantified

        headers["#Identified MS2 Spectra"] = {
            "description": "Total number of MS/MS spectra identified",
        }
        headers["%Identified MS2 Spectra"] = {
            "description": "Percentage of Identified MS/MS Spectra",
            "format": "{:,.2f}",
            "suffix": "%",
        }
        col_header = "#MS2 Spectra"

    pconfig = {
        "id": "identification_summary_table",
        "title": "Summary Table",
        "save_file": False,
        "raw_data_fn": "multiqc_summary_table_table",
        "sort_rows": False,
        "only_defined_headers": False,
        "col1_header": col_header,
        "scale": "Set1",
    }
    table_html = table.plot(summary_table, headers, pconfig)

    if enable_mzid or enable_dia:
        description_str = "This table shows the summary statistics of the submitted data."
        helptext_str = """
            This table shows the summary statistics of the submitted data.
          """
    else:
        description_str = "This table shows the quantms pipeline summary statistics."
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


def draw_quantms_identi_num(
        sub_sections,
        enable_exp=False,
        enable_sdrf=False,
        is_multi_conditions=False,
        sample_df=None,
        file_df=None,
        cal_num_table_data=None
):
    rows_by_group: Dict[SampleGroup, List[InputRow]] = {}

    if enable_exp or enable_sdrf:
        if is_multi_conditions:
            for sample in sorted(
                    sample_df["Sample"].tolist(),
                    key=lambda x: (str(x).isdigit(), int(x) if str(x).isdigit() else str(x).lower()),
            ):
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
                headers["MSstats_Condition_" + str(k)] = {
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
            for sample in sorted(
                    sample_df["Sample"].tolist(),
                    key=lambda x: (str(x).isdigit(), int(x) if str(x).isdigit() else str(x).lower()),
            ):

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

            <p>The plot will show percentages, i.e. is normalized by the total number of peptide sequences 
            (where different charge state counts as a separate peptide) per Raw file.
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
            If a modification, e.g. Oxidation(M), occurs multiple times in a single peptide it's listed as a separate modification 
            (e.g. '2 Oxidation (M)' for double oxidation of a single peptide).</p>
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


def draw_num_pep_per_protein(
        sub_sections,
        pep_plot,
        enable_mzid: bool = False
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

    if enable_mzid:
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