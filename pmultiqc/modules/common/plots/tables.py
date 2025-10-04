"""
Table plotting functions for PMultiQC

This module contains functions for plotting various tables and data summaries.
"""

import itertools
import numpy as np
import re
from collections import OrderedDict
from typing import Dict, List

from multiqc.plots import table
from multiqc.plots.table_object import InputRow
from multiqc.types import SampleGroup, SampleName
from multiqc import config

from pmultiqc.modules.core.section_groups import add_sub_section
from pmultiqc.modules.common.utils import read_openms_design, condition_split


def draw_peptides_table(sub_section, table_data, headers, report_type):
    """
    Draw peptides quantification table.

    Args:
        sub_section: MultiQC sub-section for adding plots
        table_data: Dictionary containing peptide table data
        headers: Dictionary containing table headers configuration
        report_type: String indicating the type of report
    """
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


def draw_protein_table(sub_sections, table_data, headers, report_type):
    """
    Draw protein quantification table.

    Args:
        sub_sections: MultiQC sub-sections for adding plots
        table_data: Dictionary containing protein table data
        headers: Dictionary containing table headers configuration
        report_type: String indicating the type of report
    """
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


def draw_exp_design(sub_sections, exp_design):
    """
    Draw experimental design table.

    Args:
        sub_sections: MultiQC sub-sections for adding plots
        exp_design: Path to experimental design file

    Returns:
        Tuple of (sample_df, file_df, exp_design_runs, is_bruker, is_multi_conditions)
    """
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


def draw_summary_protein_ident_table(sub_sections, enable_dia: bool = False, total_peptide_count: int = 0,
                                   total_protein_quantified: int = 0, total_ms2_spectra: int = 0,
                                   total_ms2_spectral_identified: int = 0, total_protein_identified: int = 0):
    """
    Draw summary protein identification table.

    Args:
        sub_sections: MultiQC sub-sections for adding plots
        enable_dia: Boolean indicating if DIA is enabled
        total_peptide_count: Total number of peptides
        total_protein_quantified: Total number of quantified proteins
        total_ms2_spectra: Total number of MS2 spectra
        total_ms2_spectral_identified: Total number of identified MS2 spectra
        total_protein_identified: Total number of identified proteins
    """
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
        summary_table[total_ms2_spectra]["#Peptides Identified"] = total_peptide_count
        summary_table[total_ms2_spectra]["#Proteins Identified"] = total_protein_identified

        if not config.kwargs.get("mzid_plugin", False):
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

    if config.kwargs.get("mzid_plugin", False) or config.kwargs.get("diann_plugin", False):
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


def draw_quantms_identi_num(sub_sections, enable_exp, enable_sdrf, is_multi_conditions,
                          sample_df, file_df, cal_num_table_data):
    """
    Draw quantms identification numbers table.

    Args:
        sub_sections: MultiQC sub-sections for adding plots
        enable_exp: Boolean indicating if experimental design is enabled
        enable_sdrf: Boolean indicating if SDRF is enabled
        is_multi_conditions: Boolean indicating if multiple conditions exist
        sample_df: DataFrame containing sample information
        file_df: DataFrame containing file information
        cal_num_table_data: Dictionary containing calculated number table data
    """
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
