import itertools
import logging

import numpy as np
import pandas as pd
from multiqc.plots import table

from pmultiqc.modules.common.file_utils import drop_empty_row
from pmultiqc.modules.maxquant.maxquant_utils import evidence_rt_count
from ..core.section_groups import add_sub_section

from pmultiqc.modules.common.logging import get_logger
log = get_logger("pmultiqc.modules.common.mzidentml_utils")


def get_mzidentml_mzml_df(mzid_psm, mzml_ms_df):

    # check filename
    mzid_filenames = set(mzid_psm["filename"].unique())
    mzml_filenames = set(mzml_ms_df["filename"].unique())
    only_in_mzid = mzid_filenames - mzml_filenames
    only_in_mzml = mzml_filenames - mzid_filenames

    if len(only_in_mzid) + len(only_in_mzml) != 0:
        log.warning(
            f"Please check the filename. Unique to mzIdentML: {only_in_mzid}, unique to mzML: {only_in_mzml}"
        )
        return pd.DataFrame()
    else:
        if "retention_time" in mzid_psm.columns:
            mzml_mzid_df = pd.merge(
                mzid_psm,
                mzml_ms_df[["spectrumID", "intensity", "filename"]],
                on=["spectrumID", "filename"],
                how="inner",
            )
        else:
            mzml_mzid_df = pd.merge(
                mzid_psm, mzml_ms_df, on=["spectrumID", "filename"], how="inner"
            )

        return mzml_mzid_df


# ProteinGroups Count / Peptide ID Count
def get_mzid_num_data(df):

    num_data = dict()
    identified_spectra = dict()
    for file_name, group in df.groupby("filename"):

        num_data[file_name] = {
            "protein_num": group["accession_group"].nunique(),
            "peptide_num": len(group[["Modifications", "PeptideSequence"]].drop_duplicates()),
        }

        identified_spectra[file_name] = {"Identified": len(set(group["spectrumID"]))}

    return num_data, identified_spectra


# Charge-state of Per File
def get_mzidentml_charge(df):
    charge_state_df = df.groupby(["filename", "chargeState"]).size().unstack(fill_value=0)
    charge_state_df.rename(columns=lambda x: str(x), inplace=True)

    charge_state = {
        "plot_data": charge_state_df.to_dict(orient="index"),
        "cats": charge_state_df.columns.tolist(),
    }

    return charge_state


# IDs over RT
def get_mzid_rt_id(df):

    rt_file_df = df[["filename", "retention_time"]].copy()
    rt_file_df.rename(
        columns={"filename": "raw file", "retention_time": "retention time"}, inplace=True
    )
    result = evidence_rt_count(rt_file_df)

    return result


# mzIdentML: Peptides Quantification Table
def create_peptides_table(mzid_spectrum_df):

    # Validation: remove rows with 0 or NA intensity values
    report_data = mzid_spectrum_df[mzid_spectrum_df["intensity"] > 0].copy()
    report_data = drop_empty_row(report_data, ["accession_group", "PeptideSequence"])

    table_dict = dict()
    for sequence_protein, group in report_data.groupby(["PeptideSequence", "accession_group"]):

        table_dict[sequence_protein] = {
            "ProteinName": sequence_protein[1],
            "PeptideSequence": sequence_protein[0],
            "BestSearchScore": group["search_engine_score"].max(),
            "Average Intensity": np.log10(group["intensity"].mean()),
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
            "format": "{:,.4f}",
        },
    }

    result_dict = {i: v for i, (_, v) in enumerate(table_dict.items(), start=1)}

    return result_dict, headers


# mzIdentML: Protein Quantification Table
def create_protein_table(mzid_spectrum_df):

    # Validation: remove rows with 0 or NA intensity values
    report_data = mzid_spectrum_df[mzid_spectrum_df["intensity"] > 0].copy()
    report_data = drop_empty_row(report_data, ["accession_group", "PeptideSequence"])

    table_dict = dict()
    for protein_name, group in report_data.groupby("accession_group"):

        table_dict[protein_name] = {
            "ProteinName": protein_name,
            "Peptides_Number": group["PeptideSequence"].nunique(),
            "Average Intensity": np.log10(group["intensity"].mean()),
        }

    headers = {
        "ProteinName": {
            "title": "Protein Name",
            "description": "Name/Identifier(s) of the protein (group)",
        },
        "Peptides_Number": {
            "title": "Number of Peptides",
            "description": "Number of peptides per proteins",
            "format": "{:,.0f}",
        },
        "Average Intensity": {
            "title": "Average Intensity",
            "format": "{:,.4f}",
        },
    }

    result_dict = {i: v for i, (_, v) in enumerate(table_dict.items(), start=1)}

    return result_dict, headers


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
        helptext=helptext_text,
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

    # only use the first 50 lines for the table
    display_rows = 50
    table_html = table.plot(
        dict(itertools.islice(table_data.items(), display_rows)),
        headers=headers,
        pconfig=draw_config,
    )

    if report_type == "DIA-NN":
        description_text = """
            This plot shows the quantification information of proteins in the final result (DIA-NN report).
            """
        helptext_text = """
            The quantification information of proteins is obtained from the DIA-NN output file. 
            The table shows the quantitative level and distribution of proteins in different study variables, 
            run and peptiforms. The distribution show all the intensity values in a bar plot above and below 
            the average intensity for all the fractions, runs and peptiforms.

            * Peptides_Number: Number of peptides per protein.
            * Average Intensity: Average intensity of each protein sequence across all conditions (0 or NA ignored).
            * Protein intensity in each condition (Eg. `CT=Mixture;CN=UPS1;QY=0.1fmol`).
            """
    elif report_type == "mzIdentML":
        description_text = """
            This plot shows the quantification information of proteins in the final result (mzIdentML).
            """
        helptext_text = """
            The quantification information of proteins is obtained from the mzIdentML. 
            The table shows the quantitative level and distribution of proteins in different study variables, 
            run and peptiforms. The distribution show all the intensity values in a bar plot above and below 
            the average intensity for all the fractions, runs and peptiforms.

            * Peptides_Number: Number of peptides per protein.
            * Average Intensity: Average intensity of each protein sequence (0 or NA ignored).
            """
    else:
        description_text = ""
        helptext_text = ""

    add_sub_section(
        sub_section=sub_sections,
        plot=table_html,
        order=2,
        description=description_text,
        helptext=helptext_text,
    )


# mzIdentML: Quantification Table
def draw_mzid_quant_table(sub_section, mzid_mzml_df):

    # Peptides Quantification Table
    peptides_table, peptides_headers = create_peptides_table(mzid_mzml_df)
    draw_peptides_table(sub_section, peptides_table, peptides_headers, "mzIdentML")

    # Protein Quantification Table
    protein_table, protein_headers = create_protein_table(mzid_mzml_df)
    draw_protein_table(sub_section, protein_table, protein_headers, "mzIdentML")