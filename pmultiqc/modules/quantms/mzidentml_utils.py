import logging
import pandas as pd
from pmultiqc.modules.common.file_utils import drop_empty_row
from pmultiqc.modules.maxquant.maxquant_utils import evidence_rt_count
import numpy as np

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

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
                how="inner"
            )
        else:
            mzml_mzid_df = pd.merge(
                mzid_psm,
                mzml_ms_df,
                on=["spectrumID", "filename"],
                how="inner"
            )
        
        return mzml_mzid_df

# Charge-state of Per File
def get_mzid_mzml_charge(df):
    charge_state_df = df.groupby(["filename", "chargeState"]).size().unstack(fill_value=0)
    charge_state_df.rename(columns=lambda x: str(x), inplace=True)

    charge_state = {
        "plot_data": charge_state_df.to_dict(orient="index"),
        "cats": charge_state_df.columns.tolist()
    }

    return charge_state

# IDs over RT
def get_mzid_rt_id(df):

    rt_file_df = df[["filename", "retention_time"]].copy()
    rt_file_df.rename(columns={
        "filename": "raw file",
        "retention_time": "retention time"
    }, inplace=True)
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
            "Average Intensity": np.log10(
                group["intensity"].mean()
            )
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


# mzIdentML: Peptides Quantification Table
def create_protein_table(mzid_spectrum_df):

    # Validation: remove rows with 0 or NA intensity values
    report_data = mzid_spectrum_df[mzid_spectrum_df["intensity"] > 0].copy()
    report_data = drop_empty_row(report_data, ["accession_group", "PeptideSequence"])

    table_dict = dict()
    for protein_name, group in report_data.groupby("accession_group"):

        table_dict[protein_name] = {
            "ProteinName": protein_name,
            "Peptides_Number": group["PeptideSequence"].nunique(),
            "Average Intensity": np.log10(
                group["intensity"].mean()
            ),
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
