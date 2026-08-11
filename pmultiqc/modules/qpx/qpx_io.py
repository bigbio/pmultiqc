from __future__ import absolute_import

from pmultiqc.modules.common.logging import get_logger
from multiqc import config

import pandas as pd
from pmultiqc.modules.common.file_utils import file_prefix


# Initialise the module logger via centralized logger
log = get_logger("pmultiqc.modules.qpx.qpx_io")


QPX_COLUMNS = {
    "psm": [
        "sequence", "peptidoform",
        "modifications", "charge",
        "posterior_error_probability", "is_decoy",
        "calculated_mz", "observed_mz",
        "mass_error_ppm", "run_file_name",
        "scan", "rt", "protein_accessions"
    ],
    "pg": [
        "pg_accessions", "anchor_protein",
        "grouped_runs", "global_qvalue",
        "intensity", "is_decoy"
    ],
    "feature": [
        "feature_id", "sequence", "peptidoform",
        "charge", "is_decoy", "run_file_name",
        "intensities", "anchor_protein"
    ]
}


def parse_qpx_parquet(file_path, qpx_type):

    req_columns = QPX_COLUMNS[qpx_type]
    if "is_decoy" not in req_columns:
        req_columns = list(req_columns) + ["is_decoy"]

    remove_decoy = config.kwargs.get("remove_decoy", False)
    parquet_filters = [("is_decoy", "==", False)] if remove_decoy else None

    df = pd.read_parquet(
        path=file_path,
        columns=req_columns,
        engine="pyarrow",
        filters=parquet_filters
    )

    df = df.drop(columns=["is_decoy"], errors="ignore")

    # psm/feature: run_file_name --> run
    if "run_file_name" in df.columns:
        df["run"] = df["run_file_name"].apply(file_prefix)

        df = df.drop(
            columns=["run_file_name"],
            errors="ignore"
        )

    # pg: grouped_runs --> run
    if "grouped_runs" in df.columns:
        df = df.explode("grouped_runs")
        df["run"] = df["grouped_runs"].apply(file_prefix)

        df = df.drop(
            columns=["grouped_runs"],
            errors="ignore"
        )

    log.info(f"[Loaded data: {qpx_type}.parquet] {file_path}, Data shape: {df.shape}.")

    return df
