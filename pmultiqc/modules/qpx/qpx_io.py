from __future__ import absolute_import

from pmultiqc.modules.common.logging import get_logger
from multiqc import config

import pandas as pd
import pyarrow.parquet as pq

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
        "intensity", "is_decoy",
        # Optional (absent in older writers); skipped automatically when missing.
        "contaminant"
    ],
    "feature": [
        "feature_id", "sequence", "peptidoform",
        "charge", "is_decoy", "run_file_name",
        "intensities", "anchor_protein",
        # Needed so peptide counts key on the same protein-group identity as pg.parquet;
        # keying on anchor_protein alone silently loses every multi-accession group.
        # Optional (absent in older writers); skipped automatically when missing.
        "pg_accessions", "modifications", "rt",
        "calculated_mz", "observed_mz", "unique"
    ]
}


def has_data(df, *columns):
    """True only if every named column exists AND holds at least one non-null value.

    quantms.io declares many fields as optional, and different producers populate
    different subsets -- a DIA-NN-sourced project fills fields a DDA/OpenMS one leaves
    null, and vice versa. Several columns are present in the schema but entirely empty
    in real files (mass_error_ppm, missed_cleavages, unique, ...). Presence is therefore
    not usability: every section must test the data, not the schema.
    """
    if df is None or getattr(df, "empty", True):
        return False

    for column in columns:
        if column not in df.columns:
            return False
        try:
            if not df[column].notna().any():
                return False
        except (TypeError, ValueError):
            return False

    return True


def select_columns(df, columns, context):
    """Return ``df[columns]`` only if every column is present, else ``None``.

    Because the reader skips columns a file does not have, a loaded table may
    legitimately lack one. Selecting a fixed list would raise KeyError, which the
    module's _safe_draw swallows -- silently dropping a whole section. Callers use
    this to skip that plot deliberately, with a log line saying why.
    """
    if df is None or getattr(df, "empty", True):
        return None

    missing = [column for column in columns if column not in df.columns]
    if missing:
        log.warning(
            f"[{context}] Skipped: required column(s) not present in the parquet data: "
            f"{', '.join(missing)}."
        )
        return None

    return df[list(columns)].copy()


def parse_qpx_parquet(file_path, qpx_type):

    req_columns = QPX_COLUMNS[qpx_type]
    if "is_decoy" not in req_columns:
        req_columns = list(req_columns) + ["is_decoy"]

    # Only request columns the file actually has: quantms.io writers vary by version,
    # and pyarrow raises rather than ignoring an absent column.
    available = set(pq.ParquetFile(file_path).schema_arrow.names)
    missing = [column for column in req_columns if column not in available]
    if missing:
        log.debug(f"{qpx_type}.parquet has no {', '.join(missing)}; reading without them.")
    req_columns = [column for column in req_columns if column in available]

    remove_decoy = config.kwargs.get("remove_decoy", False)
    parquet_filters = (
        [("is_decoy", "==", False)] if remove_decoy and "is_decoy" in available else None
    )

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
