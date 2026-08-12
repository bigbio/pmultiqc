"""Reconstruct an OpenMS-style experimental design from quantms.io parquet tables.

A quantms.io project is self-describing: ``run.parquet`` carries the run -> sample
mapping (including label, fraction and replicate numbers) and ``sample.parquet``
carries the sample characteristics. This module turns those two tables into the
same ``(sample_df, file_df)`` pair that :func:`read_openms_design` produces from an
external experimental-design TSV, so the sample-level sections of the report work
without one.

Target schema (matching ``read_openms_design``)::

    sample_df: Sample, Condition, BioReplicate
    file_df:   FractionGroup, Fraction, Label, Sample, Filename, Run
"""

from __future__ import absolute_import

import re

import pandas as pd

from pmultiqc.modules.common.common_utils import file_prefix
from pmultiqc.modules.common.logging import get_logger


log = get_logger("pmultiqc.modules.qpx.qpx_design")

# Isobaric channel -> numeric label, matching the ordering used elsewhere in pmultiqc.
CHANNEL_TO_LABEL = {
    "TMT126": 1, "TMT127N": 2, "TMT127C": 3, "TMT128N": 4, "TMT128C": 5,
    "TMT129N": 6, "TMT129C": 7, "TMT130N": 8, "TMT130C": 9, "TMT131": 10,
    "TMT131C": 11, "TMT132N": 12, "TMT132C": 13, "TMT133N": 14, "TMT133C": 15,
    "TMT134N": 16,
    "ITRAQ113": 1, "ITRAQ114": 2, "ITRAQ115": 3, "ITRAQ116": 4,
    "ITRAQ117": 5, "ITRAQ118": 6, "ITRAQ119": 7, "ITRAQ121": 8,
}

# Sample characteristics that never distinguish samples, so are never a Condition.
_NON_CONDITION_COLUMNS = {"sample_accession", "sample_description"}


def build_design_from_parquet(run_df, sample_df_raw):
    """Build ``(sample_df, file_df)`` from the run and sample parquet tables.

    Returns ``(None, None)`` when the tables cannot supply a usable design, so the
    caller can fall back to its existing behaviour rather than emit a broken table.
    """
    if run_df is None or getattr(run_df, "empty", True):
        return None, None
    if "samples" not in run_df.columns or "run_file_name" not in run_df.columns:
        log.debug("run.parquet lacks 'samples'/'run_file_name'; cannot derive a design")
        return None, None

    accessions = _sample_accessions(run_df, sample_df_raw)
    if not accessions:
        log.debug("no sample accessions found in run.parquet; cannot derive a design")
        return None, None
    sample_numbers = _number_samples(accessions)

    file_rows = []
    bio_replicates = {}
    for run in run_df.itertuples():
        run_file_name = getattr(run, "run_file_name", None)
        if not run_file_name:
            continue

        filename = getattr(run, "file_name", None) or run_file_name
        fraction = _as_int(getattr(run, "fraction", None), default=1)

        for entry in _iter_samples(getattr(run, "samples", None)):
            accession = entry.get("sample_accession")
            number = sample_numbers.get(accession)
            if number is None:
                continue

            bio_replicates.setdefault(
                number, _as_int(entry.get("biological_replicate"), default=1)
            )
            file_rows.append(
                {
                    # quantms.io has no explicit fraction group; runs of one sample
                    # form the group, so the sample number stands in for it.
                    "FractionGroup": number,
                    "Fraction": fraction,
                    "Label": _label_number(entry.get("label")),
                    "Sample": number,
                    "Filename": filename,
                    "Run": file_prefix(run_file_name),
                }
            )

    if not file_rows:
        log.debug("run.parquet produced no run/sample rows; cannot derive a design")
        return None, None

    file_df = pd.DataFrame(file_rows).drop_duplicates().reset_index(drop=True)

    conditions = _conditions(sample_df_raw, accessions)
    sample_df = pd.DataFrame(
        [
            {
                "Sample": number,
                "Condition": conditions.get(accession, "not available"),
                "BioReplicate": bio_replicates.get(number, 1),
            }
            for accession, number in sorted(sample_numbers.items(), key=lambda kv: kv[1])
            if number in set(file_df["Sample"])
        ]
    ).reset_index(drop=True)

    log.info(
        "Derived experimental design from quantms.io parquet: "
        f"{len(sample_df)} sample(s), {len(file_df)} run/sample row(s)."
    )
    return sample_df, file_df


def _iter_samples(cell):
    """Yield the per-sample dicts held in a run's 'samples' list column."""
    if cell is None:
        return
    try:
        entries = list(cell)
    except TypeError:
        return
    for entry in entries:
        if isinstance(entry, dict):
            yield entry


def _sample_accessions(run_df, sample_df_raw):
    """Collect sample accessions, preferring sample.parquet's order when available."""
    ordered = []
    if sample_df_raw is not None and not getattr(sample_df_raw, "empty", True):
        if "sample_accession" in sample_df_raw.columns:
            ordered = [str(a) for a in sample_df_raw["sample_accession"].dropna()]

    seen = set(ordered)
    for cell in run_df["samples"]:
        for entry in _iter_samples(cell):
            accession = entry.get("sample_accession")
            if accession is not None and str(accession) not in seen:
                seen.add(str(accession))
                ordered.append(str(accession))
    return ordered


def _number_samples(accessions):
    """Map each accession to the integer Sample id the report tables require.

    Downstream code does ``file_df["Sample"].astype(int)``, so a numeric id is
    mandatory. A trailing number in the accession ("Sample 3" -> 3) is used when it
    is unique across samples; otherwise samples are numbered by their order.
    """
    trailing = {}
    for accession in accessions:
        match = re.search(r"(\d+)\s*$", str(accession))
        trailing[accession] = int(match.group(1)) if match else None

    values = list(trailing.values())
    if all(v is not None for v in values) and len(set(values)) == len(values):
        return {a: int(v) for a, v in trailing.items()}

    return {accession: i for i, accession in enumerate(accessions, start=1)}


def _label_number(label):
    """Map a quantms.io label ('LFQ', 'TMT126', ...) to its numeric channel."""
    if label is None:
        return 1
    key = str(label).strip().upper().replace("-", "").replace("_", "")
    if key in CHANNEL_TO_LABEL:
        return CHANNEL_TO_LABEL[key]
    if key.isdigit():
        return int(key)
    # LFQ (and anything unrecognised) is a single channel.
    return 1


def _as_int(value, default):
    """Coerce a value to int, falling back to a default for null or malformed input."""
    try:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return default
        return int(value)
    except (TypeError, ValueError):
        return default


def _conditions(sample_df_raw, accessions):
    """Derive a Condition per sample from the characteristics that actually vary.

    quantms design files put a plain characteristic value in MSstats_Condition (e.g.
    the cell line). We mirror that: a single varying characteristic is used as-is,
    and several are emitted as ``key=value;key=value``, the multi-condition form
    ``draw_exp_design`` already understands.
    """
    if sample_df_raw is None or getattr(sample_df_raw, "empty", True):
        return {}
    if "sample_accession" not in sample_df_raw.columns:
        return {}

    df = sample_df_raw.drop_duplicates(subset=["sample_accession"])
    df = df[df["sample_accession"].notna()]
    if df.empty:
        return {}

    varying = [
        column
        for column in df.columns
        if column not in _NON_CONDITION_COLUMNS
        and df[column].notna().any()
        and df[column].astype(str).nunique() > 1
    ]

    if len(accessions) <= 1 or not varying:
        # Nothing distinguishes the samples; fall back to a constant characteristic
        # so the column still carries something meaningful.
        constant = [
            column
            for column in df.columns
            if column not in _NON_CONDITION_COLUMNS and df[column].notna().any()
        ]
        if not constant:
            return {}
        column = constant[0]
        return {
            str(r["sample_accession"]): str(r[column])
            for _, r in df.iterrows()
            if pd.notna(r[column])
        }

    conditions = {}
    for _, row in df.iterrows():
        if len(varying) == 1:
            value = row[varying[0]]
            conditions[str(row["sample_accession"])] = (
                str(value) if pd.notna(value) else "not available"
            )
        else:
            parts = [
                f"{column}={row[column]}"
                for column in varying
                if pd.notna(row[column]) and "=" not in str(row[column])
                and ";" not in str(row[column])
            ]
            conditions[str(row["sample_accession"])] = (
                ";".join(parts) if parts else "not available"
            )
    return conditions
