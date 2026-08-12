"""Protein-level quantification for the QPX module.

``pg.parquet`` carries an intensity and a q-value per (protein group, run), but the
QPX module previously drew nothing from it -- its Quantification Analysis section held
only the peptide intensity box. This module turns that table into the same Protein
Quantification Table the DIA-NN and mzIdentML reports draw, plus a protein intensity
distribution matching the peptide one.
"""

from __future__ import absolute_import

import numpy as np
import pandas as pd

from multiqc.plots import box

from pmultiqc.modules.common.logging import get_logger
from pmultiqc.modules.common.plots.general import (
    plot_html_check,
    stat_pep_intensity,
)
from pmultiqc.modules.core.section_groups import add_sub_section


log = get_logger("pmultiqc.modules.qpx.qpx_quant")

# Matches the row cap used by the DIA-NN/mzIdentML protein tables.
MAX_TABLE_ROWS = 50


def create_qpx_protein_table(pg_df, feature_df, sample_df, file_df):
    """Build ``(table_dict, headers)`` for :func:`draw_protein_table`.

    Returns ``(None, None)`` when pg.parquet cannot supply protein quantification.
    """
    if pg_df is None or getattr(pg_df, "empty", True):
        return None, None
    if "intensity" not in pg_df.columns:
        log.debug("pg.parquet lacks intensity; no protein table.")
        return None, None

    df = pg_df.copy()
    df["intensity"] = pd.to_numeric(df["intensity"], errors="coerce")
    df = df[df["intensity"] > 0]

    # The spec is explicit that anchor_protein is a display label only and must not key
    # per-protein statistics -- protein groups are not guaranteed disjoint, so several
    # groups can share an anchor. Group on the full accession set instead and keep
    # anchor_protein purely for the displayed name.
    df["_group_key"] = protein_group_key(df)
    df = df.dropna(subset=["_group_key"])
    if df.empty:
        log.debug("pg.parquet has no positive intensities; no protein table.")
        return None, None

    peptides_per_protein = _peptides_per_protein(feature_df)

    table_dict = {}
    for group_key, group in df.groupby("_group_key"):
        protein = _display_name(group)
        entry = {
            "ProteinName": protein,
            "Average Intensity": float(np.log10(group["intensity"].mean())),
        }
        if group_key in peptides_per_protein:
            entry["Peptides_Number"] = int(peptides_per_protein[group_key])
        if "global_qvalue" in group.columns and group["global_qvalue"].notna().any():
            entry["Global Q-value"] = float(group["global_qvalue"].min())
        table_dict[group_key] = entry

    headers = {
        "ProteinName": {
            "title": "Protein Name",
            "description": "Name/Identifier(s) of the protein (group)",
        },
    }
    if peptides_per_protein:
        headers["Peptides_Number"] = {
            "title": "Number of Peptides",
            "description": "Number of distinct peptide sequences per protein",
            "format": "{:,.0f}",
        }
    headers["Average Intensity"] = {
        "title": "Average Intensity",
        "description": "log10 of the average protein group intensity across runs",
        "format": "{:,.4f}",
    }
    if any("Global Q-value" in entry for entry in table_dict.values()):
        headers["Global Q-value"] = {
            "title": "Global Q-value",
            "description": "Best (lowest) global q-value reported for the protein group",
            "format": "{:,.2e}",
        }

    condition_columns = _add_condition_intensities(table_dict, df, sample_df, file_df)
    for condition in condition_columns:
        headers[condition] = {
            "title": condition,
            "description": f"log10 average protein intensity in condition '{condition}'",
            "format": "{:,.4f}",
        }

    # Sort by intensity so the capped table shows the most abundant proteins.
    ordered = sorted(
        table_dict.values(), key=lambda e: e.get("Average Intensity", 0), reverse=True
    )
    result_dict = {i: entry for i, entry in enumerate(ordered, start=1)}

    log.info(
        f"[Quantification] Protein quantification table: {len(table_dict)} protein group(s)"
        f"{f', showing top {MAX_TABLE_ROWS}' if len(table_dict) > MAX_TABLE_ROWS else ''}."
    )
    return result_dict, headers


def _accession_key(cell):
    """Stable key for a list<string> accession set, or None when unusable."""
    if cell is None:
        return None
    if isinstance(cell, str):
        return cell or None
    try:
        parts = sorted({str(a) for a in cell if a is not None and str(a) != ""})
    except TypeError:
        text = str(cell)
        return text or None
    return ";".join(parts) or None


def protein_group_key(df):
    """Key protein groups by their accession set, never by anchor_protein.

    Protein groups are not guaranteed disjoint and anchor_protein is documented as a
    display label, so two distinct groups can share one anchor. pg_accessions is the
    identifying set; anchor_protein is only a fallback when it is unavailable.
    """
    if "pg_accessions" in df.columns:
        keys = df["pg_accessions"].apply(_accession_key)
        if keys.notna().any():
            return keys
    if "anchor_protein" in df.columns:
        return df["anchor_protein"].astype(str).where(df["anchor_protein"].notna())
    return pd.Series([None] * len(df), index=df.index)


def _display_name(group):
    """Human-readable label for a protein group: the anchor when present."""
    if "anchor_protein" in group.columns and group["anchor_protein"].notna().any():
        return str(group["anchor_protein"].dropna().iloc[0])
    return str(group.name if hasattr(group, "name") else group.index[0])


def _peptides_per_protein(feature_df):
    """Count distinct peptide sequences per protein group, keyed like the pg table."""
    if feature_df is None or getattr(feature_df, "empty", True):
        return {}
    if "sequence" not in feature_df.columns:
        return {}

    df = feature_df.copy()
    # feature.pg_accessions is a list<struct> (accession/start/end); reduce to accessions.
    if "pg_accessions" in df.columns:
        df["_group_key"] = df["pg_accessions"].apply(_feature_accession_key)
    elif "anchor_protein" in df.columns:
        df["_group_key"] = df["anchor_protein"].astype(str)
    else:
        return {}

    df = df.dropna(subset=["_group_key", "sequence"])
    if df.empty:
        return {}

    counts = df.groupby("_group_key")["sequence"].nunique()
    return {str(k): v for k, v in counts.items()}


def _feature_accession_key(cell):
    """feature.pg_accessions holds structs; pull out the accession strings."""
    if cell is None:
        return None
    try:
        accessions = []
        for entry in cell:
            if isinstance(entry, dict):
                value = entry.get("accession")
            else:
                value = entry
            if value is not None and str(value) != "":
                accessions.append(str(value))
    except TypeError:
        return _accession_key(cell)
    return ";".join(sorted(set(accessions))) or None


def _add_condition_intensities(table_dict, df, sample_df, file_df):
    """Add a per-condition average intensity column, when a design is available."""
    if sample_df is None or getattr(sample_df, "empty", True):
        return []
    if file_df is None or getattr(file_df, "empty", True):
        return []
    if "run" not in df.columns:
        return []
    if not {"Sample", "Condition"}.issubset(sample_df.columns):
        return []
    if not {"Sample", "Run"}.issubset(file_df.columns):
        return []

    run_to_sample = file_df[["Sample", "Run"]].drop_duplicates()
    sample_to_condition = sample_df[["Sample", "Condition"]].drop_duplicates()

    merged = df.merge(run_to_sample, left_on="run", right_on="Run", how="inner")
    merged = merged.merge(sample_to_condition, on="Sample", how="inner")
    if merged.empty:
        return []

    conditions = []
    for condition, condition_group in merged.groupby("Condition"):
        condition = str(condition)
        conditions.append(condition)
        for group_key, group in condition_group.groupby("_group_key"):
            entry = table_dict.get(group_key)
            if entry is not None:
                entry[condition] = float(np.log10(group["intensity"].mean()))

    return conditions


def get_protein_intensity(pg_df, sdrf_file_df):
    """Protein group intensity distributions as ``[by_run, by_sample]``.

    Mirrors :func:`get_pep_intensity` so the protein plot matches the peptide one.
    """
    by_run = {}
    by_sample = {}

    if pg_df is None or getattr(pg_df, "empty", True):
        return [by_run, by_sample]
    if "run" not in pg_df.columns or "intensity" not in pg_df.columns:
        return [by_run, by_sample]

    df = pg_df[["run", "intensity"]].copy()
    df["intensity"] = pd.to_numeric(df["intensity"], errors="coerce")
    df = df[df["intensity"] > 0]
    if df.empty:
        return [by_run, by_sample]

    for run, group in df.groupby("run"):
        by_run[str(run)] = stat_pep_intensity(group["intensity"])

    if sdrf_file_df is not None and not getattr(sdrf_file_df, "empty", True):
        if {"Sample", "Run"}.issubset(sdrf_file_df.columns):
            merged = df.merge(
                right=sdrf_file_df[["Sample", "Run"]].drop_duplicates(),
                left_on="run",
                right_on="Run",
                how="left",
            ).dropna(subset=["Sample"])

            for sample, group in merged.groupby("Sample", sort=True):
                by_sample[f"Sample {str(sample)}"] = stat_pep_intensity(group["intensity"])

    return [by_run, by_sample]


def draw_protein_intensity(sub_section, plot_data):
    """Box plot of protein group intensity, by run and (when known) by sample."""
    if not plot_data or not plot_data[0]:
        return

    draw_config = {
        "id": "protein_intensity_distribution_box",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Protein Intensity Distribution",
        "tt_decimals": 2,
        "xlab": "log2(Intensity)",
        "sort_samples": False,
        "save_data_file": False,
    }

    if len(plot_data) > 1 and plot_data[1]:
        draw_config["data_labels"] = ["by Run", "by Sample"]
    else:
        plot_data = [plot_data[0]]

    box_html = plot_html_check(box.plot(plot_data, pconfig=draw_config))

    add_sub_section(
        sub_section=sub_section,
        plot=box_html,
        order=6,
        description="Protein group intensity per Run.",
        helptext="""
            The intensity of each protein group is taken from the `intensity` column of
            pg.parquet and log2-transformed. Zero and missing intensities are ignored.
            """,
    )


def create_qpx_peptide_table(feature_df, sample_df, file_df):
    """Peptide-level quantification table, mirroring the protein one.

    Returns ``(table_dict, headers)`` or ``(None, None)``.
    """
    from pmultiqc.modules.qpx.qpx_io import has_data

    if not has_data(feature_df, "peptidoform", "intensities"):
        log.debug("feature.parquet lacks peptidoform/intensities; no peptide table.")
        return None, None

    df = feature_df.dropna(subset=["intensities"]).explode("intensities")
    df = df.dropna(subset=["intensities"])
    if df.empty:
        return None, None

    df["intensity"] = pd.to_numeric(df["intensities"].str.get("intensity"), errors="coerce")
    df = df[df["intensity"] > 0]
    if df.empty:
        return None, None

    df["_group_key"] = protein_group_key(df) if "pg_accessions" in df.columns else None

    table_dict = {}
    for peptidoform, group in df.groupby("peptidoform"):
        entry = {
            "PeptideID": str(peptidoform),
            "Average Intensity": float(np.log10(group["intensity"].mean())),
        }
        if "anchor_protein" in group.columns and group["anchor_protein"].notna().any():
            entry["ProteinName"] = str(group["anchor_protein"].dropna().iloc[0])
        if "charge" in group.columns and group["charge"].notna().any():
            entry["Charge"] = int(group["charge"].dropna().iloc[0])
        table_dict[str(peptidoform)] = entry

    headers = {"PeptideID": {"title": "Peptidoform", "description": "Peptide sequence with modifications"}}
    if any("ProteinName" in e for e in table_dict.values()):
        headers["ProteinName"] = {"title": "Protein Name", "description": "Anchor protein of the group"}
    if any("Charge" in e for e in table_dict.values()):
        headers["Charge"] = {"title": "Charge", "description": "Precursor charge", "format": "{:,.0f}"}
    headers["Average Intensity"] = {
        "title": "Average Intensity",
        "description": "log10 of the average peptide intensity across runs",
        "format": "{:,.4f}",
    }

    ordered = sorted(table_dict.values(), key=lambda e: e.get("Average Intensity", 0), reverse=True)
    result = {i: entry for i, entry in enumerate(ordered, start=1)}
    log.info(f"[Quantification] Peptide quantification table: {len(table_dict)} peptidoform(s).")
    return result, headers


def protein_intensity_pca(pg_df):
    """PCA of the protein-group x run intensity matrix. None if not computable."""
    try:
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
    except ImportError:
        log.debug("[Quantification] scikit-learn unavailable; skipping PCA.")
        return None

    from pmultiqc.modules.qpx.qpx_io import has_data

    if not has_data(pg_df, "run", "intensity"):
        return None

    df = pg_df.copy()
    df["intensity"] = pd.to_numeric(df["intensity"], errors="coerce")
    df["_group_key"] = protein_group_key(df)
    df = df.dropna(subset=["intensity", "_group_key"])
    df = df[df["intensity"] > 0]
    if df.empty:
        return None

    matrix = df.pivot_table(
        index="_group_key", columns="run", values="intensity", aggfunc="mean"
    )
    # PCA needs at least two runs, and rows complete across them.
    matrix = matrix.dropna()
    if matrix.shape[1] < 2 or matrix.shape[0] < 2:
        log.info(
            "[Quantification] PCA skipped: needs >=2 runs and >=2 proteins quantified in "
            f"all of them (got {matrix.shape[0]} x {matrix.shape[1]})."
        )
        return None

    scaled = StandardScaler().fit_transform(np.log2(matrix).T)
    result = PCA(n_components=2).fit_transform(scaled)

    log.info(
        f"[Quantification] PCA over {matrix.shape[0]} protein groups x {matrix.shape[1]} runs."
    )
    return {
        str(run): {"x": float(result[i, 0]), "y": float(result[i, 1])}
        for i, run in enumerate(matrix.columns)
    }


def calculate_intensity_std(feature_df, sample_df, file_df):
    """Per-condition standard deviation of log2 peptide intensity."""
    from pmultiqc.modules.qpx.qpx_io import has_data

    if not has_data(feature_df, "peptidoform", "intensities", "run"):
        return {}
    if sample_df is None or getattr(sample_df, "empty", True):
        return {}
    if file_df is None or getattr(file_df, "empty", True):
        return {}
    if not {"Sample", "Condition"}.issubset(sample_df.columns):
        return {}
    if not {"Sample", "Run"}.issubset(file_df.columns):
        return {}

    df = feature_df[["run", "peptidoform", "intensities"]].dropna(subset=["intensities"])
    df = df.explode("intensities").dropna(subset=["intensities"])
    df["intensity"] = pd.to_numeric(df["intensities"].str.get("intensity"), errors="coerce")
    df = df[df["intensity"] > 0]
    if df.empty:
        return {}

    df = df.merge(
        file_df[["Sample", "Run"]].drop_duplicates(), left_on="run", right_on="Run", how="inner"
    ).merge(sample_df[["Sample", "Condition"]].drop_duplicates(), on="Sample", how="inner")
    if df.empty:
        return {}

    df["log2_intensity"] = np.log2(df["intensity"])

    result = {}
    for condition, group in df.groupby("Condition"):
        stds = group.groupby("peptidoform")["log2_intensity"].std().dropna()
        if not stds.empty:
            result[str(condition)] = stds.tolist()

    return result


def draw_intensity_std(sub_section, box_data):
    """Box plot of per-condition standard deviation of log2 peptide intensity."""
    if not box_data:
        return

    draw_config = {
        "id": "qpx_std_intensity_box",
        "title": "Standard Deviation of Intensity",
        "cpswitch": False,
        "tt_decimals": 5,
        "xlab": "Standard Deviation of log2(Intensity)",
        "save_data_file": False,
    }

    box_html = plot_html_check(box.plot(list_of_data_by_sample=box_data, pconfig=draw_config))

    add_sub_section(
        sub_section=sub_section,
        plot=box_html,
        order=8,
        description="Spread of each peptide's log2 intensity within a condition.",
        helptext="""
            For every peptidoform the standard deviation of its log2 intensity is computed
            across the runs of one condition; the box shows the distribution of those values.
            A tight distribution indicates reproducible quantification within the condition.
            """,
    )
