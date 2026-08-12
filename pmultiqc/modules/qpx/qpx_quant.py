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
    if "anchor_protein" not in pg_df.columns or "intensity" not in pg_df.columns:
        log.debug("pg.parquet lacks anchor_protein/intensity; no protein table.")
        return None, None

    df = pg_df.copy()
    df["intensity"] = pd.to_numeric(df["intensity"], errors="coerce")
    df = df[df["intensity"] > 0].dropna(subset=["anchor_protein"])
    if df.empty:
        log.debug("pg.parquet has no positive intensities; no protein table.")
        return None, None

    peptides_per_protein = _peptides_per_protein(feature_df)

    table_dict = {}
    for protein, group in df.groupby("anchor_protein"):
        protein = str(protein)
        entry = {
            "ProteinName": protein,
            "Average Intensity": float(np.log10(group["intensity"].mean())),
        }
        if protein in peptides_per_protein:
            entry["Peptides_Number"] = int(peptides_per_protein[protein])
        if "global_qvalue" in group.columns and group["global_qvalue"].notna().any():
            entry["Global Q-value"] = float(group["global_qvalue"].min())
        table_dict[protein] = entry

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


def _peptides_per_protein(feature_df):
    """Count distinct peptide sequences per anchor protein, if features are available."""
    if feature_df is None or getattr(feature_df, "empty", True):
        return {}
    if "anchor_protein" not in feature_df.columns or "sequence" not in feature_df.columns:
        return {}

    counts = (
        feature_df.dropna(subset=["anchor_protein", "sequence"])
        .groupby("anchor_protein")["sequence"]
        .nunique()
    )
    return {str(k): v for k, v in counts.items()}


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
        for protein, group in condition_group.groupby("anchor_protein"):
            entry = table_dict.get(str(protein))
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
