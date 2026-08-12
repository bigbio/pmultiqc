"""QC heatmap scores for the QPX module, computed from the quantms.io parquet tables.

Mirrors the DDA heatmap that quantms/maxquant/mzIdentML already draw, using the same
formulas so scores are comparable across report types. Every metric is optional: one
whose input column is absent or empty is dropped from the heatmap rather than emitted
as a misleading flat row.
"""

from __future__ import absolute_import

import numpy as np
import pandas as pd

from pmultiqc.modules.common.logging import get_logger
from pmultiqc.modules.common.stats import cal_hm_charge, nanmedian, qual_uniform
from pmultiqc.modules.qpx.qpx_sections import resolve_contaminant_flags


log = get_logger("pmultiqc.modules.qpx.qpx_heatmap")

# quantms scales the median intensity against this reference before clamping to 1.
INTENSITY_REFERENCE = 2 ** 23

# Display order, matching the quantms DDA heatmap.
METRIC_ORDER = [
    "Contaminants",
    "Peptide Intensity",
    "Charge",
    "Missed Cleavages",
    "Missed Cleavages Var",
    "ID rate over RT",
    "MS2 OverSampling",
    "Pep Missing Values",
]


def calculate_qpx_heatmap(psm_df, pg_df, feature_df, missed_cleavages_by_run, contaminant_affix):
    """Build ``(heat_map_score, xnames, ynames)`` for :func:`draw_heatmap`.

    Returns ``([], [], [])`` when nothing can be computed, so the caller can skip the
    plot entirely.
    """
    if psm_df is None or getattr(psm_df, "empty", True) or "run" not in psm_df.columns:
        return [], [], []

    runs = sorted({str(r) for r in psm_df["run"].dropna()})
    if not runs:
        return [], [], []

    scores = {
        "Contaminants": _contaminant_scores(pg_df, contaminant_affix),
        "Peptide Intensity": _peptide_intensity_scores(feature_df),
        "Charge": _charge_scores(psm_df),
        "Missed Cleavages": _missed_cleavage_scores(missed_cleavages_by_run),
        "ID rate over RT": _id_rate_over_rt_scores(psm_df),
        "MS2 OverSampling": _oversampling_scores(psm_df),
        "Pep Missing Values": _pep_missing_value_scores(psm_df),
    }
    scores["Missed Cleavages Var"] = _variance_scores(scores["Missed Cleavages"])

    # Keep only metrics that cover every run, so no cell is silently invented.
    xnames = [
        metric
        for metric in METRIC_ORDER
        if scores.get(metric) and all(run in scores[metric] for run in runs)
    ]
    dropped = [
        metric
        for metric in METRIC_ORDER
        if metric not in xnames
    ]
    if dropped:
        log.info(
            "[HeatMap] Metrics not available for this project and omitted: "
            f"{', '.join(dropped)}."
        )
    if not xnames:
        log.warning("[HeatMap] No metrics could be computed; skipping the heatmap.")
        return [], [], []

    heat_map_score = [[float(scores[metric][run]) for metric in xnames] for run in runs]
    return heat_map_score, xnames, runs


def _contaminant_scores(pg_df, contaminant_affix):
    """1 - (contaminant intensity / total intensity) per run.

    Prefers the explicit ``contaminant`` flag, falling back to matching the affix
    against the protein accession the way the other plugins do.
    """
    if pg_df is None or getattr(pg_df, "empty", True):
        return {}
    if "run" not in pg_df.columns or "intensity" not in pg_df.columns:
        return {}

    df = pg_df[["run", "intensity"]].copy()
    df["intensity"] = pd.to_numeric(df["intensity"], errors="coerce")

    # Resolved by the shared helper so this metric and the Contaminants section can
    # never disagree: the producer flag unioned with a case-insensitive affix match.
    flags = resolve_contaminant_flags(pg_df, contaminant_affix)
    if flags is None:
        log.info(
            "[HeatMap] Contaminants omitted: the 'contaminant' flag is entirely null "
            "(unknown for this producer) and no accession column matches the affix."
        )
        return {}
    df["is_contaminant"] = flags

    df = df.dropna(subset=["intensity"])
    if df.empty:
        return {}

    result = {}
    for run, group in df.groupby("run"):
        total = group["intensity"].sum()
        if total <= 0:
            continue
        contaminated = group.loc[group["is_contaminant"], "intensity"].sum()
        result[str(run)] = float(max(0.0, 1.0 - contaminated / total))
    return result


def _accessions_match(accessions, affix):
    """Match an affix against a list<string> accession column."""
    def matches(cell):
        if cell is None:
            return False
        try:
            return any(affix in str(a) for a in cell)
        except TypeError:
            return affix in str(cell)

    return accessions.apply(matches)


def _peptide_intensity_scores(feature_df):
    """min(1, median(feature intensity) / 2**23) per run."""
    if feature_df is None or getattr(feature_df, "empty", True):
        return {}
    if "run" not in feature_df.columns or "intensities" not in feature_df.columns:
        return {}

    exploded = feature_df[["run", "intensities"]].dropna(subset=["intensities"])
    exploded = exploded.explode("intensities").dropna(subset=["intensities"])
    if exploded.empty:
        return {}

    exploded["intensity"] = pd.to_numeric(
        exploded["intensities"].str.get("intensity"), errors="coerce"
    )
    exploded = exploded[exploded["intensity"] > 0]
    if exploded.empty:
        return {}

    result = {}
    for run, group in exploded.groupby("run"):
        median = nanmedian(group["intensity"].to_numpy(), np.float64(0.0))
        result[str(run)] = float(min(1.0, median / INTENSITY_REFERENCE))
    return result


def _charge_scores(psm_df):
    """Charge metric: deviation of each run's charge-2 fraction from the median."""
    if "charge" not in psm_df.columns or psm_df["charge"].isna().all():
        return {}
    return {str(k): float(v) for k, v in cal_hm_charge(psm_df, "run", "charge").items()}


def _missed_cleavage_scores(missed_cleavages_by_run):
    """Fraction of PSMs with zero missed cleavages, per run."""
    if not missed_cleavages_by_run:
        return {}

    result = {}
    for run, counts in missed_cleavages_by_run.items():
        if not isinstance(counts, dict):
            continue
        total = sum(counts.values())
        if total <= 0:
            continue
        # Keys may arrive as int or str depending on how the counts were built.
        zero = counts.get(0, counts.get("0", 0))
        result[str(run)] = float(zero / total)
    return result


def _variance_scores(base_scores):
    """1 - |value - median| across runs, as used for Missed Cleavages Var."""
    if not base_scores:
        return {}
    median = np.median(list(base_scores.values()))
    return {
        run: float(max(0.0, 1.0 - abs(value - median)))
        for run, value in base_scores.items()
    }


def _id_rate_over_rt_scores(psm_df):
    """ID-rate-over-RT metric: how uniformly identifications spread over the gradient."""
    if "rt" not in psm_df.columns or psm_df["rt"].isna().all():
        return {}
    return {
        str(run): float(qual_uniform(group["rt"]))
        for run, group in psm_df[["run", "rt"]].groupby("run")
    }


def _oversampling_scores(psm_df):
    """Fraction of precursors identified exactly once (MS2 oversampling)."""
    if "peptidoform" not in psm_df.columns or "charge" not in psm_df.columns:
        return {}

    df = psm_df[["run", "peptidoform", "charge"]].dropna(subset=["peptidoform"])
    if df.empty:
        return {}

    result = {}
    for run, group in df.groupby("run"):
        counts = group.groupby(["peptidoform", "charge"]).size()
        if counts.empty:
            continue
        result[str(run)] = float(min(1.0, (counts == 1).sum() / len(counts)))
    return result


def _pep_missing_value_scores(psm_df):
    """Share of the project's peptidoforms that each run identifies."""
    if "peptidoform" not in psm_df.columns:
        return {}

    df = psm_df[["run", "peptidoform"]].dropna(subset=["peptidoform"])
    if df.empty:
        return {}

    global_peptidoforms = set(df["peptidoform"])
    if not global_peptidoforms:
        return {}

    return {
        str(run): float(
            min(1.0, len(set(group["peptidoform"]) & global_peptidoforms)
                / len(global_peptidoforms))
        )
        for run, group in df.groupby("run")
    }
