"""Identification, contaminant and search-engine sections for the QPX module.

Each builder returns None (or an empty structure) when the QPX fields it needs are
absent or entirely null, so a section appears only for producers that actually supply
its data. See :func:`pmultiqc.modules.qpx.qpx_io.has_data` for why presence alone is
not enough.
"""

from __future__ import absolute_import

from collections import OrderedDict

import pandas as pd

from pmultiqc.modules.common.common_utils import cal_contaminant_percent
from pmultiqc.modules.common.histogram import Histogram
from pmultiqc.modules.common.logging import get_logger
from pmultiqc.modules.qpx.qpx_io import has_data
from pmultiqc.modules.qpx.qpx_quant import protein_group_key


log = get_logger("pmultiqc.modules.qpx.qpx_sections")

TOP_N_CONTAMINANTS = 5

# Sentinels for the contaminant label column. cal_contaminant_percent flags rows by
# substring, so the "clean" sentinel must not contain the contaminant marker -- naming
# them CONTAMINANT / NOT_CONTAMINANT would match every row and report 100%.
CONTAMINANT_MARKER = "__CONTAM__"
CLEAN_LABEL = "clean"


# ---------------------------------------------------------------- peptides per protein

def build_peptides_per_protein(feature_df, peptides_per_group):
    """Histogram of distinct peptides per protein group, as draw_num_pep_per_protein wants."""
    if not peptides_per_group:
        log.debug("[Identification] No peptides-per-protein counts available.")
        return None

    pep_plot = Histogram("number of peptides per proteins", plot_category="frequency")
    for count in peptides_per_group.values():
        pep_plot.add_value(int(count))

    categories = OrderedDict()
    categories["Frequency"] = {
        "name": "Frequency",
        "description": "number of peptides per proteins",
    }
    pep_plot.to_dict(percentage=True, cats=categories)

    return pep_plot


# ------------------------------------------------------------------------ oversampling

def calculate_oversampling(id_df):
    """MS/MS counts per 3D-peak: how often each precursor was identified per run.

    Returns ``(plot_data, cats)`` or ``(None, None)``. A precursor is a
    (peptidoform, charge) pair, matching the definition the other plugins use.
    """
    if not has_data(id_df, "run", "peptidoform", "charge"):
        log.debug("[MS2] Oversampling needs run/peptidoform/charge; skipping.")
        return None, None

    df = id_df[["run", "peptidoform", "charge"]].dropna(subset=["peptidoform", "charge"])
    if df.empty:
        return None, None

    plot_data = {}
    seen_bins = set()
    for run, group in df.groupby("run"):
        counts = group.groupby(["peptidoform", "charge"]).size()
        if counts.empty:
            continue
        # Bin as 1, 2, >=3, the convention used by the MaxQuant/FragPipe plots.
        binned = counts.clip(upper=3).value_counts()
        total = binned.sum()
        if total <= 0:
            continue
        run_data = {}
        for value, count in binned.items():
            label = ">=3" if value >= 3 else str(int(value))
            run_data[label] = float(count / total * 100)
            seen_bins.add(label)
        plot_data[str(run)] = run_data

    if not plot_data:
        return None, None

    cats = OrderedDict()
    for label in ["1", "2", ">=3"]:
        if label in seen_bins:
            cats[label] = {"name": label}

    return plot_data, cats


# ------------------------------------------------------------------------ contaminants

def calculate_contaminants(pg_df, contaminant_affix):
    """Per-run contaminant intensity share, plus the top-N contaminant breakdown.

    Returns ``(per_run_percent, top_n_data)``; either may be None. Distinguishes an
    all-null 'contaminant' flag (unknown -> skip) from an all-False one (known clean).
    """
    if pg_df is None or getattr(pg_df, "empty", True):
        return None, None
    if not has_data(pg_df, "run", "intensity"):
        log.debug("[Contaminants] Needs run/intensity in pg.parquet; skipping.")
        return None, None

    df = pg_df.copy()
    df["intensity"] = pd.to_numeric(df["intensity"], errors="coerce")
    df = df.dropna(subset=["intensity"])
    if df.empty:
        return None, None

    labels = _contaminant_labels(df, contaminant_affix)
    if labels is None:
        log.info(
            "[Contaminants] Skipped: the 'contaminant' flag is entirely null for this "
            "producer and no accession column is available to match the affix."
        )
        return None, None

    df["_contaminant"] = labels

    per_run = cal_contaminant_percent(
        df=df.assign(_protein=df["_contaminant"]),
        protein_col="_protein",
        intensity_col="intensity",
        run_col="run",
        contam_affix=CONTAMINANT_MARKER,
    )

    # An all-zero series means the producer flagged no contaminants at all. That is a
    # real answer, but a bar chart of zeros shows nothing (and MultiQC rejects it), so
    # report it in the log and skip the section.
    if per_run and not any(
        v.get("Potential Contaminants", 0) > 0 for v in per_run.values()
    ):
        log.info(
            "[Contaminants] No contaminant signal in any run (0% throughout); "
            "skipping the contaminants section."
        )
        return None, None

    top_n = _top_n_contaminants(df)
    return per_run, top_n


def _contaminant_labels(df, affix):
    """Label each row 'CONTAMINANT:<name>' or 'NOT_CONTAMINANT'; None if undeterminable."""
    names = _protein_names(df)

    if "contaminant" in df.columns and df["contaminant"].notna().any():
        flags = df["contaminant"].fillna(False).astype(bool)
    elif affix:
        if "pg_accessions" in df.columns:
            flags = df["pg_accessions"].apply(
                lambda cell: _cell_contains(cell, affix)
            )
        elif "anchor_protein" in df.columns:
            flags = df["anchor_protein"].astype(str).str.contains(affix, regex=False)
        else:
            return None
    else:
        return None

    return pd.Series(
        [CONTAMINANT_MARKER + n if f else CLEAN_LABEL for n, f in zip(names, flags)],
        index=df.index,
    )


def _protein_names(df):
    if "anchor_protein" in df.columns:
        return df["anchor_protein"].astype(str).fillna("unknown")
    return protein_group_key(df).astype(str)


def _cell_contains(cell, affix):
    if cell is None:
        return False
    try:
        return any(affix in str(a) for a in cell)
    except TypeError:
        return affix in str(cell)


def _top_n_contaminants(df):
    """Per-run intensity share of the N most abundant contaminants."""
    contaminated = df[df["_contaminant"] != CLEAN_LABEL]
    if contaminated.empty:
        return None

    totals = df.groupby("run")["intensity"].sum()
    ranked = (
        contaminated.groupby("_contaminant")["intensity"].sum().sort_values(ascending=False)
    )
    top = list(ranked.index[:TOP_N_CONTAMINANTS])

    result = {}
    for run, group in df.groupby("run"):
        total = totals.get(run, 0)
        if total <= 0:
            continue
        run_data = {}
        for name in top:
            share = group.loc[group["_contaminant"] == name, "intensity"].sum()
            run_data[name.replace(CONTAMINANT_MARKER, "")] = float(share / total * 100)
        other = group.loc[
            (group["_contaminant"] != CLEAN_LABEL)
            & (~group["_contaminant"].isin(top)),
            "intensity",
        ].sum()
        if other > 0:
            run_data["Other"] = float(other / total * 100)
        result[str(run)] = run_data

    return result or None


# ------------------------------------------------------------- search engine scores

def calculate_search_engine_scores(id_df):
    """Per-run distribution of a search-engine score, binned adaptively.

    QPX score names are an open, free-form set (snake_case by convention), and their
    ranges differ wildly -- a q-value or posterior error probability lives in [0, 1]
    while a raw engine score can reach the hundreds. So the score is chosen from what
    the project actually carries and the bins are derived from its observed range,
    rather than assuming MaxQuant's fixed 0-300 Andromeda bins.

    Returns ``(plot_data, score_name)`` or ``(None, None)``.
    """
    if id_df is None or getattr(id_df, "empty", True) or "run" not in id_df.columns:
        return None, None

    scores, score_name = _pick_score(id_df)
    if scores is None:
        log.debug("[Search Engine Scores] No usable score column; skipping.")
        return None, None

    frame = pd.DataFrame({"run": id_df["run"], "score": scores}).dropna(subset=["score"])
    frame = frame[pd.to_numeric(frame["score"], errors="coerce").notna()]
    if frame.empty:
        return None, None

    edges = _bin_edges(frame["score"])
    if edges is None:
        log.debug(f"[Search Engine Scores] '{score_name}' has no spread; skipping.")
        return None, None

    labels = [f"{edges[i]:.4g} ~ {edges[i + 1]:.4g}" for i in range(len(edges) - 1)]
    frame["bin"] = pd.cut(
        frame["score"], bins=edges, labels=labels, include_lowest=True, right=False
    )

    plot_data = {}
    for run, group in frame.groupby("run"):
        counts = group["bin"].value_counts()
        plot_data[str(run)] = {label: int(counts.get(label, 0)) for label in labels}

    log.info(
        f"[Search Engine Scores] Using '{score_name}' over {len(labels)} bins "
        f"across {len(plot_data)} run(s)."
    )
    return plot_data, score_name


# Preferred in order: a real engine score, then generic confidence measures.
_SCORE_PREFERENCE = ["andromeda_score", "hyperscore", "xcorr", "sage_hyperscore",
                     "msgf_raw_score", "comet_xcorr", "diann_cscore"]


def _pick_score(id_df):
    """Choose a score column: a named engine score if present, else PEP, else q-value."""
    if has_data(id_df, "additional_scores"):
        available = _score_names(id_df["additional_scores"])
        for preferred in _SCORE_PREFERENCE:
            for name in available:
                if name.lower() == preferred:
                    return _extract_score(id_df["additional_scores"], name), name

    if has_data(id_df, "posterior_error_probability"):
        return pd.to_numeric(
            id_df["posterior_error_probability"], errors="coerce"
        ), "posterior_error_probability"

    if has_data(id_df, "additional_scores"):
        available = _score_names(id_df["additional_scores"])
        for name in available:
            if "qvalue" in name.lower().replace("-", "").replace("_", ""):
                return _extract_score(id_df["additional_scores"], name), name
        if available:
            name = sorted(available)[0]
            return _extract_score(id_df["additional_scores"], name), name

    return None, None


def _score_names(series, sample_size=5000):
    names = set()
    for cell in series.head(sample_size):
        if cell is None:
            continue
        try:
            for entry in cell:
                if isinstance(entry, dict) and entry.get("score_name"):
                    names.add(str(entry["score_name"]))
        except TypeError:
            continue
    return names


def _extract_score(series, name):
    def pull(cell):
        if cell is None:
            return None
        try:
            for entry in cell:
                if isinstance(entry, dict) and entry.get("score_name") == name:
                    return entry.get("score_value")
        except TypeError:
            return None
        return None

    return pd.to_numeric(series.apply(pull), errors="coerce")


def _bin_edges(scores, target_bins=20):
    """Bin edges spanning the observed range, clipped to the 1st-99th percentile."""
    numeric = pd.to_numeric(scores, errors="coerce").dropna()
    if numeric.empty or numeric.nunique() < 2:
        return None

    low = float(numeric.quantile(0.01))
    high = float(numeric.quantile(0.99))
    if high <= low:
        low, high = float(numeric.min()), float(numeric.max())
    if high <= low:
        return None

    step = (high - low) / target_bins
    return [low + step * i for i in range(target_bins + 1)]


def draw_qpx_search_engine_scores(sub_section, plot_data, score_name):
    """Bar plot of the chosen search-engine score distribution per run.

    A dedicated drawer rather than general.draw_search_engine_scores, because that one
    hard-codes a title per tool ("Andromeda", "Hyperscore") while the QPX score is
    chosen at runtime from whatever the project carries.
    """
    from multiqc.plots import bargraph

    from pmultiqc.modules.common.plots.general import plot_html_check
    from pmultiqc.modules.core.section_groups import add_sub_section

    if not plot_data:
        return

    pretty = score_name.replace("_", " ").title()
    draw_config = {
        "id": "qpx_search_engine_scores",
        "title": f"Summary of {pretty}",
        "cpswitch": True,
        "cpswitch_c_active": False,
        "tt_decimals": 0,
        "ylab": "Count",
        "save_data_file": False,
    }

    bar_html = plot_html_check(bargraph.plot(data=plot_data, pconfig=draw_config))

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=1,
        description=f"Distribution of <code>{score_name}</code> per run.",
        helptext=f"""
            QPX stores engine scores as a free-form <code>additional_scores</code> list, so the
            score shown is whichever the project actually carries -- here
            <code>{score_name}</code>. Bin edges span the observed 1st-99th percentile
            range rather than a fixed scale, because score ranges differ by orders of
            magnitude between engines and between q-values, PEPs and raw scores.
            """,
    )
