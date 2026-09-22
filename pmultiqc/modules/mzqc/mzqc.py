from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from math import ceil, comb, log10
from pathlib import Path
from statistics import median
from typing import Any

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, heatmap, scatter, table


@dataclass(frozen=True, slots=True)
class MzQCMetric:
    """One mzQC qualityMetric with its original metadata preserved."""

    accession: str
    name: str
    value: Any
    unit_accession: str | None = None
    unit_name: str | None = None
    description: str | None = None

    @property
    def scalar_value(self) -> int | float | None:
        """Return a numeric scalar, or ``None`` for non-scalar values."""
        if isinstance(self.value, bool):
            return None
        if isinstance(self.value, (int, float)):
            return self.value
        return None


@dataclass(frozen=True, slots=True)
class MzQCRun:
    """A single mzQC runQuality record, independent of its source document."""

    sample_name: str
    source_path: Path
    metadata: dict[str, Any]
    metrics: tuple[MzQCMetric, ...]

    @property
    def input_files(self) -> list[dict[str, Any]]:
        """Return the run metadata's declared input-file objects."""
        value = self.metadata.get("inputFiles", [])
        return value if isinstance(value, list) else []

    @property
    def input_file_names(self) -> list[str]:
        """Return display names for all declared input files."""
        names: list[str] = []
        for item in self.input_files:
            if not isinstance(item, dict):
                continue
            name = item.get("name")
            location = item.get("location")
            if isinstance(name, str) and name:
                names.append(name)
            elif isinstance(location, str) and location:
                names.append(Path(location.split("?")[0]).name)
        return names

    @property
    def instrument(self) -> str | None:
        """Extract instrument labels from input-file properties when present."""
        values: list[str] = []
        for item in self.input_files:
            if not isinstance(item, dict):
                continue
            properties = item.get("fileProperties", [])
            if not isinstance(properties, list):
                continue
            for prop in properties:
                if not isinstance(prop, dict):
                    continue
                identity = " ".join(
                    str(prop.get(key, "")) for key in ("name", "accession")
                ).strip()
                if "instrument" in identity.casefold():
                    value = prop.get("value") or prop.get("name")
                    if value is not None:
                        values.append(str(value))
        return "; ".join(dict.fromkeys(values)) or None

    @property
    def analysis_software(self) -> str | None:
        """Return a compact label for the software declared by the run."""
        software = self.metadata.get("analysisSoftware", [])
        if not isinstance(software, list):
            return None
        labels: list[str] = []
        for item in software:
            if not isinstance(item, dict):
                continue
            name = item.get("name")
            version = item.get("version")
            if isinstance(name, str) and name:
                labels.append(f"{name} {version}" if version else name)
        return "; ".join(dict.fromkeys(labels)) or None

    @property
    def acquisition_method(self) -> str | None:
        """Return acquisition-method evidence represented by run metrics."""
        values: list[str] = []
        for metric in self.metrics:
            name = metric.name.casefold()
            if (
                name == "acquisition method"
                or "data-dependent acquisition" in name
                or "data-independent acquisition" in name
            ):
                values.append(_compact_value(metric.value))
        return " | ".join(dict.fromkeys(values)) or None

    @property
    def provenance(self) -> str | None:
        """Summarize provenance or evidence metrics without inventing metadata."""
        values: list[str] = []
        tokens = ("provenance", "evidence", "inferred", "unavailable")
        for metric in self.metrics:
            haystack = f"{metric.name} {metric.accession}".casefold()
            if any(token in haystack for token in tokens):
                value = _compact_value(metric.value)
                values.append(f"{metric.name}: {value}")
        return " | ".join(dict.fromkeys(values)) or None

    def metrics_by_category(self) -> dict[str, list[MzQCMetric]]:
        """Group the run's metrics by the reporting categories used below."""
        categories: dict[str, list[MzQCMetric]] = defaultdict(list)
        for metric in self.metrics:
            category = classify_metric(metric)
            categories[category].append(metric)
        return dict(categories)


def _metric_from_json(item: Any) -> MzQCMetric:
    """Validate and convert one JSON qualityMetric object."""
    if not isinstance(item, dict):
        raise ValueError("qualityMetrics entries must be JSON objects")
    accession = item.get("accession")
    name = item.get("name")
    if not isinstance(accession, str) or not accession:
        raise ValueError("qualityMetric is missing a string accession")
    if not isinstance(name, str) or not name:
        name = accession
    unit = item.get("unit")
    unit_accession = unit.get("accession") if isinstance(unit, dict) else None
    unit_name = unit.get("name") if isinstance(unit, dict) else None
    return MzQCMetric(
        accession=accession,
        name=name,
        value=item.get("value"),
        unit_accession=unit_accession if isinstance(unit_accession, str) else None,
        unit_name=unit_name if isinstance(unit_name, str) else None,
        description=item.get("description") if isinstance(item.get("description"), str) else None,
    )


def _run_label(run_quality: dict[str, Any], source_path: Path, index: int) -> str:
    """Derive a stable display label from run metadata and source identity."""
    metadata = run_quality.get("metadata", {})
    if isinstance(metadata, dict):
        label = metadata.get("label")
        if isinstance(label, str) and label.strip():
            return label.strip()
        input_files = metadata.get("inputFiles", [])
        if isinstance(input_files, list) and input_files:
            first = input_files[0]
            if isinstance(first, dict):
                name = first.get("name")
                if isinstance(name, str) and name.strip():
                    return Path(name).name
    return f"{source_path.stem} run {index}"


def _unique_name(base: str, used: set[str]) -> str:
    """Return a collision-free sample name while preserving the original base."""
    if base not in used:
        used.add(base)
        return base
    index = 2
    while f"{base} [{index}]" in used:
        index += 1
    result = f"{base} [{index}]"
    used.add(result)
    return result


def parse_mzqc_document(path: str | Path, used_names: set[str] | None = None) -> list[MzQCRun]:
    """Parse every ``runQuality`` object in one mzQC document.

    The parser uses the mzQC JSON representation directly. Optional ``pymzqc`` is
    deliberately not required, keeping the reporting path lightweight and usable
    wherever MultiQC itself is installed.
    """

    source_path = Path(path)
    with source_path.open("r", encoding="utf-8") as handle:
        document = json.load(handle)

    root = document.get("mzQC") if isinstance(document, dict) else None
    if not isinstance(root, dict):
        raise ValueError("missing mzQC root object")

    run_qualities = root.get("runQualities")
    if run_qualities is None:
        # Compatibility with pmultiqc's historical exporter and checked-in
        # reports, which used a non-standard singular key. New mzQC files
        # should use the specification-defined plural ``runQualities`` key.
        run_qualities = root.get("runQuality")
    if not isinstance(run_qualities, list):
        raise ValueError("mzQC document has no runQualities/runQuality array")

    names = used_names if used_names is not None else set()
    parsed: list[MzQCRun] = []
    for index, run_quality in enumerate(run_qualities, start=1):
        if not isinstance(run_quality, dict):
            raise ValueError(f"runQuality #{index} is not a JSON object")
        metadata = run_quality.get("metadata")
        if not isinstance(metadata, dict):
            metadata = {}
        quality_metrics = run_quality.get("qualityMetrics", [])
        if not isinstance(quality_metrics, list):
            raise ValueError(f"runQuality #{index} has invalid qualityMetrics")
        metrics = tuple(_metric_from_json(item) for item in quality_metrics)
        label = _unique_name(_run_label(run_quality, source_path, index), names)
        parsed.append(
            MzQCRun(
                sample_name=label,
                source_path=source_path,
                metadata=metadata,
                metrics=metrics,
            )
        )
    return parsed


def _compact_value(value: Any) -> str:
    """Serialize a metric value compactly for metadata-table display."""
    return json.dumps(value, separators=(",", ":"), ensure_ascii=False)


def _slug(value: str) -> str:
    """Convert a metric identity into a stable MultiQC column key."""
    value = re.sub(r"[^A-Za-z0-9]+", "_", value.strip()).strip("_")
    return value or "metric"


def classify_metric(metric: MzQCMetric) -> str:
    """Classify by the metric's declared name/accession without inventing values."""
    text = f"{metric.name} {metric.accession}".casefold()
    if any(
        token in text
        for token in ("mass error", "mass accuracy", "precision", "tolerance", "delta m")
    ):
        return "mass"
    if any(
        token in text
        for token in ("acquisition", "isolation", "duty cycle", "target", "charge", "polarity")
    ):
        return "acquisition"
    if any(token in text for token in ("chromat", "tic", "retention", "rt_", " rt", "xic")):
        return "chromatography"
    if any(
        token in text
        for token in ("ms1", "ms2", "spectrum", "spectra", "peak", "m/z", "mz range", "scan rate")
    ):
        return "spectra"
    return "other"


def _is_two_number_tuple(value: Any) -> bool:
    """Return whether a JSON value is a two-endpoint numeric range."""
    return isinstance(value, list) and len(value) == 2 and all(
        isinstance(item, (int, float)) and not isinstance(item, bool) for item in value
    )


def _metric_key(metric: MzQCMetric, occurrence: int = 1) -> str:
    """Build a stable column key from a metric accession and occurrence."""
    base = _slug(metric.accession)
    return base if occurrence == 1 else f"{base}_{occurrence}"


def _metric_display_title(metric: MzQCMetric, endpoint: str | None = None) -> str:
    """Return compact user-facing titles while preserving mzQC semantics."""
    name = metric.name.strip()
    mapping = {
        "IsolationWidth_MS2_Median": "MS2 isolation width — median",
        "IsolationWidth_MS2_Min": "MS2 isolation width — min",
        "IsolationWidth_MS2_Max": "MS2 isolation width — max",
        "IsolationWidth_MS2_Count": "MS2 isolation width — observations",
        "IsolationWidth_MS2_Fraction": "MS2 isolation width — fraction",
    }
    if name in mapping:
        return mapping[name]
    if name == "ScanWindow_MS1" and endpoint:
        return f"MS1 scan window {endpoint} limit"
    if name == "ObservedMzRange_MS1" and endpoint:
        return f"Observed MS1 m/z range — {endpoint}"
    if name == "ObservedMzRange_MS2" and endpoint:
        return f"Observed MS2 m/z range — {endpoint}"
    if name == "PrecursorMzRange_MS2" and endpoint:
        return f"MS2 precursor m/z range — {endpoint}"
    if name in {"RTRange_MS1", "RetentionTimeRange_MS1"} and endpoint:
        return f"MS1 retention-time range — {endpoint}"
    if name in {"RTRange_MS2", "RetentionTimeRange_MS2"} and endpoint:
        return f"MS2 retention-time range — {endpoint}"
    if endpoint:
        return f"{name} {endpoint}"
    return name


def _metric_unit_suffix(metric: MzQCMetric) -> str | None:
    """Return a display unit, restoring units for local prideQC aggregate metrics."""
    if metric.unit_name:
        return metric.unit_name
    name = metric.name
    if name.startswith("IsolationWidth_MS2_"):
        return "Th"
    if name.startswith(("ScanWindow_", "ObservedMzRange_", "PrecursorMzRange_")):
        return "m/z"
    if name.startswith(("RTRange_", "RetentionTimeRange_")):
        return "second"
    return None


def _metric_data(
    runs: list[MzQCRun],
    categories: set[str] | None = None,
    numeric_only: bool = False,
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    """Create table/general-stat data with the same columns across all runs."""
    data: dict[str, dict[str, Any]] = {}
    headers: dict[str, dict[str, Any]] = {}

    def matching(metric: MzQCMetric) -> bool:
        """Return whether a metric belongs in the requested category set."""
        return categories is None or classify_metric(metric) in categories

    # Determine a stable occurrence number within each run. mzQC accessions are
    # the primary semantic identity; an occurrence suffix handles repeated terms
    # without letting the ordering of other runs change the columns.
    for run in runs:
        occurrences: dict[str, int] = defaultdict(int)
        row: dict[str, Any] = {}
        for metric in run.metrics:
            if not matching(metric):
                continue
            occurrences[metric.accession] += 1
            occurrence = occurrences[metric.accession]
            key = _metric_key(metric, occurrence)
            unit = _metric_unit_suffix(metric)
            description = metric.description or ""
            if _is_two_number_tuple(metric.value):
                if numeric_only:
                    continue
                min_key = f"{key}_min"
                max_key = f"{key}_max"
                row[min_key] = metric.value[0]
                row[max_key] = metric.value[1]
                headers[min_key] = {
                    "title": _metric_display_title(metric, "lower"),
                    "description": description,
                }
                headers[max_key] = {
                    "title": _metric_display_title(metric, "upper"),
                    "description": description,
                }
                if unit:
                    headers[min_key]["suffix"] = f" {unit}"
                    headers[max_key]["suffix"] = f" {unit}"
            elif metric.scalar_value is not None:
                row[key] = metric.scalar_value
                headers[key] = {
                    "title": _metric_display_title(metric),
                    "description": description,
                }
                if unit:
                    headers[key]["suffix"] = f" {unit}"
        data[run.sample_name] = row
    return data, headers


def _all_metric_data(runs: list[MzQCRun]) -> dict[str, dict[str, Any]]:
    """Create the raw MultiQC data-file payload without discarding non-scalars."""
    result: dict[str, dict[str, Any]] = {}
    for run in runs:
        row: dict[str, Any] = {}
        counts: dict[str, int] = defaultdict(int)
        for metric in run.metrics:
            base = f"{metric.accession}::{metric.name}"
            counts[base] += 1
            key = base if counts[base] == 1 else f"{base} [{counts[base]}]"
            row[key] = metric.value
        result[run.sample_name] = row
    return result


MASS_ERROR_METRIC_NAMES = {
    "precursor_precision_ppm": "estimated precursor mass error precision",
    "precursor_precision_da": "estimated precursor mass error precision in daltons",
    "fragment_precision_ppm": "estimated fragment mass error precision",
    "fragment_precision_da": "estimated fragment mass error precision in daltons",
    "precursor_tolerance_ppm": "suggested precursor search tolerance",
    "fragment_tolerance_ppm": "suggested fragment search tolerance",
    "fragment_tolerance_da": "suggested fragment search tolerance in daltons",
    "diagnostics": "mass error estimator diagnostics",
}
MASS_SHIFT_METRIC_NAME = "putative modification mass shifts"
MASS_SHIFT_DIAGNOSTICS_NAME = "mass shift scout diagnostics"
MASS_SHIFT_CLASS_ORDER = (
    "putative-ptm",
    "sample-prep-modification",
    "putative-modification",
    "mass-compatible-other",
    "unknown",
    "isotope-like",
    "adduct-like",
)
MASS_SHIFT_CLASS_LABELS = {
    "putative-ptm": "PTM-compatible",
    "sample-prep-modification": "Sample-prep / artifact",
    "putative-modification": "Modification-compatible",
    "mass-compatible-other": "Other mass-compatible",
    "unknown": "Unknown",
    "isotope-like": "Isotope-like",
    "adduct-like": "Adduct-like",
}


def _metric_named(run: MzQCRun, name: str) -> MzQCMetric | None:
    """Return the first metric with the exact report-contract name."""
    wanted = name.casefold()
    return next((metric for metric in run.metrics if metric.name.casefold() == wanted), None)


def _structured_value(run: MzQCRun, name: str, expected: type) -> Any:
    """Return one structured metric only when it has the expected JSON shape."""
    metric = _metric_named(run, name)
    return metric.value if metric is not None and isinstance(metric.value, expected) else None


def _finite_number(value: Any) -> int | float | None:
    """Return an ordinary JSON number without accepting booleans."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return value


def _nested_number(value: dict[str, Any] | None, key: str) -> int | float | None:
    if value is None:
        return None
    return _finite_number(value.get(key))


def _mass_accuracy_report_data(
    runs: list[MzQCRun],
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    """Extract report-ready mass precision/tolerance values from prideQC mzQC annotations."""
    data: dict[str, dict[str, Any]] = {}
    ppm_plot: dict[str, dict[str, Any]] = {}
    for run in runs:
        precursor_precision = _structured_value(
            run, MASS_ERROR_METRIC_NAMES["precursor_precision_ppm"], dict
        )
        fragment_precision = _structured_value(
            run, MASS_ERROR_METRIC_NAMES["fragment_precision_ppm"], dict
        )
        precursor_tolerance = _structured_value(
            run, MASS_ERROR_METRIC_NAMES["precursor_tolerance_ppm"], dict
        )
        fragment_tolerance = _structured_value(
            run, MASS_ERROR_METRIC_NAMES["fragment_tolerance_ppm"], dict
        )
        fragment_tolerance_da = _structured_value(
            run, MASS_ERROR_METRIC_NAMES["fragment_tolerance_da"], dict
        )
        diagnostics = _structured_value(run, MASS_ERROR_METRIC_NAMES["diagnostics"], dict)

        row: dict[str, Any] = {}
        fields = {
            "precursor_sigma_ppm": _nested_number(precursor_precision, "single_measurement_sigma"),
            "precursor_tolerance_ppm": _nested_number(
                precursor_tolerance, "suggested_tolerance"
            ),
            "fragment_sigma_ppm": _nested_number(fragment_precision, "single_measurement_sigma"),
            "fragment_tolerance_ppm": _nested_number(fragment_tolerance, "suggested_tolerance"),
            "fragment_tolerance_da": _nested_number(fragment_tolerance_da, "suggested_tolerance"),
            "precursor_repeat_pairs": _nested_number(diagnostics, "precursor_paired_spectra"),
            "precursor_clusters": _nested_number(diagnostics, "precursor_clusters_used"),
            "fragment_pairs": _nested_number(diagnostics, "fragment_pairs"),
        }
        for key, value in fields.items():
            if value is not None:
                row[key] = value

        regime = None if diagnostics is None else diagnostics.get("fragment_resolution_regime")
        if isinstance(regime, str) and regime:
            row["fragment_resolution_regime"] = regime
        confidence = None if precursor_tolerance is None else precursor_tolerance.get("confidence")
        if isinstance(confidence, str) and confidence:
            row["precursor_tolerance_confidence"] = confidence
        fragment_payload = fragment_tolerance or fragment_tolerance_da
        fragment_confidence = (
            None if fragment_payload is None else fragment_payload.get("confidence")
        )
        if isinstance(fragment_confidence, str) and fragment_confidence:
            row["fragment_tolerance_confidence"] = fragment_confidence

        if row:
            data[run.sample_name] = row

        ppm_row = {
            key: value
            for key, value in (
                ("Precursor tolerance", fields["precursor_tolerance_ppm"]),
                ("Fragment tolerance", fields["fragment_tolerance_ppm"]),
            )
            if value is not None
        }
        if ppm_row:
            ppm_plot[run.sample_name] = ppm_row
    return data, ppm_plot


def _mass_accuracy_detail_description(data: dict[str, dict[str, Any]]) -> str:
    """Explain unit choice and make estimator abstention explicit in the report."""
    description = (
        "Run-level precision, suggested tolerances and estimator support from prideQC. "
        "High-resolution fragment recommendations are reported in ppm; low-resolution "
        "fragment recommendations are reported in Da. Unavailable or abstained estimates "
        "remain blank."
    )
    any_suggested_tolerance = any(
        any(
            key in row
            for key in (
                "precursor_tolerance_ppm",
                "fragment_tolerance_ppm",
                "fragment_tolerance_da",
            )
        )
        for row in data.values()
    )
    if not any_suggested_tolerance:
        description += (
            " No runs met the estimator criteria for a suggested search tolerance; "
            "precision diagnostics remain available below."
        )
    return description


def _mass_accuracy_summary_text(data: dict[str, dict[str, Any]], total_runs: int) -> str:
    """Create a compact cohort summary without inventing unavailable tolerances."""
    precursor = [
        float(row["precursor_tolerance_ppm"])
        for row in data.values()
        if "precursor_tolerance_ppm" in row
    ]
    fragment = [
        float(row["fragment_tolerance_ppm"])
        for row in data.values()
        if "fragment_tolerance_ppm" in row
    ]
    clauses = []
    if precursor:
        clauses.append(
            f"precursor ppm tolerance available for {len(precursor)}/{total_runs} runs "
            f"(median {median(precursor):.3g} ppm; range "
            f"{min(precursor):.3g}-{max(precursor):.3g})"
        )
    if fragment:
        clauses.append(
            f"high-resolution fragment ppm tolerance available for {len(fragment)}/{total_runs} "
            f"runs (median {median(fragment):.3g} ppm; range "
            f"{min(fragment):.3g}-{max(fragment):.3g})"
        )
    if not clauses:
        return "No supported ppm tolerance estimates were present."
    return "; ".join(clauses) + "."


def _metric_scalar_named(run: MzQCRun, name: str) -> float | None:
    metric = _metric_named(run, name)
    value = None if metric is None else metric.scalar_value
    return None if value is None else float(value)


def _run_group_assignments(
    runs: list[MzQCRun],
) -> tuple[dict[str, str], dict[str, dict[str, Any]], list[str]]:
    """Detect conservative potential acquisition/run groups for one report cohort."""
    if not runs:
        return {}, {}, []
    accuracy, _ = _mass_accuracy_report_data(runs)
    by_name = {run.sample_name: run for run in runs}
    features: dict[str, dict[str, Any]] = {}
    for run in runs:
        row = accuracy.get(run.sample_name, {})
        fragment_unit = "ppm" if "fragment_tolerance_ppm" in row else (
            "Da" if "fragment_tolerance_da" in row else ""
        )
        fragment_value = row.get("fragment_tolerance_ppm", row.get("fragment_tolerance_da"))
        features[run.sample_name] = {
            "precursor": row.get("precursor_tolerance_ppm"),
            "fragment": fragment_value,
            "fragment_unit": fragment_unit,
            "fragment_regime": row.get("fragment_resolution_regime", ""),
            "rt": _metric_scalar_named(run, "chromatography duration"),
            "isolation": _metric_scalar_named(run, "IsolationWidth_MS2_Median"),
            "instrument": run.instrument or "",
            "acquisition": run.acquisition_method or "",
        }

    names = [run.sample_name for run in runs]
    min_group = max(3, ceil(len(names) * 0.10))
    groups: list[list[str]] = [names]
    evidence: list[str] = []

    def supported_partition(group: list[str], key: str) -> list[list[str]] | None:
        values = {name: str(features[name].get(key) or "") for name in group}
        if any(not value for value in values.values()):
            return None
        buckets: dict[str, list[str]] = defaultdict(list)
        for name, value in values.items():
            buckets[value].append(name)
        if len(buckets) <= 1 or any(len(bucket) < min_group for bucket in buckets.values()):
            return None
        return [buckets[key] for key in sorted(buckets)]

    # Hard incompatibilities only split when every resulting bucket is well supported.
    for key in ("fragment_unit", "fragment_regime", "instrument", "acquisition"):
        changed = True
        while changed and len(groups) < 4:
            changed = False
            new_groups: list[list[str]] = []
            for group in groups:
                partition = supported_partition(group, key)
                if partition is not None and len(new_groups) + len(partition) <= 4:
                    new_groups.extend(partition)
                    evidence.append(f"categorical {key} split")
                    changed = True
                else:
                    new_groups.append(group)
            groups = new_groups

    thresholds = {
        "precursor": 1.75,
        "fragment": 1.75,
        "rt": 1.20,
        "isolation": 1.35,
    }

    def best_gap(group: list[str]) -> tuple[float, str, int, list[tuple[float, str]]] | None:
        best = None
        for key, threshold in thresholds.items():
            pairs = []
            for name in group:
                value = features[name].get(key)
                if isinstance(value, (int, float)) and value > 0:
                    pairs.append((float(value), name))
            if len(pairs) != len(group):
                continue
            pairs.sort()
            for index in range(min_group, len(pairs) - min_group + 1):
                low = pairs[index - 1][0]
                high = pairs[index][0]
                ratio = high / low if low > 0 else 1.0
                score = ratio / threshold
                if ratio >= threshold and (best is None or score > best[0]):
                    best = (score, key, index, pairs)
        return best

    while len(groups) < 4:
        choices = []
        for group_index, group in enumerate(groups):
            gap = best_gap(group)
            if gap is not None:
                choices.append((gap[0], group_index, gap))
        if not choices:
            break
        _, group_index, gap = max(choices, key=lambda item: item[0])
        _, key, index, pairs = gap
        left = [name for _, name in pairs[:index]]
        right = [name for _, name in pairs[index:]]
        low = pairs[index - 1][0]
        high = pairs[index][0]
        evidence.append(f"{key} gap {high / low:.2f}x")
        groups[group_index : group_index + 1] = [left, right]

    groups.sort(key=lambda group: min(names.index(name) for name in group))
    assignments: dict[str, str] = {}
    summaries: dict[str, dict[str, Any]] = {}
    for index, group in enumerate(groups, start=1):
        label = f"Group {index}"
        for name in group:
            assignments[name] = label
        rows = [accuracy.get(name, {}) for name in group]
        precursor = [float(row["precursor_tolerance_ppm"]) for row in rows if "precursor_tolerance_ppm" in row]
        fragment_ppm = [float(row["fragment_tolerance_ppm"]) for row in rows if "fragment_tolerance_ppm" in row]
        fragment_da = [float(row["fragment_tolerance_da"]) for row in rows if "fragment_tolerance_da" in row]
        rt = [float(features[name]["rt"]) for name in group if isinstance(features[name]["rt"], (int, float))]
        isolation = [float(features[name]["isolation"]) for name in group if isinstance(features[name]["isolation"], (int, float))]
        summary: dict[str, Any] = {"files": len(group)}
        if precursor:
            summary["precursor_median_ppm"] = median(precursor)
            summary["precursor_common_ppm"] = max(precursor)
        if fragment_ppm:
            summary["fragment_regime"] = "high-resolution / ppm"
            summary["fragment_median"] = median(fragment_ppm)
            summary["fragment_common"] = max(fragment_ppm)
            summary["fragment_unit"] = "ppm"
        elif fragment_da:
            summary["fragment_regime"] = "low-resolution / Da"
            summary["fragment_median"] = median(fragment_da)
            summary["fragment_common"] = max(fragment_da)
            summary["fragment_unit"] = "Da"
        if rt:
            summary["rt_median_min"] = median(rt) / 60.0
            summary["rt_min_min"] = min(rt) / 60.0
            summary["rt_max_min"] = max(rt) / 60.0
        if isolation:
            summary["isolation_median_th"] = median(isolation)
        instruments = sorted({by_name[name].instrument for name in group if by_name[name].instrument})
        acquisitions = sorted({by_name[name].acquisition_method for name in group if by_name[name].acquisition_method})
        if instruments:
            summary["instrument"] = "; ".join(instruments)
        if acquisitions:
            summary["acquisition"] = "; ".join(acquisitions)
        summaries[label] = summary
    return assignments, summaries, evidence


def _beta_prevalence_probability(hits: int, runs: int, threshold: float = 0.10) -> float:
    """Return P(prevalence > threshold | hits, runs, Beta(1,1))."""
    if runs <= 0 or hits < 0 or hits > runs:
        return 0.0
    # For Beta(hits+1, runs-hits+1), 1-I_x equals a binomial lower tail.
    total = runs + 1
    probability = sum(
        comb(total, index) * threshold**index * (1.0 - threshold) ** (total - index)
        for index in range(hits + 1)
    )
    return min(1.0, max(0.0, probability))


def _mass_shift_family_context(
    runs: list[MzQCRun], assignments: dict[str, str]
) -> tuple[dict[int, dict[str, Any]], list[dict[str, Any]]]:
    """Aggregate 0.02-Da recurrent mass families within inferred run groups."""
    group_members: dict[str, list[str]] = defaultdict(list)
    for run in runs:
        group_members[assignments.get(run.sample_name, "Group 1")].append(run.sample_name)

    record_context: dict[int, dict[str, Any]] = {}
    family_rows: list[dict[str, Any]] = []
    for group, members in group_members.items():
        observations: list[tuple[float, str, dict[str, Any]]] = []
        for run in runs:
            if run.sample_name not in members:
                continue
            for record in _mass_shift_records(run):
                delta = _finite_number(record.get("delta_mass_da"))
                if delta is not None:
                    observations.append((float(delta), run.sample_name, record))
        observations.sort(key=lambda item: item[0])
        families: list[list[tuple[float, str, dict[str, Any]]]] = []
        for observation in observations:
            if not families:
                families.append([observation])
                continue
            center = median([item[0] for item in families[-1]])
            if abs(observation[0] - center) <= 0.02:
                families[-1].append(observation)
            else:
                families.append([observation])
        for family in families:
            run_hits = sorted({item[1] for item in family})
            probability = _beta_prevalence_probability(len(run_hits), len(members))
            prevalence = len(run_hits) / len(members) if members else 0.0
            representative = max(
                family, key=lambda item: int(item[2].get("pair_support", 0) or 0)
            )
            rep_record = representative[2]
            candidates = _mass_shift_candidate_names(rep_record)
            high_support_runs = {
                run_name
                for _, run_name, record in family
                if str(record.get("confidence") or "") == "high-support"
            }
            context = {
                "run_group": group,
                "family_runs": len(run_hits),
                "group_runs": len(members),
                "run_prevalence": prevalence,
                "raw_prevalence_probability": probability,
                "mass_identity_ambiguous": len(candidates) > 1,
                "alternative_candidates": "; ".join(candidates[1:]),
            }
            for _, _, record in family:
                record_context[id(record)] = context
            if (
                str(rep_record.get("classification") or "") == "putative-ptm"
                and candidates
                and len(run_hits) >= 2
                and len(high_support_runs) >= 2
                and probability >= 0.95
            ):
                supports = [int(item[2].get("pair_support", 0) or 0) for item in family]
                family_rows.append(
                    {
                        "run_group": group,
                        "candidate": candidates[0],
                        "delta_mass_da": median([item[0] for item in family]),
                        "family_runs": len(run_hits),
                        "run_prevalence": prevalence,
                        "raw_prevalence_probability": probability,
                        "median_pair_support": median(supports) if supports else 0,
                        "mass_identity_ambiguous": len(candidates) > 1,
                        "alternative_candidates": "; ".join(candidates[1:]),
                    }
                )
    family_rows.sort(
        key=lambda row: (float(row["raw_prevalence_probability"]), int(row["family_runs"])),
        reverse=True,
    )
    return record_context, family_rows


def _mass_shift_records(run: MzQCRun) -> list[dict[str, Any]]:
    """Return the bounded QC-facing mass-shift records for one run."""
    value = _structured_value(run, MASS_SHIFT_METRIC_NAME, list)
    if value is None:
        return []
    return [item for item in value if isinstance(item, dict)]


def _mass_shift_candidate_names(record: dict[str, Any]) -> list[str]:
    """Return displayed mass-compatible candidate labels in prideQC order."""
    names: list[str] = []
    candidates = record.get("unimod_candidates")
    if isinstance(candidates, list):
        for candidate in candidates:
            if not isinstance(candidate, dict):
                continue
            name = candidate.get("name") or candidate.get("unimod_name")
            if isinstance(name, str) and name and name not in names:
                names.append(name)
    if names:
        return names
    artifact = record.get("artifact_candidate")
    if isinstance(artifact, dict):
        name = artifact.get("artifact_name") or artifact.get("name")
        if isinstance(name, str) and name:
            return [name]
    return []


def _mass_shift_report_data(runs: list[MzQCRun]) -> dict[str, Any]:
    """Build bounded tables and plots from prideQC's reported mass-shift evidence."""
    by_run: dict[str, list[dict[str, Any]]] = {}
    diagnostics: dict[str, dict[str, Any]] = {}
    classification_plot: dict[str, dict[str, int]] = {}
    scatter_data: dict[str, list[dict[str, float]]] = defaultdict(list)
    all_rows: list[tuple[str, dict[str, Any]]] = []
    class_counts: Counter[str] = Counter()
    confidence_counts: Counter[str] = Counter()

    family_total_support: Counter[float] = Counter()
    family_best: dict[float, tuple[int, dict[str, Any]]] = {}
    family_by_run: dict[str, Counter[float]] = defaultdict(Counter)

    for run in runs:
        records = _mass_shift_records(run)
        if records:
            by_run[run.sample_name] = records
        diag = _structured_value(run, MASS_SHIFT_DIAGNOSTICS_NAME, dict)
        if diag is not None:
            diagnostics[run.sample_name] = diag

        row_counts = Counter(
            str(record.get("classification") or "unknown") for record in records
        )
        if row_counts:
            classification_plot[run.sample_name] = {
                MASS_SHIFT_CLASS_LABELS[category]: row_counts[category]
                for category in MASS_SHIFT_CLASS_ORDER
                if row_counts[category]
            }

        for record in records:
            classification = str(record.get("classification") or "unknown")
            confidence = str(record.get("confidence") or "")
            class_counts[classification] += 1
            if confidence:
                confidence_counts[confidence] += 1
            all_rows.append((run.sample_name, record))

            delta = _finite_number(record.get("delta_mass_da"))
            support = _finite_number(record.get("pair_support"))
            if delta is None or support is None:
                continue
            scatter_data[classification].append({"x": float(delta), "y": float(support)})

            family = round(float(delta), 2)
            support_i = int(support)
            family_total_support[family] += support_i
            family_by_run[run.sample_name][family] += support_i
            previous = family_best.get(family)
            if previous is None or support_i > previous[0]:
                family_best[family] = (support_i, record)

    top_families = [family for family, _ in family_total_support.most_common(12)]
    family_labels: dict[float, str] = {}
    for family in top_families:
        record = family_best[family][1]
        representative = float(record.get("delta_mass_da", family))
        family_labels[family] = f"{representative:+.3f} Da"

    heatmap_data: dict[str, dict[str, float]] = {}
    if top_families:
        for run in runs:
            counts = family_by_run.get(run.sample_name, Counter())
            heatmap_data[run.sample_name] = {
                family_labels[family]: log10(1.0 + counts[family])
                for family in top_families
            }

    assignments, _, _ = _run_group_assignments(runs)
    family_context, high_support_families = _mass_shift_family_context(runs, assignments)
    return {
        "by_run": by_run,
        "diagnostics": diagnostics,
        "classification_plot": classification_plot,
        "scatter_data": dict(scatter_data),
        "all_rows": all_rows,
        "class_counts": class_counts,
        "confidence_counts": confidence_counts,
        "heatmap_data": heatmap_data,
        "run_groups": assignments,
        "family_context": family_context,
        "high_support_families": high_support_families,
    }


def _mass_shift_summary_description(report: dict[str, Any]) -> str:
    """Describe the scout conservatively and call out evidence-free cohorts."""
    description = (
        "Identification-free recurrent neutral precursor-mass differences from related "
        "MS2 spectra. Candidate annotations are mass-compatible interpretations, not "
        "peptide- or site-localized PTM identifications."
    )
    if not report["by_run"]:
        description += (
            " No recurrent modification-compatible mass-shift clusters were detected in "
            "the analysed runs."
        )
    return description


def _mass_shift_summary_table(report: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Create one report-level summary row for the bounded mass-shift scout output."""
    diagnostics = report["diagnostics"]
    raw_clusters = sum(
        int(value.get("raw_recurrent_clusters", 0) or 0) for value in diagnostics.values()
    )
    profile_ms2 = sum(
        int(value.get("explicit_profile_ms2", 0) or 0)
        for value in diagnostics.values()
    )
    profile_failures = sum(
        int(value.get("profile_peak_pick_failures", 0) or 0) for value in diagnostics.values()
    )
    classes: Counter[str] = report["class_counts"]
    confidence: Counter[str] = report["confidence_counts"]
    reported = sum(classes.values())
    return {
        "Cohort": {
            "runs_analyzed": len(diagnostics) or len(report["by_run"]),
            "runs_with_mass_shifts": len(report["by_run"]),
            "reported_clusters": reported,
            "raw_recurrent_clusters": raw_clusters,
            "high_support": confidence["high-support"],
            "putative_ptm": classes["putative-ptm"],
            "sample_prep": classes["sample-prep-modification"],
            "artifacts": classes["isotope-like"] + classes["adduct-like"],
            "unknown": classes["unknown"],
            "profile_ms2": profile_ms2,
            "profile_peak_pick_failures": profile_failures,
        }
    }


def _mass_shift_cluster_table(
    report: dict[str, Any], limit: int = 250
) -> dict[str, dict[str, Any]]:
    """Return the strongest reported clusters for an interactive report table."""
    rows = sorted(
        report["all_rows"],
        key=lambda item: int(item[1].get("pair_support", 0) or 0),
        reverse=True,
    )[:limit]
    result: dict[str, dict[str, Any]] = {}
    for index, (run_name, record) in enumerate(rows, start=1):
        delta = _finite_number(record.get("delta_mass_da"))
        candidate_names = _mass_shift_candidate_names(record)
        candidate = candidate_names[0] if candidate_names else ""
        context = report.get("family_context", {}).get(id(record), {})
        result[f"{index:03d} · {run_name}"] = {
            "run": run_name,
            "run_group": context.get("run_group", "Group 1"),
            "delta_mass_da": delta,
            "classification": MASS_SHIFT_CLASS_LABELS.get(
                str(record.get("classification") or "unknown"),
                str(record.get("classification") or "Unknown"),
            ),
            "candidate": candidate,
            "alternative_candidates": "; ".join(candidate_names[1:]),
            "pair_support": int(record.get("pair_support", 0) or 0),
            "unique_spectra": int(record.get("unique_spectrum_support", 0) or 0),
            "spectral_similarity": _finite_number(record.get("median_spectral_similarity")),
            "confidence": str(record.get("confidence") or ""),
            "family_runs": context.get("family_runs"),
            "run_prevalence": context.get("run_prevalence"),
            "raw_prevalence_probability": context.get("raw_prevalence_probability"),
            "mass_identity_ambiguous": context.get("mass_identity_ambiguous"),
            "mass_residual_da": (
                _finite_number(
                    record["unimod_candidates"][0].get(
                        "residual_da",
                        record["unimod_candidates"][0].get("unimod_residual_da"),
                    )
                )
                if candidate
                and isinstance(record.get("unimod_candidates"), list)
                and record["unimod_candidates"]
                and isinstance(record["unimod_candidates"][0], dict)
                else None
            ),
        }
    return result


class MzQCModule(BaseMultiqcModule):
    """
    MultiQC module for mzQC 1.0 run-level quality-control documents.

    The module reads ``.mzQC`` JSON directly. Each ``runQuality`` object is treated
    as one MultiQC sample, so a directory of documents and a single multi-run document
    both aggregate naturally into one dataset report. Aggregate mzQC metrics are never
    expanded into synthetic per-spectrum observations.

    The module reports only values actually present in the input mzQC documents. Scalar
    metrics feed the General Statistics table and category-specific tables; two-value
    range metrics are displayed as their declared endpoints. Non-scalar metrics remain
    available in the ``multiqc_mzqc`` data file rather than being silently discarded.
    """

    def __init__(self) -> None:
        """Discover indexed mzQC files, parse their runs, and render the module."""
        super().__init__(
            name="mzQC",
            target="mzQC",
            anchor="mzqc",
            href="https://github.com/HUPO-PSI/mzQC",
            info=(
                "A standards-based format for reporting and exchanging mass spectrometry "
                "quality-control metrics."
            ),
        )

        used_names: set[str] = set()
        runs: list[MzQCRun] = []
        errors: list[str] = []

        for file_info in self.find_log_files("mzqc", filecontents=False):
            path = (Path(file_info["root"]) / file_info["fn"]).resolve()
            try:
                runs.extend(parse_mzqc_document(path, used_names))
                self.add_data_source(file_info)
            except (OSError, ValueError) as exc:
                errors.append(f"{path.name}: {exc}")

        if errors:
            for message in errors:
                self.log.warning("Skipping invalid mzQC input: %s", message)

        self.runs = runs
        self._draw_report()

    def _draw_report(self) -> None:
        """Render overview, category, general-statistics, and raw-data outputs."""
        self.add_software_version(None)
        if not self.runs:
            raise ModuleNoSamplesFound

        self._run_groups, self._run_group_summaries, self._run_group_evidence = (
            _run_group_assignments(self.runs)
        )
        if len(self._run_group_summaries) > 1:
            details = "; ".join(
                f"{group}: {summary['files']} files"
                for group, summary in self._run_group_summaries.items()
            )
            evidence = ", ".join(self._run_group_evidence)
            self.add_section(
                name="Potential Acquisition Heterogeneity",
                anchor="mzqc-potential-acquisition-heterogeneity",
                description=(
                    f"Detected {len(self._run_group_summaries)} potential acquisition/run "
                    f"groups ({details}). This is a conservative QC warning, not a claim of "
                    f"distinct biological experiments. Separation evidence: {evidence}."
                ),
            )

        run_data = self._run_overview_data()
        run_headers = {
            "source_file": {"title": "mzQC file"},
            "input_file": {"title": "Input file"},
            "instrument": {"title": "Instrument"},
            "acquisition_method": {"title": "Acquisition method"},
            "software": {"title": "Analysis software"},
            "provenance": {"title": "Provenance / evidence"},
        }
        self.add_section(
            name="Run Overview",
            anchor="mzqc-run-overview",
            description=(
                "One row per mzQC runQuality object, preserving source and metadata identity."
            ),
            plot=table.plot(
                run_data,
                headers=run_headers,
                pconfig={
                    "id": "mzqc_run_overview_table",
                    "title": "Run Overview",
                    "save_file": False,
                    "raw_data_fn": "mzqc_run_overview",
                    "only_defined_headers": True,
                    "sort_rows": False,
                    "no_violin": True,
                    "save_data_file": False,
                },
            ),
        )

        self._add_category_table("Chromatography", "mzqc-chromatography", {"chromatography"})
        self._add_category_table("MS1 / MS2 Summary", "mzqc-spectra", {"spectra"})
        self._add_category_table("Acquisition", "mzqc-acquisition", {"acquisition"})
        self._add_mass_accuracy_sections()
        self._add_mass_shift_sections()
        self._add_category_table(
            "Other Mass Accuracy / Precision Metrics",
            "mzqc-mass-other",
            {"mass"},
        )
        self._add_category_table("Other Scalar QC Metrics", "mzqc-other", {"other"})

        general_data, headers = _metric_data(self.runs, numeric_only=True)
        if general_data:
            visible_keys = self._preferred_general_stat_keys(headers)
            for key, header in headers.items():
                header["hidden"] = key not in visible_keys
            self.general_stats_addcols(general_data, headers)

        self._add_mass_accuracy_general_stats()
        self.write_data_file(_all_metric_data(self.runs), "multiqc_mzqc")

    def _run_overview_data(self) -> dict[str, dict[str, Any]]:
        """Build one metadata overview row per parsed mzQC run."""
        result: dict[str, dict[str, Any]] = {}
        for run in self.runs:
            result[run.sample_name] = {
                "source_file": run.source_path.name,
                "input_file": ", ".join(run.input_file_names),
                "instrument": run.instrument or "",
                "acquisition_method": run.acquisition_method or "",
                "software": run.analysis_software or "",
                "provenance": run.provenance or "",
            }
        return result

    def _add_category_table(self, title: str, anchor: str, categories: set[str]) -> None:
        """Add a metric table when at least one run has data for the category."""
        data, headers = _metric_data(self.runs, categories)
        if not any(row for row in data.values()):
            return
        self.add_section(
            name=title,
            anchor=anchor,
            description=(
                "Scalar and two-endpoint values directly represented in the input mzQC documents."
            ),
            plot=table.plot(
                data,
                headers=headers,
                pconfig={
                    "id": f"{anchor}_table",
                    "title": title,
                    "save_file": False,
                    "raw_data_fn": _slug(anchor),
                    "sort_rows": False,
                    "only_defined_headers": True,
                    "no_violin": True,
                    "save_data_file": False,
                },
            ),
        )

    def _add_mass_accuracy_sections(self) -> None:
        """Render prideQC's frozen mass-precision/tolerance evidence when present."""
        data, ppm_plot = _mass_accuracy_report_data(self.runs)
        if not data:
            return

        if self._run_group_summaries:
            summary_headers = {
                "files": {"title": "Files", "format": "{:,.0f}"},
                "precursor_median_ppm": {"title": "Precursor median", "suffix": " ppm", "format": "{:,.2f}"},
                "precursor_common_ppm": {"title": "Precursor common/max", "suffix": " ppm", "format": "{:,.2f}"},
                "fragment_median": {"title": "Fragment median", "format": "{:,.3f}"},
                "fragment_common": {"title": "Fragment common/max", "format": "{:,.3f}"},
                "fragment_unit": {"title": "Fragment unit"},
                "fragment_regime": {"title": "Fragment regime"},
                "rt_median_min": {"title": "RT median", "suffix": " min", "format": "{:,.1f}"},
                "rt_min_min": {"title": "RT min", "suffix": " min", "format": "{:,.1f}"},
                "rt_max_min": {"title": "RT max", "suffix": " min", "format": "{:,.1f}"},
                "isolation_median_th": {"title": "MS2 isolation width — median", "suffix": " Th", "format": "{:,.2f}"},
                "instrument": {"title": "Instrument"},
                "acquisition": {"title": "Acquisition method"},
            }
            self.add_section(
                name="Reanalysis Tolerance Summary",
                anchor="mzqc-reanalysis-tolerance-summary",
                description=(
                    "One conservative common search-tolerance recommendation per potential "
                    "acquisition/run group. Common values use the maximum supported run-level "
                    "recommendation so one setting covers every run in that group; medians are "
                    "shown for context."
                ),
                plot=table.plot(
                    self._run_group_summaries,
                    headers=summary_headers,
                    pconfig={
                        "id": "mzqc_reanalysis_tolerance_summary",
                        "title": "Reanalysis Tolerance Summary",
                        "only_defined_headers": True,
                        "sort_rows": False,
                        "no_violin": True,
                        "save_data_file": False,
                    },
                ),
            )

        headers = {
            "precursor_sigma_ppm": {
                "title": "Precursor precision σ",
                "description": "Robust single-measurement precursor precision estimate.",
                "suffix": " ppm",
                "format": "{:,.3f}",
            },
            "precursor_tolerance_ppm": {
                "title": "Precursor search tolerance",
                "description": (
                    "Precision-derived suggested precursor search-tolerance starting point."
                ),
                "suffix": " ppm",
                "format": "{:,.3f}",
            },
            "fragment_sigma_ppm": {
                "title": "Fragment precision σ",
                "description": "Robust single-measurement fragment precision estimate.",
                "suffix": " ppm",
                "format": "{:,.3f}",
            },
            "fragment_tolerance_ppm": {
                "title": "Fragment tolerance (ppm)",
                "description": (
                    "Suggested fragment search tolerance for high-resolution fragment data."
                ),
                "suffix": " ppm",
                "format": "{:,.3f}",
            },
            "fragment_tolerance_da": {
                "title": "Fragment tolerance (Da)",
                "description": (
                    "Suggested fragment search tolerance for low-resolution fragment data."
                ),
                "suffix": " Da",
                "format": "{:,.4f}",
            },
            "fragment_resolution_regime": {"title": "Fragment regime"},
            "precursor_tolerance_confidence": {"title": "Precursor confidence"},
            "fragment_tolerance_confidence": {"title": "Fragment confidence"},
            "precursor_repeat_pairs": {
                "title": "Precursor repeat pairs",
                "format": "{:,.0f}",
            },
            "precursor_clusters": {"title": "Precursor clusters", "format": "{:,.0f}"},
            "fragment_pairs": {"title": "Fragment pairs", "format": "{:,.0f}"},
        }

        if ppm_plot:
            self.add_section(
                name="Estimated Search Tolerances",
                anchor="mzqc-estimated-search-tolerances",
                description=(
                    _mass_accuracy_summary_text(data, len(self.runs))
                    + " Values are measurement-precision-derived starting points, not recovered "
                    "historical database-search settings."
                ),
                plot=bargraph.plot(
                    data=ppm_plot,
                    pconfig={
                        "id": "mzqc_estimated_search_tolerances_ppm",
                        "title": "Estimated Search Tolerances",
                        "ylab": "Tolerance (ppm)",
                        "tt_decimals": 3,
                        "cpswitch": False,
                        "stacking": "group",
                        "sort_samples": False,
                        "save_data_file": False,
                    },
                ),
            )

        self.add_section(
            name="Mass Accuracy & Tolerance Details",
            anchor="mzqc-mass-accuracy-tolerance-details",
            description=_mass_accuracy_detail_description(data),
            plot=table.plot(
                data,
                headers=headers,
                pconfig={
                    "id": "mzqc_mass_accuracy_tolerance_table",
                    "title": "Mass Accuracy & Tolerance Details",
                    "save_file": False,
                    "raw_data_fn": "mzqc_mass_accuracy_tolerance",
                    "sort_rows": False,
                    "only_defined_headers": True,
                    "no_violin": True,
                    "save_data_file": False,
                },
            ),
        )

    def _add_mass_accuracy_general_stats(self) -> None:
        """Expose common group tolerances instead of different per-run recommendations."""
        data, _ = _mass_accuracy_report_data(self.runs)
        assignments = getattr(self, "_run_groups", {})
        summaries = getattr(self, "_run_group_summaries", {})
        general_data: dict[str, dict[str, Any]] = {}
        for sample, row in data.items():
            group = assignments.get(sample, "Group 1")
            summary = summaries.get(group, {})
            values: dict[str, Any] = {}
            if "precursor_common_ppm" in summary:
                values["prideqc_precursor_tolerance_ppm"] = summary["precursor_common_ppm"]
            if summary.get("fragment_unit") == "ppm" and "fragment_common" in summary:
                values["prideqc_fragment_tolerance_ppm"] = summary["fragment_common"]
            if len(summaries) > 1:
                values["prideqc_run_group"] = group
            if values:
                general_data[sample] = values
        if not general_data:
            return
        headers: dict[str, dict[str, Any]] = {
            "prideqc_precursor_tolerance_ppm": {
                "title": "Precursor tol.",
                "description": "Common/max prideQC precursor tolerance for this run group",
                "suffix": " ppm",
                "format": "{:,.2f}",
                "hidden": False,
            },
            "prideqc_fragment_tolerance_ppm": {
                "title": "Fragment tol.",
                "description": "Common/max high-resolution fragment tolerance for this run group",
                "suffix": " ppm",
                "format": "{:,.2f}",
                "hidden": False,
            },
        }
        if len(summaries) > 1:
            headers["prideqc_run_group"] = {
                "title": "Run group",
                "description": "Potential acquisition/run group detected from mzQC evidence",
                "hidden": False,
            }
        self.general_stats_addcols(general_data, headers)

    def _add_mass_shift_sections(self) -> None:
        """Render bounded recurrent mass-shift evidence without implying PTM identification."""
        report = _mass_shift_report_data(self.runs)
        if not report["by_run"] and not report["diagnostics"]:
            return

        summary_data = _mass_shift_summary_table(report)
        summary_headers = {
            "runs_analyzed": {"title": "Runs analysed", "format": "{:,.0f}"},
            "runs_with_mass_shifts": {"title": "Runs with shifts", "format": "{:,.0f}"},
            "reported_clusters": {"title": "Reported clusters", "format": "{:,.0f}"},
            "raw_recurrent_clusters": {
                "title": "Raw recurrent clusters",
                "format": "{:,.0f}",
            },
            "high_support": {"title": "High-support", "format": "{:,.0f}"},
            "putative_ptm": {"title": "PTM-compatible", "format": "{:,.0f}"},
            "sample_prep": {"title": "Sample-prep / artifact", "format": "{:,.0f}"},
            "artifacts": {"title": "Isotope / adduct", "format": "{:,.0f}"},
            "unknown": {"title": "Unknown", "format": "{:,.0f}"},
            "profile_ms2": {"title": "Profile MS2", "format": "{:,.0f}"},
            "profile_peak_pick_failures": {
                "title": "Peak-pick failures",
                "format": "{:,.0f}",
            },
        }
        self.add_section(
            name="Putative Modification Mass Shifts",
            anchor="mzqc-putative-modification-mass-shifts",
            description=_mass_shift_summary_description(report),
            plot=table.plot(
                summary_data,
                headers=summary_headers,
                pconfig={
                    "id": "mzqc_mass_shift_summary",
                    "title": "Mass-shift Scout Summary",
                    "only_defined_headers": True,
                    "sort_rows": False,
                    "no_violin": True,
                    "save_data_file": False,
                },
            ),
        )

        high_support_rows = {
            f"{index:03d}": row
            for index, row in enumerate(report.get("high_support_families", []), start=1)
        }
        if high_support_rows:
            self.add_section(
                name="High-support Modification Families",
                anchor="mzqc-high-support-modification-families",
                description=(
                    "Recurrent PTM-compatible mass-shift families with at least two supporting "
                    "runs, at least two high-support runs, and P(prevalence > 10%) >= 0.95. "
                    "These are mass-compatible families, not peptide- or site-localized PTM "
                    "identifications."
                ),
                plot=table.plot(
                    high_support_rows,
                    headers={
                        "run_group": {"title": "Run group"},
                        "candidate": {"title": "Mass-compatible candidate"},
                        "delta_mass_da": {"title": "Median Δ mass", "suffix": " Da", "format": "{:,.5f}"},
                        "family_runs": {"title": "Supporting runs", "format": "{:,.0f}"},
                        "run_prevalence": {"title": "Run prevalence", "format": "{:,.1%}"},
                        "raw_prevalence_probability": {"title": "P(prevalence > 10%)", "format": "{:,.3f}"},
                        "median_pair_support": {"title": "Median pair support", "format": "{:,.0f}"},
                        "mass_identity_ambiguous": {"title": "Mass ambiguous"},
                        "alternative_candidates": {"title": "Alternative candidates"},
                    },
                    pconfig={
                        "id": "mzqc_high_support_modification_families",
                        "title": "High-support Modification Families",
                        "only_defined_headers": True,
                        "sort_rows": False,
                        "no_violin": True,
                        "save_data_file": False,
                    },
                ),
            )

        if report["classification_plot"]:
            self.add_section(
                name="Mass-shift Classification by Run",
                anchor="mzqc-mass-shift-classification",
                description=(
                    "QC-facing reported clusters grouped by prideQC classification. Raw and "
                    "suppressed clusters remain available in the mzQC diagnostics."
                ),
                plot=bargraph.plot(
                    data=report["classification_plot"],
                    pconfig={
                        "id": "mzqc_mass_shift_classification",
                        "title": "Mass-shift Classification by Run",
                        "ylab": "Reported clusters",
                        "tt_decimals": 0,
                        "cpswitch": False,
                        "sort_samples": False,
                        "save_data_file": False,
                    },
                ),
            )

        if report["heatmap_data"]:
            self.add_section(
                name="Recurrent Mass-shift Families",
                anchor="mzqc-mass-shift-family-heatmap",
                description=(
                    "Top 12 report-wide mass-shift families grouped into 0.01-Da display bins. "
                    "Columns are labelled only by observed delta mass; candidate chemistry "
                    "remains available in the candidate table. Cell intensity is log10(1 + pair "
                    "support), so lower-support recurrent families remain visible alongside "
                    "dominant chemistry."
                ),
                plot=heatmap.plot(
                    data=report["heatmap_data"],
                    pconfig={
                        "id": "mzqc_mass_shift_family_heatmap",
                        "title": "Recurrent Mass-shift Families",
                        "xlab": "Mass-shift family",
                        "ylab": "Run",
                        "zlab": "log10(1 + pair support)",
                        "tt_decimals": 3,
                        "display_values": False,
                        "xcats_samples": False,
                        "square": False,
                        "cluster_rows": False,
                        "cluster_cols": False,
                        "save_data_file": False,
                    },
                ),
            )

        if report["scatter_data"]:
            self.add_section(
                name="Mass-shift Landscape",
                anchor="mzqc-mass-shift-landscape",
                description=(
                    "Each point is one bounded reported recurrent cluster. The x-axis is absolute "
                    "neutral delta mass and the y-axis is related-spectrum pair support on a "
                    "logarithmic scale. Classification is available from the companion bar chart."
                ),
                plot=scatter.plot(
                    data=report["scatter_data"],
                    pconfig={
                        "id": "mzqc_mass_shift_landscape",
                        "title": "Mass-shift Landscape",
                        "xlab": "Delta mass (Da)",
                        "ylab": "Pair support (log scale)",
                        "ylog": True,
                        "showlegend": False,
                        "save_data_file": False,
                    },
                ),
            )

        cluster_table = _mass_shift_cluster_table(report)
        if cluster_table:
            self.add_section(
                name="Top Mass-shift Candidates",
                anchor="mzqc-top-mass-shift-candidates",
                description=(
                    "Strongest 250 reported clusters by pair support. The complete bounded "
                    "structured evidence remains available in the mzQC data export."
                ),
                plot=table.plot(
                    cluster_table,
                    headers={
                        "run": {"title": "Run"},
                        "run_group": {"title": "Run group"},
                        "delta_mass_da": {
                            "title": "Δ mass",
                            "suffix": " Da",
                            "format": "{:,.5f}",
                        },
                        "classification": {"title": "Classification"},
                        "candidate": {"title": "Mass-compatible candidate"},
                        "alternative_candidates": {"title": "Alternative candidates"},
                        "pair_support": {"title": "Pair support", "format": "{:,.0f}"},
                        "unique_spectra": {"title": "Unique spectra", "format": "{:,.0f}"},
                        "spectral_similarity": {
                            "title": "Median spectral similarity",
                            "format": "{:,.3f}",
                        },
                        "confidence": {"title": "Support tier"},
                        "family_runs": {"title": "Family runs", "format": "{:,.0f}"},
                        "run_prevalence": {"title": "Run prevalence", "format": "{:,.1%}"},
                        "raw_prevalence_probability": {
                            "title": "P(prevalence > 10%)",
                            "description": "RAW family recurrence probability; not PTM identity probability.",
                            "format": "{:,.3f}",
                        },
                        "mass_identity_ambiguous": {"title": "Mass ambiguous"},
                        "mass_residual_da": {
                            "title": "Candidate residual",
                            "suffix": " Da",
                            "format": "{:,.5f}",
                        },
                    },
                    pconfig={
                        "id": "mzqc_top_mass_shift_candidates",
                        "title": "Top Mass-shift Candidates",
                        "only_defined_headers": True,
                        "sort_rows": False,
                        "no_violin": True,
                        "save_data_file": False,
                    },
                ),
            )

    @staticmethod
    def _preferred_general_stat_keys(headers: dict[str, dict[str, Any]]) -> list[str]:
        """Choose up to eight preferred General Statistics columns deterministically."""
        preferred: list[str] = []
        patterns = (
            "number of ms1 spectra",
            "number of ms2 spectra",
            "number of spectra",
            "chromatography duration",
            "scan rate",
            "acquisition cycle count",
            "mass error",
            "precision",
            "tolerance",
        )
        for key, header in headers.items():
            title = str(header.get("title", "")).casefold()
            if any(pattern in title for pattern in patterns):
                preferred.append(key)
        if len(preferred) < 2:
            for key in headers:
                if key not in preferred:
                    preferred.append(key)
                if len(preferred) >= 2:
                    break
        return preferred[:8]
