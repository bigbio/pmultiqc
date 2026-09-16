from __future__ import annotations

import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table


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
        value = self.metadata.get("inputFiles", [])
        return value if isinstance(value, list) else []

    @property
    def input_file_names(self) -> list[str]:
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
        values: list[str] = []
        tokens = ("provenance", "evidence", "inferred", "unavailable")
        for metric in self.metrics:
            haystack = f"{metric.name} {metric.accession}".casefold()
            if any(token in haystack for token in tokens):
                value = _compact_value(metric.value)
                values.append(f"{metric.name}: {value}")
        return " | ".join(dict.fromkeys(values)) or None

    def metrics_by_category(self) -> dict[str, list[MzQCMetric]]:
        categories: dict[str, list[MzQCMetric]] = defaultdict(list)
        for metric in self.metrics:
            category = classify_metric(metric)
            categories[category].append(metric)
        return dict(categories)


def _metric_from_json(item: Any) -> MzQCMetric:
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
    if not isinstance(run_qualities, list):
        raise ValueError("mzQC document has no runQualities array")

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
    return json.dumps(value, separators=(",", ":"), ensure_ascii=False)


def _slug(value: str) -> str:
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
    return isinstance(value, list) and len(value) == 2 and all(
        isinstance(item, (int, float)) and not isinstance(item, bool) for item in value
    )


def _metric_key(metric: MzQCMetric, occurrence: int = 1) -> str:
    """Build a stable column key from a metric accession and occurrence."""
    base = _slug(metric.accession)
    return base if occurrence == 1 else f"{base}_{occurrence}"


def _metric_data(
    runs: list[MzQCRun],
    categories: set[str] | None = None,
    numeric_only: bool = False,
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    """Create table/general-stat data with the same columns across all runs."""
    data: dict[str, dict[str, Any]] = {}
    headers: dict[str, dict[str, Any]] = {}

    def matching(metric: MzQCMetric) -> bool:
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
            unit = metric.unit_name
            description = metric.description or ""
            if _is_two_number_tuple(metric.value):
                min_key = f"{key}_min"
                max_key = f"{key}_max"
                row[min_key] = metric.value[0]
                row[max_key] = metric.value[1]
                headers[min_key] = {
                    "title": f"{metric.name} min",
                    "description": description,
                }
                headers[max_key] = {
                    "title": f"{metric.name} max",
                    "description": description,
                }
                if unit:
                    headers[min_key]["suffix"] = f" {unit}"
                    headers[max_key]["suffix"] = f" {unit}"
            elif metric.scalar_value is not None:
                row[key] = metric.scalar_value
                headers[key] = {"title": metric.name, "description": description}
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

        matched_paths: set[Path] = set()
        for file_info in self.find_log_files("mzqc", filecontents=False):
            path = (Path(file_info["root"]) / file_info["fn"]).resolve()
            matched_paths.add(path)
            try:
                runs.extend(parse_mzqc_document(path, used_names))
                self.add_data_source(file_info)
            except (OSError, ValueError) as exc:
                errors.append(f"{path.name}: {exc}")

        # pmultiqc 0.0.48 is a plugin for MultiQC 1.35. In this compatibility
        # path, third-party search patterns can be registered after MultiQC has
        # already built its file index. When that happens find_log_files() is
        # empty even though the user explicitly supplied .mzQC inputs. Fall
        # back to the analysis paths already selected by MultiQC so that the
        # mzQC module remains usable without requiring a separate conversion or
        # configuration file. Once the core index supplies matches, this path
        # is not used and normal MultiQC filtering / data-source handling wins.
        if not matched_paths:
            analysis_dirs = getattr(config, "analysis_dir", []) or []
            fallback_paths: set[Path] = set()
            for analysis_dir in analysis_dirs:
                candidate = Path(analysis_dir).expanduser()
                if candidate.is_file() and candidate.suffix == ".mzQC":
                    fallback_paths.add(candidate.resolve())
                elif candidate.is_dir():
                    fallback_paths.update(
                        path.resolve()
                        for path in candidate.rglob("*.mzQC")
                        if path.is_file()
                    )
            for path in sorted(fallback_paths):
                try:
                    runs.extend(parse_mzqc_document(path, used_names))
                except (OSError, ValueError) as exc:
                    errors.append(f"{path.name}: {exc}")

        if errors:
            for message in errors:
                self.log.warning("Skipping invalid mzQC input: %s", message)

        self.runs = runs
        self._draw_report()

    def _draw_report(self) -> None:
        self.add_software_version(None)
        if not self.runs:
            raise ModuleNoSamplesFound

        run_data = self._run_overview_data()
        run_headers = {
            "source_file": {"title": "mzQC file"},
            "input_file": {"title": "Input file"},
            "instrument": {"title": "Instrument"},
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
        self._add_category_table("Mass Accuracy / Precision", "mzqc-mass", {"mass"})
        self._add_category_table("Other Scalar QC Metrics", "mzqc-other", {"other"})

        general_data, headers = _metric_data(self.runs, numeric_only=True)
        if general_data:
            visible_keys = self._preferred_general_stat_keys(headers)
            for key, header in headers.items():
                header["hidden"] = key not in visible_keys
            self.general_stats_addcols(general_data, headers)

        self.write_data_file(_all_metric_data(self.runs), "multiqc_mzqc")

    def _run_overview_data(self) -> dict[str, dict[str, Any]]:
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

    @staticmethod
    def _preferred_general_stat_keys(headers: dict[str, dict[str, Any]]) -> set[str]:
        preferred: set[str] = set()
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
                preferred.add(key)
        if len(preferred) < 2:
            preferred.update(list(headers)[:2])
        return set(list(preferred)[:8])
