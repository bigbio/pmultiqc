"""Tests for direct mzQC -> MultiQC input handling."""

from __future__ import annotations

import json
import os
from html import unescape
import shutil
import subprocess
from pathlib import Path

import pytest

from pmultiqc.modules.mzqc import MzQCMetric, MzQCModule, MzQCRun, parse_mzqc_document
from pmultiqc.modules.mzqc.mzqc import (
    _mass_accuracy_report_data,
    _mass_shift_cluster_table,
    _mass_shift_report_data,
    _mass_shift_summary_table,
    _metric_data,
)


FIXTURE_DIR = Path("tests/resources/mzqc")


def test_parse_single_file_returns_one_run():
    """Parse one standards-shaped document and expose its run metadata."""
    runs = parse_mzqc_document(FIXTURE_DIR / "run1.mzQC")
    assert len(runs) == 1
    assert isinstance(runs[0], MzQCRun)
    assert runs[0].sample_name == "synthetic_run_1"
    assert len(runs[0].metrics) == 7
    assert runs[0].instrument == "Orbitrap Eclipse"
    assert runs[0].acquisition_method == "\"Data-dependent acquisition\""
    assert runs[0].provenance == "acquisition method provenance: \"inferred\""


def test_parse_single_document_with_multiple_run_qualities():
    """Preserve every runQuality object in a multi-run mzQC document."""
    used = set()
    runs = parse_mzqc_document(FIXTURE_DIR / "multi_run.mzQC", used)

    assert [run.sample_name for run in runs] == [
        "embedded_run_1",
        "embedded_run_2",
        "embedded_run_3",
    ]
    assert [
        next(metric.value for metric in run.metrics if metric.accession == "MS:4000060")
        for run in runs
    ] == [100, 200, 300]


def test_parse_multiple_files_keeps_runs_independent():
    """Keep source identity intact when aggregating separate mzQC files."""
    used = set()
    runs = []
    for path in sorted(FIXTURE_DIR.glob("run[12].mzQC")):
        runs.extend(parse_mzqc_document(path, used))

    assert [run.sample_name for run in runs] == ["synthetic_run_1", "synthetic_run_2"]
    assert [run.source_path.name for run in runs] == ["run1.mzQC", "run2.mzQC"]


def test_legacy_singular_run_quality_is_supported(tmp_path):
    """Accept pmultiqc outputs that used the legacy singular runQuality key."""
    source = json.loads((FIXTURE_DIR / "run1.mzQC").read_text())
    source["mzQC"]["runQuality"] = source["mzQC"].pop("runQualities")
    path = tmp_path / "legacy.mzQC"
    path.write_text(json.dumps(source), encoding="utf-8")

    runs = parse_mzqc_document(path)

    assert [run.sample_name for run in runs] == ["synthetic_run_1"]


def test_general_statistics_exclude_two_endpoint_ranges():
    """Keep range endpoints out of the scalar-only General Statistics payload."""
    run = MzQCRun(
        sample_name="run",
        source_path=Path("run.mzQC"),
        metadata={},
        metrics=(
            MzQCMetric("MS:range", "mass error range", [1.0, 2.0]),
            MzQCMetric("MS:scalar", "number of MS2 spectra", 123),
        ),
    )

    data, headers = _metric_data([run], numeric_only=True)

    assert data == {"run": {"MS_scalar": 123}}
    assert set(headers) == {"MS_scalar"}


def test_prideqc_mass_accuracy_annotations_feed_reporting_helpers():
    """Extract ppm tolerances and support diagnostics without recomputing the estimator."""
    run = MzQCRun(
        sample_name="run.raw",
        source_path=Path("run.mzQC"),
        metadata={},
        metrics=(
            MzQCMetric(
                "QCPRIDE:PRECURSOR_PRECISION",
                "estimated precursor mass error precision",
                {"single_measurement_sigma": 1.25, "unit": "ppm"},
            ),
            MzQCMetric(
                "QCPRIDE:PRECURSOR_TOLERANCE",
                "suggested precursor search tolerance",
                {"suggested_tolerance": 7.5, "confidence": "high"},
            ),
            MzQCMetric(
                "QCPRIDE:FRAGMENT_PRECISION",
                "estimated fragment mass error precision",
                {"single_measurement_sigma": 2.0, "resolution_regime": "high-resolution"},
            ),
            MzQCMetric(
                "QCPRIDE:FRAGMENT_TOLERANCE",
                "suggested fragment search tolerance",
                {"suggested_tolerance": 12.0, "confidence": "moderate"},
            ),
            MzQCMetric(
                "QCPRIDE:MASS_ERROR_DIAGNOSTICS",
                "mass error estimator diagnostics",
                {
                    "precursor_paired_spectra": 1200,
                    "precursor_clusters_used": 300,
                    "fragment_pairs": 12000,
                    "fragment_resolution_regime": "high-resolution",
                },
            ),
        ),
    )

    data, ppm_plot = _mass_accuracy_report_data([run])

    assert data["run.raw"]["precursor_sigma_ppm"] == 1.25
    assert data["run.raw"]["precursor_tolerance_ppm"] == 7.5
    assert data["run.raw"]["fragment_tolerance_ppm"] == 12.0
    assert data["run.raw"]["precursor_repeat_pairs"] == 1200
    assert data["run.raw"]["fragment_resolution_regime"] == "high-resolution"
    assert ppm_plot == {
        "run.raw": {"Precursor tolerance": 7.5, "Fragment tolerance": 12.0}
    }


def test_prideqc_mass_shift_annotations_feed_summary_heatmap_and_table():
    """Build report visuals from the bounded v22.1 mzQC evidence only."""
    records = [
        {
            "delta_mass_da": 79.966331,
            "pair_support": 78,
            "unique_spectrum_support": 119,
            "median_spectral_similarity": 0.84,
            "classification": "putative-ptm",
            "confidence": "high-support",
            "unimod_candidates": [
                {
                    "name": "Phosphorylation",
                    "residual_da": -0.0006,
                }
            ],
        },
        {
            "delta_mass_da": 15.994915,
            "pair_support": 300,
            "unique_spectrum_support": 420,
            "median_spectral_similarity": 0.91,
            "classification": "putative-ptm",
            "confidence": "high-support",
            "unimod_candidates": [
                {
                    "name": "Oxidation or Hydroxylation",
                    "residual_da": 0.0001,
                }
            ],
        },
    ]
    run = MzQCRun(
        sample_name="run.raw",
        source_path=Path("run.mzQC"),
        metadata={},
        metrics=(
            MzQCMetric(
                "QCPRIDE:MASS_SHIFTS",
                "putative modification mass shifts",
                records,
            ),
            MzQCMetric(
                "QCPRIDE:MASS_SHIFT_DIAGNOSTICS",
                "mass shift scout diagnostics",
                {
                    "raw_recurrent_clusters": 120,
                    "reported_clusters": 2,
                    "explicit_profile_ms2": 1000,
                    "profile_peak_pick_failures": 1,
                },
            ),
        ),
    )

    report = _mass_shift_report_data([run])
    summary = _mass_shift_summary_table(report)["Cohort"]
    table_data = _mass_shift_cluster_table(report)

    assert report["classification_plot"] == {"run.raw": {"putative-ptm": 2}}
    assert len(report["heatmap_data"]["run.raw"]) == 2
    assert summary["reported_clusters"] == 2
    assert summary["raw_recurrent_clusters"] == 120
    assert summary["high_support"] == 2
    assert summary["putative_ptm"] == 2
    assert summary["profile_ms2"] == 1000
    assert len(table_data) == 2
    strongest = next(iter(table_data.values()))
    assert strongest["candidate"] == "Oxidation or Hydroxylation"
    assert strongest["pair_support"] == 300


def test_preferred_general_stat_selection_preserves_header_order():
    """Select preferred General Statistics columns in deterministic encounter order."""
    headers = {
        f"metric_{index}": {"title": f"number of spectra {index}"}
        for index in range(10)
    }

    assert MzQCModule._preferred_general_stat_keys(headers) == [
        f"metric_{index}" for index in range(8)
    ]


def test_duplicate_run_labels_are_disambiguated(tmp_path):
    """Disambiguate repeated run labels without collapsing separate runs."""
    source = json.loads((FIXTURE_DIR / "run2.mzQC").read_text())
    source["mzQC"]["runQualities"][0]["metadata"]["label"] = "same"
    first = tmp_path / "a.mzQC"
    second = tmp_path / "b.mzQC"
    first.write_text(json.dumps(source))
    second.write_text(json.dumps(source))

    used = set()
    runs = parse_mzqc_document(first, used) + parse_mzqc_document(second, used)
    assert [run.sample_name for run in runs] == ["same", "same [2]"]


def test_invalid_mzqc_is_rejected(tmp_path):
    """Reject JSON documents without an mzQC root object."""
    path = tmp_path / "bad.mzQC"
    path.write_text('{"not_mzqc": {}}')
    with pytest.raises(ValueError, match="missing mzQC root"):
        parse_mzqc_document(path)


def test_module_imports_with_multiqc():
    """Expose the module when MultiQC is installed."""
    pytest.importorskip("multiqc")
    assert MzQCModule is not None


def _run_multiqc(
    tmp_path: Path, input_path: Path
) -> tuple[subprocess.CompletedProcess[str], Path]:
    """Run the mzQC module through the real MultiQC CLI for integration tests."""
    multiqc_exe = shutil.which("multiqc")
    if not multiqc_exe:
        pytest.skip("MultiQC executable not installed")
    output = tmp_path / "report"
    result = subprocess.run(
        [
            multiqc_exe,
            "--strict",
            "--module",
            "mzqc",
            "--require-logs",
            str(input_path),
            "-o",
            str(output),
        ],
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
    )
    return result, output / "multiqc_report.html"


def test_multiqc_cli_aggregates_multiple_mzqc_documents(tmp_path):
    """Aggregate multiple mzQC files into one MultiQC report."""
    pytest.importorskip("multiqc")
    result, report = _run_multiqc(tmp_path, FIXTURE_DIR)
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert report.exists()
    html = report.read_text(errors="replace")
    assert "Acquisition method" in html
    assert "Data-dependent acquisition" in html
    assert "run1.mzQC" in html
    for sample in (
        "synthetic_run_1",
        "synthetic_run_2",
        "embedded_run_1",
        "embedded_run_2",
        "embedded_run_3",
    ):
        assert sample in html


def test_multiqc_cli_renders_prideqc_mass_accuracy_and_mass_shift_sections(tmp_path):
    """Render v23 prideQC structured mzQC annotations through the real MultiQC CLI."""
    pytest.importorskip("multiqc")
    source_dir = tmp_path / "input"
    source_dir.mkdir()
    source = source_dir / "prideqc-v23.mzQC"
    source.write_text(
        json.dumps(
            {
                "mzQC": {
                    "version": "1.0.0",
                    "creationDate": "2026-09-18T12:00:00+00:00",
                    "runQualities": [
                        {
                            "metadata": {
                                "label": "prideqc_demo.raw",
                                "inputFiles": [{"name": "prideqc_demo.raw"}],
                            },
                            "qualityMetrics": [
                                {
                                    "accession": "QCPRIDE:PRECURSOR_PRECISION",
                                    "name": "estimated precursor mass error precision",
                                    "value": {"single_measurement_sigma": 1.25, "unit": "ppm"},
                                },
                                {
                                    "accession": "QCPRIDE:PRECURSOR_TOLERANCE",
                                    "name": "suggested precursor search tolerance",
                                    "value": {
                                        "suggested_tolerance": 7.5,
                                        "confidence": "high",
                                    },
                                },
                                {
                                    "accession": "QCPRIDE:FRAGMENT_PRECISION",
                                    "name": "estimated fragment mass error precision",
                                    "value": {
                                        "single_measurement_sigma": 2.0,
                                        "resolution_regime": "high-resolution",
                                    },
                                },
                                {
                                    "accession": "QCPRIDE:FRAGMENT_TOLERANCE",
                                    "name": "suggested fragment search tolerance",
                                    "value": {
                                        "suggested_tolerance": 12.0,
                                        "confidence": "moderate",
                                    },
                                },
                                {
                                    "accession": "QCPRIDE:MASS_ERROR_DIAGNOSTICS",
                                    "name": "mass error estimator diagnostics",
                                    "value": {
                                        "precursor_paired_spectra": 1200,
                                        "precursor_clusters_used": 300,
                                        "fragment_pairs": 12000,
                                        "fragment_resolution_regime": "high-resolution",
                                    },
                                },
                                {
                                    "accession": "QCPRIDE:MASS_SHIFTS",
                                    "name": "putative modification mass shifts",
                                    "value": [
                                        {
                                            "delta_mass_da": 79.966331,
                                            "pair_support": 78,
                                            "unique_spectrum_support": 119,
                                            "median_spectral_similarity": 0.84,
                                            "classification": "putative-ptm",
                                            "confidence": "high-support",
                                            "unimod_candidates": [
                                                {
                                                    "name": "Phosphorylation",
                                                    "residual_da": -0.0006,
                                                },
                                                {
                                                    "name": "O-Sulfonation",
                                                    "residual_da": 0.0089,
                                                },
                                            ],
                                        }
                                    ],
                                },
                                {
                                    "accession": "QCPRIDE:MASS_SHIFT_DIAGNOSTICS",
                                    "name": "mass shift scout diagnostics",
                                    "value": {
                                        "raw_recurrent_clusters": 120,
                                        "reported_clusters": 1,
                                        "explicit_profile_ms2": 1000,
                                        "profile_peak_pick_failures": 0,
                                    },
                                },
                            ],
                        }
                    ],
                    "controlledVocabularies": [],
                }
            }
        ),
        encoding="utf-8",
    )

    result, report = _run_multiqc(tmp_path, source_dir)

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    html = report.read_text(errors="replace")
    visible_html = unescape(html)
    for label in (
        "Estimated Search Tolerances",
        "Mass Accuracy & Tolerance Details",
        "Putative Modification Mass Shifts",
        "Mass-shift Classification by Run",
        "Recurrent Mass-shift Families",
        "Mass-shift Landscape",
        "Top Mass-shift Candidates",
        "Phosphorylation",
    ):
        assert label in visible_html

    # Check the stable section anchors separately from user-visible copy so this
    # integration test is insensitive to Jinja / MultiQC HTML entity escaping.
    for anchor in (
        "mzqc-estimated-search-tolerances",
        "mzqc-mass-accuracy-tolerance-details",
        "mzqc-putative-modification-mass-shifts",
        "mzqc-mass-shift-classification",
        "mzqc-mass-shift-family-heatmap",
        "mzqc-mass-shift-landscape",
        "mzqc-top-mass-shift-candidates",
    ):
        assert anchor in html


def test_multiqc_cli_aggregates_multiple_run_quality_objects(tmp_path):
    """Render all runQuality objects from a single multi-run document."""
    pytest.importorskip("multiqc")
    result, report = _run_multiqc(tmp_path, FIXTURE_DIR / "multi_run.mzQC")
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert report.exists()
    html = report.read_text(errors="replace")
    assert "embedded_run_1" in html
    assert "embedded_run_2" in html
    assert "embedded_run_3" in html


def test_real_prideqc_fixtures_can_be_run_when_provided(tmp_path):
    """Run the real three-file integration only when the user supplies the local fixtures."""
    fixture_dir = os.environ.get("PRIDEQC_MZQC_FIXTURE_DIR")
    multiqc_exe = shutil.which("multiqc")
    if not fixture_dir or not multiqc_exe:
        pytest.skip(
            "set PRIDEQC_MZQC_FIXTURE_DIR and install MultiQC to run the real integration test"
        )

    source_dir = Path(fixture_dir).expanduser()
    expected = {
        "Natalia_TMT0_07_120m_1pt5.raw.mzQC",
        "Prosser_1004.raw.mzQC",
        "TDM_M1808_198.raw.mzQC",
    }
    present = {path.name for path in source_dir.glob("*.mzQC")}
    assert expected <= present

    result = subprocess.run(
        [
            multiqc_exe,
            "--strict",
            "--module",
            "mzqc",
            "--require-logs",
            str(source_dir),
            "-o",
            str(tmp_path / "report"),
        ],
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    report = tmp_path / "report" / "multiqc_report.html"
    assert report.exists()
    html = report.read_text(errors="replace")
    for sample in sorted(expected):
        assert sample.replace(".mzQC", "") in html


def test_scalar_metric_columns_are_consistent_across_runs():
    """Use stable metric keys across runs with overlapping scalar terms."""
    used = set()
    runs = parse_mzqc_document(FIXTURE_DIR / "run1.mzQC", used)
    runs.extend(parse_mzqc_document(FIXTURE_DIR / "run2.mzQC", used))
    data, headers = _metric_data(runs)

    assert set(data["synthetic_run_1"]) == {
        "MS_4000053",
        "MS_4000059",
        "MS_4000060",
        "QCPRIDE_MS1_RANGE_min",
        "QCPRIDE_MS1_RANGE_max",
        "QCPRIDE_PRECISION",
    }
    assert set(data["synthetic_run_2"]) == {
        "MS_4000053",
        "MS_4000059",
        "MS_4000060",
        "QCPRIDE_PRECISION",
    }
    assert data["synthetic_run_1"]["MS_4000059"] == 12000
    assert data["synthetic_run_2"]["MS_4000059"] == 12100
    assert data["synthetic_run_1"]["MS_4000060"] == 48000
    assert data["synthetic_run_2"]["MS_4000060"] == 48200
    assert data["synthetic_run_1"]["QCPRIDE_MS1_RANGE_min"] == 350.0
    assert data["synthetic_run_1"]["QCPRIDE_MS1_RANGE_max"] == 1800.0
    assert headers["MS_4000059"]["title"] == "number of MS1 spectra"
