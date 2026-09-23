"""Tests for direct mzQC -> MultiQC input handling."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from html import unescape
from pathlib import Path

import pytest

from pmultiqc.modules.mzqc import MzQCMetric, MzQCModule, MzQCRun, parse_mzqc_document
from pmultiqc.modules.mzqc.mzqc import (
    _beta_prevalence_probability,
    _candidate_family_counts,
    _experiment_group_pca_data,
    _mass_accuracy_detail_description,
    _mass_accuracy_report_data,
    _mass_shift_cluster_table,
    _mass_shift_family_context,
    _mass_shift_report_data,
    _mass_shift_summary_description,
    _mass_shift_summary_table,
    _metric_data,
    _metric_display_title,
    _prevalence_probability_histogram,
    _run_group_assignments,
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

    assert report["classification_plot"] == {"run.raw": {"PTM-compatible": 2}}
    assert set(report["heatmap_data"]["run.raw"]) == {"+15.995 Da", "+79.966 Da"}
    assert summary["reported_clusters"] == 2
    assert summary["raw_recurrent_clusters"] == 120
    assert summary["high_support"] == 2
    assert summary["putative_ptm"] == 2
    assert summary["profile_ms2"] == 1000
    assert len(table_data) == 2
    strongest = next(iter(table_data.values()))
    assert strongest["candidate"] == "Oxidation or Hydroxylation"
    assert strongest["classification"] == "PTM-compatible"
    assert strongest["pair_support"] == 300


def test_mass_shift_heatmap_is_bounded_and_uses_mass_only_axis_labels():
    """Keep dense candidate chemistry out of the heatmap axis and cap it at 12 families."""
    records = [
        {
            "delta_mass_da": float(index) + 0.123,
            "pair_support": 1000 - index,
            "unique_spectrum_support": 20,
            "median_spectral_similarity": 0.8,
            "classification": "putative-modification",
            "confidence": "moderate",
            "unimod_candidates": [{"name": f"Candidate {index}", "residual_da": 0.001}],
        }
        for index in range(15)
    ]
    run = MzQCRun(
        sample_name="run.raw",
        source_path=Path("run.mzQC"),
        metadata={},
        metrics=(
            MzQCMetric("QCPRIDE:MASS_SHIFTS", "putative modification mass shifts", records),
        ),
    )

    report = _mass_shift_report_data([run])
    labels = list(report["heatmap_data"]["run.raw"])

    assert len(labels) == 12
    assert all(label.endswith(" Da") for label in labels)
    assert all("Candidate" not in label for label in labels)
    assert report["classification_plot"] == {
        "run.raw": {"Modification-compatible": 15}
    }


def test_prideqc_no_data_messages_make_abstention_explicit():
    """Explain evidence-free cohorts rather than making intentionally absent plots look broken."""
    accuracy_data = {"run.raw": {"fragment_sigma_ppm": 140.0}}
    assert "No runs met the estimator criteria" in _mass_accuracy_detail_description(
        accuracy_data
    )

    run = MzQCRun(
        sample_name="run.raw",
        source_path=Path("run.mzQC"),
        metadata={},
        metrics=(
            MzQCMetric(
                "QCPRIDE:MASS_SHIFT_DIAGNOSTICS",
                "mass shift scout diagnostics",
                {"raw_recurrent_clusters": 0, "reported_clusters": 0},
            ),
        ),
    )
    report = _mass_shift_report_data([run])
    assert "No recurrent modification-compatible mass-shift clusters" in (
        _mass_shift_summary_description(report)
    )


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
    tmp_path: Path, input_path: Path, output_path: Path | None = None
) -> tuple[subprocess.CompletedProcess[str], Path]:
    """Run the mzQC module through the real MultiQC CLI for integration tests."""
    multiqc_exe = shutil.which("multiqc")
    if not multiqc_exe:
        pytest.skip("MultiQC executable not installed")
    output = output_path or tmp_path / "report"
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
        "Fragment tolerance (ppm)",
        "PTM-compatible",
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


def test_prideqc_mzqc_dataset_can_be_run_when_provided(tmp_path):
    """Run the full PRIDE QC mzQC dataset when CI provides the extracted archive."""
    fixture_dir = os.environ.get("PRIDEQC_MZQC_DATASET_DIR")
    if not fixture_dir:
        pytest.skip("set PRIDEQC_MZQC_DATASET_DIR to run the PRIDE QC dataset test")

    source_dir = Path(fixture_dir).expanduser()
    mzqc_files = sorted(source_dir.rglob("*.mzQC"))
    assert len(mzqc_files) == 67

    for path in mzqc_files:
        with path.open() as handle:
            json.load(
                handle,
                parse_constant=lambda value: (_ for _ in ()).throw(
                    ValueError(f"non-standard JSON constant {value}")
                ),
            )

    output_dir = os.environ.get("PRIDEQC_MZQC_OUTPUT_DIR")
    result, report = _run_multiqc(
        tmp_path,
        source_dir,
        Path(output_dir) if output_dir else None,
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert report.exists()
    html = report.read_text(errors="replace")
    assert "Run Overview" in html
    assert mzqc_files[0].stem in html
    assert mzqc_files[-1].stem in html


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


def _synthetic_tolerance_run(
    name: str, fragment_ppm: float, rt_seconds: float = 5400.0
) -> MzQCRun:
    return MzQCRun(
        sample_name=name,
        source_path=Path(f"{name}.mzQC"),
        metadata={},
        metrics=(
            MzQCMetric(
                "QCPRIDE:PRECURSOR_TOLERANCE",
                "suggested precursor search tolerance",
                {"suggested_tolerance": 8.0, "confidence": "high"},
            ),
            MzQCMetric(
                "QCPRIDE:FRAGMENT_TOLERANCE",
                "suggested fragment search tolerance",
                {"suggested_tolerance": fragment_ppm, "confidence": "high"},
            ),
            MzQCMetric(
                "QCPRIDE:MASS_ERROR_DIAGNOSTICS",
                "mass error estimator diagnostics",
                {"fragment_resolution_regime": "high-resolution"},
            ),
            MzQCMetric("MS:4000053", "chromatography duration", rt_seconds),
            MzQCMetric("QCPRIDE:ISO", "IsolationWidth_MS2_Median", 1.6),
        ),
    )


def test_run_group_detection_finds_supported_two_population_tolerance_split():
    runs = [
        *[_synthetic_tolerance_run(f"low-{index}", 6.0 + index * 0.1) for index in range(6)],
        *[_synthetic_tolerance_run(f"high-{index}", 21.0 + index * 0.5) for index in range(6)],
    ]

    assignments, summaries, evidence = _run_group_assignments(runs)

    assert len(summaries) == 2
    assert len(set(assignments[name] for name in assignments if name.startswith("low-"))) == 1
    assert len(set(assignments[name] for name in assignments if name.startswith("high-"))) == 1
    assert assignments["low-0"] != assignments["high-0"]
    assert any("fragment tolerance" in item for item in evidence)
    assert summaries[assignments["low-0"]]["fragment_common"] == pytest.approx(6.5)
    assert summaries[assignments["high-0"]]["fragment_common"] == pytest.approx(23.5)


def test_run_group_detection_does_not_promote_singleton_outlier():
    runs = [
        *[_synthetic_tolerance_run(f"base-{index}", 12.0) for index in range(9)],
        _synthetic_tolerance_run("outlier", 40.0),
    ]

    assignments, summaries, evidence = _run_group_assignments(runs)

    assert len(summaries) == 1
    assert set(assignments.values()) == {"Experiment group 1"}
    assert evidence == []


def test_experiment_group_detection_uses_chromatography_duration():
    runs = [
        *[_synthetic_tolerance_run(f"short-{index}", 12.0, 5400.0) for index in range(4)],
        *[_synthetic_tolerance_run(f"long-{index}", 12.0, 7200.0) for index in range(4)],
    ]

    assignments, summaries, evidence = _run_group_assignments(runs)

    assert len(summaries) == 2
    assert assignments["short-0"] != assignments["long-0"]
    assert any("chromatography duration" in item for item in evidence)


def test_experiment_group_detection_can_refine_tolerance_group_by_rt():
    runs = [
        *[_synthetic_tolerance_run(f"wide-{index}", 24.0, 5400.0) for index in range(6)],
        *[_synthetic_tolerance_run(f"narrow-{index}", 7.0, 5400.0) for index in range(6)],
        *[_synthetic_tolerance_run(f"long-{index}", 7.0, 7200.0) for index in range(3)],
    ]

    assignments, summaries, evidence = _run_group_assignments(runs)

    assert len(summaries) == 3
    assert assignments["wide-0"] != assignments["narrow-0"]
    assert assignments["narrow-0"] != assignments["long-0"]
    assert any("chromatography duration" in item for item in evidence)
    assert any("fragment tolerance" in item for item in evidence)


def test_experiment_group_pca_uses_multifeature_evidence_and_group_colours():
    runs = [
        *[_synthetic_tolerance_run(f"low-{index}", 6.0, 5400.0) for index in range(4)],
        *[_synthetic_tolerance_run(f"high-{index}", 22.0, 7200.0) for index in range(4)],
    ]
    assignments, _, _ = _run_group_assignments(runs)

    plot_data, features, variance = _experiment_group_pca_data(runs, assignments)

    assert set(plot_data) == set(assignments.values())
    assert {"fragment", "rt"}.issubset(features)
    assert variance is not None
    assert sum(len(points) for points in plot_data.values()) == 8
    colours = {
        group: {point["color"] for point in points}
        for group, points in plot_data.items()
    }
    assert all(len(values) == 1 for values in colours.values())
    assert len({next(iter(values)) for values in colours.values()}) == len(colours)


def _synthetic_mass_shift_run(name: str, include_family: bool) -> MzQCRun:
    records = []
    if include_family:
        records.append(
            {
                "delta_mass_da": 79.96633,
                "pair_support": 80,
                "unique_spectrum_support": 40,
                "classification": "putative-ptm",
                "confidence": "high-support",
                "unimod_candidates": [{"name": "Phosphorylation", "residual_da": 0.0}],
            }
        )
    return MzQCRun(
        sample_name=name,
        source_path=Path(f"{name}.mzQC"),
        metadata={},
        metrics=(
            MzQCMetric(
                "QCPRIDE:MASS_SHIFTS",
                "putative modification mass shifts",
                records,
            ),
        ),
    )


def test_high_support_modification_family_requires_ninety_percent_of_group_runs():
    assignments = {f"run-{index}": "Experiment group 1" for index in range(10)}
    eight_hits = [
        _synthetic_mass_shift_run(f"run-{index}", index < 8) for index in range(10)
    ]
    _, all_eight, high_eight = _mass_shift_family_context(eight_hits, assignments)
    assert len(all_eight) == 1
    assert all_eight[0]["run_prevalence"] == pytest.approx(0.8)
    assert high_eight == []

    nine_hits = [
        _synthetic_mass_shift_run(f"run-{index}", index < 9) for index in range(10)
    ]
    _, all_nine, high_nine = _mass_shift_family_context(nine_hits, assignments)
    assert len(all_nine) == 1
    assert all_nine[0]["run_prevalence"] == pytest.approx(0.9)
    assert all_nine[0]["high_support_run_fraction"] == pytest.approx(1.0)
    assert len(high_nine) == 1

    moderate_runs = []
    for index in range(10):
        run = _synthetic_mass_shift_run(f"run-{index}", index < 9)
        if index < 2:
            record = run.metrics[0].value[0]
            record["confidence"] = "moderate"
        moderate_runs.append(run)
    _, _, high_moderate = _mass_shift_family_context(moderate_runs, assignments)
    assert high_moderate == []


def test_candidate_family_counts_bounds_pie_categories():
    rows = [
        {"candidate": f"Candidate {index}"}
        for index in range(10)
    ] + [{"candidate": "Candidate 0"}]

    counts = _candidate_family_counts(rows, limit=3)

    assert len(counts) == 4
    assert counts["Candidate 0"] == 2
    assert counts["Other"] == 7


def test_prevalence_probability_histogram_spans_zero_to_one_evenly():
    rows = [
        {"raw_prevalence_probability": 0.02},
        {"raw_prevalence_probability": 0.18},
        {"raw_prevalence_probability": 0.52},
        {"raw_prevalence_probability": 0.92},
        {"raw_prevalence_probability": 1.0},
    ]

    histogram = _prevalence_probability_histogram(rows)

    assert list(histogram) == [
        "0.0–0.1",
        "0.1–0.2",
        "0.2–0.3",
        "0.3–0.4",
        "0.4–0.5",
        "0.5–0.6",
        "0.6–0.7",
        "0.7–0.8",
        "0.8–0.9",
        "0.9–1.0",
    ]
    assert histogram["0.0–0.1"]["Families"] == 1
    assert histogram["0.1–0.2"]["Families"] == 1
    assert histogram["0.5–0.6"]["Families"] == 1
    assert histogram["0.9–1.0"]["Families"] == 2


def test_beta_prevalence_probability_matches_v4_single_run_boundary():
    assert _beta_prevalence_probability(1, 1) == pytest.approx(0.99)
    assert _beta_prevalence_probability(0, 1) == pytest.approx(0.81)


def test_cv_backed_metric_display_titles_are_compact():
    isolation = MzQCMetric("QCPRIDE:ISO", "IsolationWidth_MS2_Median", 1.6)
    scan = MzQCMetric("QCPRIDE:SCAN", "ScanWindow_MS1", [350.0, 1800.0])

    assert _metric_display_title(isolation) == "MS2 isolation width — median"
    assert _metric_display_title(scan, "lower") == "MS1 scan window lower limit"
    assert _metric_display_title(scan, "upper") == "MS1 scan window upper limit"
