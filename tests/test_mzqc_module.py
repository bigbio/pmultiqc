"""Tests for direct mzQC -> MultiQC input handling."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from pmultiqc.modules.mzqc import MzQCMetric, MzQCModule, MzQCRun, parse_mzqc_document
from pmultiqc.modules.mzqc.mzqc import _metric_data


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
