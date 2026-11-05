"""Regression tests for pipeline-specific mzQC metadata (issue #730)."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("mzqc")

from pmultiqc.export.mzqc_exporter import MzQcExporter


@pytest.mark.parametrize("pipeline", ["DIA-NN", "MaxQuant", "QuantMS"])
def test_three_parameter_call_exports_nonempty_metrics(tmp_path, pipeline):
    exporter = MzQcExporter(pipeline, {}, str(tmp_path))
    assert exporter.export_to_file([]) == ""
    assert not (tmp_path / "qc_metrics.mzQC").exists()
    metric = SimpleNamespace(
        accession="MS:1003327", name="number of identified protein groups", value=7
    )
    path = exporter.export_to_file([metric])
    run = json.loads(Path(path).read_text())["mzQC"]["runQualities"][0]
    assert run["qualityMetrics"][0]["value"] == 7


def export_run(tmp_path, pipeline="DIA-NN", **kwargs):
    exporter = MzQcExporter(pipeline, {}, str(tmp_path), **kwargs)
    metric = SimpleNamespace(
        accession="MS:1003327", name="number of identified protein groups", value=7
    )
    path = exporter.export_to_file([metric])
    return json.loads(Path(path).read_text(encoding="utf-8"))["mzQC"]["runQualities"][0]


def test_diann_is_not_exported_as_maxquant(tmp_path):
    run = export_run(tmp_path)
    software = run["metadata"]["analysisSoftware"][0]
    assert software["name"] == "DIA-NN"
    assert software["accession"] == "MS:1003253"


@pytest.mark.parametrize(
    "pipeline, accession, name, uri",
    [
        ("DIA-NN", "MS:1003253", "DIA-NN", "https://github.com/vdemichev/DiaNN"),
        ("diann", "MS:1003253", "DIA-NN", "https://github.com/vdemichev/DiaNN"),
        ("MaxQuant", "MS:1001583", "MaxQuant", "https://www.maxquant.org"),
        ("QuantMS", "MS:1003425", "quantms", "https://quantms.org"),
    ],
)
def test_pipeline_identity_and_version(tmp_path, pipeline, accession, name, uri):
    run = export_run(tmp_path, pipeline, software_version="2.1.0")
    software = run["metadata"]["analysisSoftware"][0]
    assert software["name"] == name
    assert software["accession"] == accession
    assert software["version"] == "2.1.0"
    assert software["uri"] == uri
    assert run["qualityMetrics"][0]["value"] == 7


@pytest.mark.parametrize("version", [None, "", "  "])
def test_missing_version_is_explicit(tmp_path, caplog, version):
    software = export_run(tmp_path, software_version=version)["metadata"]["analysisSoftware"][0]
    assert software["version"] == "unknown"
    assert "version is unavailable" in caplog.text


def test_unsupported_pipeline_does_not_write_output(tmp_path):
    with pytest.raises(ValueError, match="Unsupported mzQC pipeline"):
        export_run(tmp_path, "not-a-supported-workflow")
    assert not (tmp_path / "qc_metrics.mzQC").exists()


def test_failed_replacement_preserves_previous_file(tmp_path, monkeypatch, capsys):
    from pmultiqc.export import mzqc_exporter

    target = tmp_path / "qc_metrics.mzQC"
    target.write_text("previous export", encoding="utf-8")

    def fail_replace(*args):
        raise PermissionError("simulated output failure")

    monkeypatch.setattr(mzqc_exporter.os, "replace", fail_replace)
    with pytest.raises(OSError, match="simulated output failure"):
        export_run(tmp_path)
    assert target.read_text(encoding="utf-8") == "previous export"
    assert "Successfully wrote" not in capsys.readouterr().out
    assert not list(tmp_path.glob("*.tmp"))


@pytest.mark.parametrize("version", ["1.8.1", None])
@pytest.mark.parametrize("extension", ["tsv", "parquet"])
def test_diann_module_passes_version(tmp_path, monkeypatch, version, extension):
    """Exercise the real module-to-exporter wiring without rendering plots."""
    from collections import defaultdict

    from multiqc import config

    from pmultiqc.modules.diann import diann

    for name in (
        "draw_diann_metadata_table",
        "draw_ms_information",
        "draw_summary_protein_ident_table",
        "draw_identi_num",
        "draw_num_pep_per_protein",
        "draw_identification",
        "add_group_modules",
    ):
        monkeypatch.setattr(diann, name, lambda *args, **kwargs: None)
    monkeypatch.setattr(diann, "aggregate_general_stats", lambda **kwargs: {})
    monkeypatch.setattr(
        diann,
        "parse_diann_report",
        lambda **kwargs: (1, 2, {}, {}, [], {}, {}, [], {}),
    )
    monkeypatch.setattr(config, "output_dir", str(tmp_path))
    module = diann.DiannModule.__new__(diann.DiannModule)
    values = {
        "diann_version": version,
        "sub_sections": defaultdict(list),
        "ms1_general_stats": {},
        "current_sum_by_run": {},
        "file_df": None,
        "sample_df": None,
        "ms1_tic": {},
        "ms1_bpc": {},
        "ms1_peaks": {},
        "diann_report_path": str(tmp_path / f"report.{extension}"),
        "heatmap_color_list": [],
        "ms_with_psm": [],
        "ms_paths": [],
        "enable_dia": True,
        "enable_sdrf": False,
        "is_multi_conditions": False,
        "ms_info_path": [],
        "long_trends": None,
    }
    for name, value in values.items():
        setattr(module, name, value)
    module.draw_plots()
    document = json.loads((tmp_path / "diann_qc.mzQC").read_text(encoding="utf-8"))
    metadata = document["mzQC"]["runQualities"][0]["metadata"]
    assert metadata["analysisSoftware"][0]["name"] == "DIA-NN"
    assert metadata["analysisSoftware"][0]["version"] == (version or "unknown")


@pytest.mark.parametrize("version_value", ["1.5.2.8", None, "", "NA", "NaN", "  "])
def test_maxquant_real_parameter_parser_reaches_export(
    tmp_path, monkeypatch, caplog, version_value
):
    import gzip
    import logging

    from multiqc import config

    from pmultiqc.modules.maxquant import maxquant

    source = Path(__file__).parent / "resources/maxquant/PXD003133-22min/parameters.txt.gz"
    text = gzip.decompress(source.read_bytes()).decode()
    lines = [line for line in text.splitlines() if not line.startswith("Version\t")]
    if version_value is not None:
        lines.append(f"Version\t{version_value}")
    text = "\n".join(lines)
    parameters = tmp_path / "parameters.txt"
    parameters.write_text(text)
    module = maxquant.MaxQuantModule.__new__(maxquant.MaxQuantModule)
    module.log = logging.getLogger("maxquant-test")
    module.sub_sections = {"experiment": []}
    module.find_log_files = lambda key, **kwargs: (
        [{"root": str(tmp_path), "fn": "parameters.txt"}]
        if key.endswith("maxquant_result")
        else []
    )
    monkeypatch.setattr(config, "output_dir", str(tmp_path))
    monkeypatch.setattr(
        maxquant.maxquant_plots, "draw_parameters", lambda *args, **kwargs: None, raising=False
    )
    for method in (
        "_process_sdrf_file",
        "_process_protein_groups_file",
        "_process_summary_file",
        "_process_evidence_file",
        "_process_msms_file",
        "_process_msms_scans_file",
        "_calculate_heatmap",
    ):
        monkeypatch.setattr(module, method, lambda *args: {})
    module.get_data()
    run = json.loads((tmp_path / "maxquant_qc.mzQC").read_text())["mzQC"]["runQualities"][0]
    assert run["metadata"]["analysisSoftware"][0]["version"] == (
        "1.5.2.8" if version_value == "1.5.2.8" else "unknown"
    )
    assert run["metadata"]["analysisSoftware"][0]["name"] == "MaxQuant"
    assert run["qualityMetrics"]
    if version_value != "1.5.2.8":
        assert "MaxQuant software version is unavailable" in caplog.text
