"""Software-only export checks with real synthetic DIA/MS-info file parsing.

Only plotting and report assembly are suppressed. These are module call-chain
tests with synthetic inputs, not full MultiQC CLI or biological-data runs.
"""

import json
from collections import defaultdict

import pandas as pd
import pytest

pytest.importorskip("mzqc")


def _write_inputs(tmp_path, extension, with_version, with_ms_info):
    report = tmp_path / f"report.{extension}"
    frame = pd.DataFrame(
        {
            "Run": ["r1", "r1", "r2", "decoy"],
            "Modified.Sequence": ["PEPTIDEK", "PEPTIDER", "PEPTIDEK", "DECOYK"],
            "Stripped.Sequence": ["PEPTIDEK", "PEPTIDER", "PEPTIDEK", "DECOYK"],
            "Protein.Group": ["P1", "P2", "P1", "DECOY"],
            "Protein.Names": ["N1", "N2", "N1", "DECOY"],
            "Precursor.Quantity": [100.0, 200.0, 300.0, 999.0],
            "Precursor.Normalised": [100.0, 200.0, 300.0, 999.0],
            "Precursor.Charge": [2, 2, 3, 2],
            "RT": [10.0, 20.0, 11.0, 12.0],
            "Decoy": [0, 0, 0, 1],
            "Proteotypic": [1, 1, 1, 1],
        }
    )
    if extension == "tsv":
        frame.to_csv(report, sep="\t", index=False)
    else:
        pytest.importorskip("pyarrow")
        frame.to_parquet(report, index=False)
    found = {f"pmultiqc/diann_report_{extension}": [report]}
    if with_version:
        log = tmp_path / "report.log.txt"
        log.write_text(
            "DIA-NN 2.2.7 (Data-Independent Acquisition by Neural Networks)\n",
            encoding="utf-8",
        )
        found["pmultiqc/diann_log_txt"] = [log]
    if with_ms_info:
        pytest.importorskip("pyarrow")
        ms_info = tmp_path / "r1_ms_info.parquet"
        pd.DataFrame(
            {
                "ms_level": [1, 2, 2],
                "precursor_charge": [0, 2, 3],
                "rt": [60.0, 61.0, 62.0],
                "summed_peak_intensities": [1000.0, 100.0, 200.0],
                "num_peaks": [50, 20, 30],
                "base_peak_intensity": [500.0, 50.0, 100.0],
                "precursor_intensity": [0.0, 50.0, 100.0],
                "acquisition_datetime": ["2026-01-01T00:00:00Z"] * 3,
            }
        ).to_parquet(ms_info, index=False)
        found["pmultiqc/ms_info"] = [ms_info]

    def find_files(key, **kwargs):
        return [{"root": str(p.parent), "fn": p.name} for p in found.get(key, [])]

    return find_files


@pytest.mark.parametrize("pipeline", ["diann", "quantms"])
@pytest.mark.parametrize("extension", ["tsv", "parquet"])
@pytest.mark.parametrize("with_version", [False, True])
@pytest.mark.parametrize("with_ms_info", [False, True])
def test_real_dia_inputs_keep_software_identity_and_metrics(
    tmp_path, monkeypatch, caplog, pipeline, extension, with_version, with_ms_info
):
    from multiqc import config

    from pmultiqc.export.mzqc_exporter import MzQcExporter, to_json_safe
    from pmultiqc.modules.common import dia_utils
    from pmultiqc.modules.diann import diann
    from pmultiqc.modules.quantms import quantms

    module_source = diann if pipeline == "diann" else quantms
    # Readers, preprocessing, statistics and version parsing all run unchanged.
    monkeypatch.setattr(dia_utils, "_draw_diann_plots", lambda *args: None)
    for name in (
        "draw_diann_metadata_table",
        "draw_ms_information",
        "draw_summary_protein_ident_table",
        "draw_identi_num",
        "draw_num_pep_per_protein",
        "draw_identification",
        "draw_long_trends",
        "draw_peptide_length_distribution",
        "draw_peaks_per_ms2",
        "draw_peak_intensity_distribution",
        "add_group_modules",
    ):
        monkeypatch.setattr(module_source, name, lambda *args, **kwargs: None)
    monkeypatch.setattr(config, "output_dir", str(tmp_path))
    monkeypatch.setattr(
        config,
        "kwargs",
        {"quantification_method": "LFQ", "disable_table": True, "ignored_idxml": True},
    )
    module_class = diann.DiannModule if pipeline == "diann" else quantms.QuantMSModule
    module = module_class(
        _write_inputs(tmp_path, extension, with_version, with_ms_info),
        defaultdict(list),
        [],
    )
    if pipeline == "quantms":
        # These methods only render extra QuantMS sections after parsing.
        for name in (
            "draw_quantms_contaminants",
            "draw_quantms_msms_section",
            "draw_quantms_time_section",
        ):
            monkeypatch.setattr(module, name, lambda *args, **kwargs: None)

    exported_metrics = []
    metric_objects = []
    original_export = MzQcExporter.export_to_file

    def capture_export(self, metrics, **kwargs):
        metric_objects.extend(metrics)
        exported_metrics.extend(
            {
                "accession": metric.accession,
                "name": metric.name,
                "value": to_json_safe(self.sanitize_value(metric.value)),
            }
            for metric in metrics
        )
        return original_export(self, metrics, **kwargs)

    monkeypatch.setattr(MzQcExporter, "export_to_file", capture_export)
    assert module.get_data()
    assert module.diann_version == ("2.2.7" if with_version else None)
    assert module.read_ms_info == with_ms_info
    if with_ms_info:
        assert module.mzml_table["r1"]["MS1_Num"] == 1
        assert module.mzml_table["r1"]["MS2_Num"] == 2
    module.draw_plots()
    assert module.total_peptide_count == 2
    assert module.total_protein_quantified == 2
    assert set(module.ms_with_psm) == {"r1", "r2"}
    document = json.loads((tmp_path / f"{pipeline}_qc.mzQC").read_text())
    run = document["mzQC"]["runQualities"][0]
    software = run["metadata"]["analysisSoftware"][0]
    expected = (
        ("DIA-NN", "MS:1003253", "2.2.7" if with_version else "unknown")
        if pipeline == "diann"
        else ("quantms", "MS:1003425", "unknown")
    )
    assert (software["name"], software["accession"], software["version"]) == expected
    assert software["uri"].startswith("https://")
    assert run["qualityMetrics"] == exported_metrics
    assert exported_metrics
    # Changing only the optional software version must not alter these metrics.
    comparison = MzQcExporter(
        "DIA-NN" if pipeline == "diann" else "QuantMS",
        {},
        str(tmp_path / "comparison"),
        software_version="9.8.7",
    )
    comparison_path = original_export(comparison, metric_objects)
    comparison_run = json.loads(open(comparison_path, encoding="utf-8").read())["mzQC"][
        "runQualities"
    ][0]
    assert comparison_run["qualityMetrics"] == run["qualityMetrics"]
    assert comparison_run["metadata"]["analysisSoftware"][0]["version"] == "9.8.7"
    assert "Metric extraction or export failed" not in caplog.text
