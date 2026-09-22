"""Regression tests for MultiQC plugin lifecycle hooks."""

from __future__ import annotations

from pmultiqc import main as pmultiqc_main


def test_before_config_registers_all_search_patterns(monkeypatch):
    """Discovery defaults must exist before MultiQC builds its file index."""
    monkeypatch.setattr(pmultiqc_main.config, "sp", {})
    monkeypatch.setattr(pmultiqc_main.config, "log_filesize_limit", 50_000_000)
    monkeypatch.setattr(pmultiqc_main.config, "thousandsSep_format", ",")

    pmultiqc_main.pmultiqc_plugin_before_config()

    assert set(pmultiqc_main.PMULTIQC_SEARCH_PATTERNS) <= set(pmultiqc_main.config.sp)
    assert (
        pmultiqc_main.config.log_filesize_limit
        == pmultiqc_main.PMULTIQC_LOG_FILESIZE_LIMIT
    )
    assert pmultiqc_main.config.thousandsSep_format == ""
    assert pmultiqc_main.config.sp["mzqc"]["fn"] == "*.mzQC"
    assert pmultiqc_main.config.sp["pmultiqc/mztab"]["fn"] == "*.mzTab"
    assert pmultiqc_main.config.sp["pmultiqc/qpx_psm"]["fn"] == "*.psm.parquet"
    assert pmultiqc_main.config.sp["pmultiqc/diann_report_parquet"]["fn"] == "report.parquet"


def test_before_config_preserves_existing_search_pattern_overrides(monkeypatch):
    """Early defaults must not overwrite a pattern already supplied by MultiQC config."""
    custom = {"pmultiqc/mztab": {"fn": "custom.mzTab", "num_lines": 0}}
    monkeypatch.setattr(pmultiqc_main.config, "sp", custom)

    pmultiqc_main.pmultiqc_plugin_before_config()

    assert pmultiqc_main.config.sp["pmultiqc/mztab"]["fn"] == "custom.mzTab"
