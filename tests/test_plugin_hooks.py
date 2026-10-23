"""Regression tests for MultiQC plugin lifecycle hooks."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

from multiqc import report

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


class TestLfqGroupsSearchPattern(unittest.TestCase):
    def test_lfq_groups_manifest_is_not_taken_by_the_generic_tsv_pattern(self):
        """Each file goes to the first search key that matches it, so *.tsv must not come first."""
        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(pmultiqc_main.config, "sp", {}):
            tmp_path = Path(tmp)
            (tmp_path / "lfq_groups.tsv").write_text("group_id\trun_id\n")
            (tmp_path / "psm.tsv").write_text("x\n")
            pmultiqc_main.pmultiqc_plugin_before_config()

            report.reset_file_search()
            report.analysis_files = [str(tmp_path)]
            try:
                report.search_files(["pmultiqc"])
                found = {key: [f["fn"] for f in files] for key, files in report.files.items()}
            finally:
                report.reset_file_search()

        self.assertEqual(found["pmultiqc/lfq_groups"], ["lfq_groups.tsv"])
        self.assertEqual(found["pmultiqc/tsv"], ["psm.tsv"])
