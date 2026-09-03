"""Tests for DIA-NN module functionality."""

import os
from pathlib import Path

import pytest

from pmultiqc.modules.common.dia_utils import parse_diann_version

TEST_DATA_DIR = Path(os.path.dirname(__file__)) / "resources" / "diann"


class TestDiann:
    """Tests for DIA-NN module functions."""

    def test_parse_diann_version_from_log(self):
        """Test parsing DIA-NN version from log file."""
        log_file = TEST_DATA_DIR / "report.log.txt"
        assert log_file.exists(), f"Test log file {log_file} does not exist"

        version = parse_diann_version(str(log_file))
        assert version is not None, "DIA-NN version should be parsed from log file"
        assert version == "1.8.1", f"Expected version '1.8.1', got '{version}'"

    def test_parse_diann_version_nonexistent_file(self):
        """Test that parse_diann_version returns None for nonexistent file."""
        version = parse_diann_version("/nonexistent/path/report.log.txt")
        assert version is None, "Should return None for nonexistent file"

    def test_parse_diann_version_no_version_line(self, tmp_path):
        """Test that parse_diann_version returns None when no version line is found."""
        log_file = tmp_path / "report.log.txt"
        log_file.write_text("This log file has no version line.\nAnother line.\n")

        version = parse_diann_version(str(log_file))
        assert version is None, "Should return None when no version line is present"

    def test_parse_diann_version_different_formats(self, tmp_path):
        """Test parsing DIA-NN version with different version string formats."""
        log_file = tmp_path / "report.log.txt"
        log_file.write_text("DIA-NN 2.0 (Data-Independent Acquisition by Neural Networks)\n")

        version = parse_diann_version(str(log_file))
        assert version == "2.0", f"Expected version '2.0', got '{version}'"


class TestDiaReportMemory:
    """The DIA plotting entry points must not duplicate the whole DIA-NN report.

    On large experiments the report has tens of millions of rows and ~50 columns,
    several of them object-dtype strings. Copying it wholesale to add one derived
    column is what makes these steps run out of memory, so the frames handed to the
    plotting helpers must carry only the columns those helpers read.
    """

    @staticmethod
    def _report_df():
        """A DIA-NN-shaped report with extra columns none of these plots use."""
        import pandas as pd

        return pd.DataFrame(
            {
                "Run": ["run1", "run1", "run2", "run2"],
                "Modified.Sequence": ["PEPTIDEK", "PEPTIDER", "PEPTIDEK", "PEPTIDER"],
                "Protein.Group": ["P1", "P2", "P1", "P2"],
                "Precursor.Quantity": [100.0, 200.0, 0.0, 400.0],
                "RT": [10.0, 20.0, 11.0, 21.0],
                "iRT": [9.0, 19.0, 10.0, 20.0],
                "Predicted.RT": [10.5, 19.5, 11.5, 20.5],
                "RT.Start": [9.5, 19.5, 10.5, 20.5],
                "RT.Stop": [10.5, 20.5, 11.5, 21.5],
                "FWHM": [0.5, 0.6, 0.5, 0.6],
                "Normalisation.Factor": [1.0, 1.1, 1.0, 1.1],
                # Columns that must never reach the plotting helpers:
                "Stripped.Sequence": ["PEPTIDEK", "PEPTIDER", "PEPTIDEK", "PEPTIDER"],
                "Precursor.Id": ["PEPTIDEK2", "PEPTIDER2", "PEPTIDEK2", "PEPTIDER2"],
                "Ms1.Area": [1.0, 2.0, 3.0, 4.0],
            }
        )

    def test_draw_dia_intensitys_narrows_columns(self, monkeypatch):
        """Only Run/Modified.Sequence/Protein.Group/intensity reach the plots."""
        from pmultiqc.modules.common import dia_utils
        import pandas as pd

        captured = []

        def capture(sub_section, df, sdrf_file_df):
            captured.append(df)

        monkeypatch.setattr(dia_utils.dia_plots, "draw_dia_intensity_dis", capture)
        monkeypatch.setattr(dia_utils.dia_plots, "draw_dia_intensity_std", capture)

        report_df = self._report_df()
        dia_utils.draw_dia_intensitys(None, report_df, pd.DataFrame())

        assert len(captured) == 2, "both intensity plots should be called"
        for df in captured:
            assert set(df.columns) == {
                "Run",
                "Modified.Sequence",
                "Protein.Group",
                "Precursor.Quantity",
                "log_intensity",
            }, f"unexpected columns carried into the plot: {sorted(df.columns)}"
            # Rows with no quantity are still filtered out.
            assert (df["Precursor.Quantity"] > 0).all()
            assert len(df) == 3

    def test_draw_dia_intensitys_does_not_mutate_report(self, monkeypatch):
        """The caller's report must not gain a log_intensity column."""
        from pmultiqc.modules.common import dia_utils
        import pandas as pd

        monkeypatch.setattr(
            dia_utils.dia_plots, "draw_dia_intensity_dis", lambda *a, **k: None
        )
        monkeypatch.setattr(
            dia_utils.dia_plots, "draw_dia_intensity_std", lambda *a, **k: None
        )

        report_df = self._report_df()
        before = list(report_df.columns)

        dia_utils.draw_dia_intensitys(None, report_df, pd.DataFrame())

        assert list(report_df.columns) == before

    def test_draw_dia_rt_qc_narrows_columns(self, monkeypatch):
        """RT QC copies only the RT-related columns, not the whole report."""
        from pmultiqc.modules.common import dia_utils

        captured = []

        def capture_ids_rt(sub_section, report_df):
            # Snapshot the columns now: draw_dia_rt_qc adds derived columns to this
            # same frame after this call, and we care about what it started from.
            captured.append(list(report_df.columns))

        monkeypatch.setattr(dia_utils, "draw_dia_ids_rt", capture_ids_rt)
        monkeypatch.setattr(dia_utils, "cal_feature_avg_rt", lambda df, col: None)
        monkeypatch.setattr(dia_utils, "cal_rt_irt_loess", lambda df: None)

        dia_utils.draw_dia_rt_qc(None, self._report_df())

        assert len(captured) == 1
        assert set(captured[0]) <= set(dia_utils.RT_QC_COLUMNS)
        assert "Stripped.Sequence" not in captured[0]
        assert "Precursor.Id" not in captured[0]

    def test_draw_dia_rt_qc_tolerates_missing_optional_columns(self, monkeypatch):
        """A report without the optional RT columns still works."""
        from pmultiqc.modules.common import dia_utils

        monkeypatch.setattr(dia_utils, "draw_dia_ids_rt", lambda *a, **k: None)
        monkeypatch.setattr(dia_utils, "cal_feature_avg_rt", lambda df, col: None)
        monkeypatch.setattr(dia_utils, "cal_rt_irt_loess", lambda df: None)

        minimal = self._report_df()[["Run", "RT"]]
        dia_utils.draw_dia_rt_qc(None, minimal)
