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
