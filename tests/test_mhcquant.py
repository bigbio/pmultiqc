import os
from pathlib import Path

import pytest

from pmultiqc.modules.mhcquant.mhcquant import (
    _parse_histogram_file,
    _parse_headerless_box_file,
    _parse_general_stats,
    _parse_percolator_median_weights,
)


MHCQUANT_DATA_DIR = Path("tests/resources/mhcquant/multiqc_data")


class TestMhcquantParsers:
    """Test the data parsing functions for the mhcquant module."""

    def test_general_stats_parsing(self):
        filepath = MHCQUANT_DATA_DIR / "multiqc_general_stats.txt"
        data, headers = _parse_general_stats(filepath)

        assert data is not None, "general_stats data should not be None"
        assert headers is not None, "general_stats headers should not be None"
        assert len(data) >= 1, "Should have at least one sample"

        sample_name = list(data.keys())[0]
        assert "# Peptides" in data[sample_name] or "# PSMs" in data[sample_name], (
            "Expected general stats columns not found"
        )

    def test_histogram_mz_parsing(self):
        filepath = MHCQUANT_DATA_DIR / "multiqc_histogram_mz.txt"
        data = _parse_histogram_file(filepath)

        assert data is not None, "histogram_mz data should not be None"
        assert len(data) >= 1, "Should have at least one sample"

        sample_name = list(data.keys())[0]
        assert isinstance(data[sample_name], dict), "Sample data should be a dict"
        assert all(
            isinstance(k, float) for k in data[sample_name].keys()
        ), "Keys should be float bin centers"
        assert all(
            isinstance(v, float) for v in data[sample_name].values()
        ), "Values should be float counts"

    def test_histogram_rt_parsing(self):
        filepath = MHCQUANT_DATA_DIR / "multiqc_histogram_rt.txt"
        data = _parse_histogram_file(filepath)

        assert data is not None
        assert len(data) >= 1

    def test_histogram_im_parsing(self):
        filepath = MHCQUANT_DATA_DIR / "multiqc_histogram_im.txt"
        data = _parse_histogram_file(filepath)

        assert data is not None, "histogram_im data should not be None for IM dataset"
        assert len(data) >= 1

    def test_histogram_scores_parsing(self):
        filepath = MHCQUANT_DATA_DIR / "multiqc_histogram_scores.txt"
        data = _parse_histogram_file(filepath)

        assert data is not None
        assert len(data) >= 1

    def test_length_dist_parsing(self):
        filepath = MHCQUANT_DATA_DIR / "multiqc_length_dist.txt"
        data = _parse_histogram_file(filepath)

        assert data is not None
        sample_name = list(data.keys())[0]
        # Length values should be peptide lengths (integers represented as float)
        assert all(v >= 0 for v in data[sample_name].values()), "Fractions should be non-negative"

    def test_chromatogram_parsing(self):
        filepath = MHCQUANT_DATA_DIR / "multiqc_chromatogram.txt"
        data = _parse_histogram_file(filepath)

        assert data is not None
        sample_name = list(data.keys())[0]
        # Chromatogram should have many RT bins
        assert len(data[sample_name]) > 10, "Chromatogram should have many RT bins"

    def test_mass_error_parsing(self):
        filepath = MHCQUANT_DATA_DIR / "multiqc_mass_error.txt"
        data = _parse_headerless_box_file(filepath)

        assert data is not None
        sample_name = list(data.keys())[0]
        assert isinstance(data[sample_name], list), "Box data should be a list of values"
        assert len(data[sample_name]) > 0, "Should have mass error values"

    def test_scores_xcorr_parsing(self):
        filepath = MHCQUANT_DATA_DIR / "multiqc_scores_xcorr.txt"
        data = _parse_headerless_box_file(filepath)

        assert data is not None
        sample_name = list(data.keys())[0]
        assert isinstance(data[sample_name], list)
        assert all(isinstance(v, float) for v in data[sample_name])

    def test_peptide_intensity_parsing(self):
        filepath = MHCQUANT_DATA_DIR / "multiqc_peptide_intensity.txt"
        data = _parse_headerless_box_file(filepath)

        assert data is not None
        sample_name = list(data.keys())[0]
        assert isinstance(data[sample_name], list)
        assert len(data[sample_name]) > 0

    def test_percolator_median_weights_parsing(self):
        filepath = MHCQUANT_DATA_DIR / "percolator_plot.txt"
        data, groups = _parse_percolator_median_weights(filepath)

        assert data is not None, "Percolator data should not be None"
        assert groups is not None, "Groups should not be None"
        assert len(groups) >= 2, "Should have at least 2 groups"
        assert len(data) > 0, "Should have at least one feature"

        # Each feature should map to at least one group with a weight value
        for feature, weights in data.items():
            assert isinstance(weights, dict), "Feature weights should be a dict"
            assert len(weights) >= 1, "Each feature should have at least one group weight"
            for group, val in weights.items():
                assert group in groups, f"Group {group} should be in groups list"
                assert isinstance(val, float), "Weight should be a float"

    def test_missing_file_returns_none(self):
        data = _parse_histogram_file("/nonexistent/file.txt")
        assert data is None

        data = _parse_headerless_box_file("/nonexistent/file.txt")
        assert data is None

        data, headers = _parse_general_stats("/nonexistent/file.txt")
        assert data is None

        data, groups = _parse_percolator_median_weights("/nonexistent/file.txt")
        assert data is None
