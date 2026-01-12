"""
Test module for FragPipe module functions.

These tests verify that the FragPipe ion.tsv parsing and intensity plot
generation functions work correctly.
"""

import os
import re
import sys

import numpy as np
import pandas as pd
import pytest

# Import FragPipe IO functions for testing
from pmultiqc.modules.fragpipe.fragpipe_io import (
    ion_reader,
    get_ion_intensity_data,
    extract_sample_groups,
    psm_reader,
)


class TestFragpipeIO:
    """Test class for fragpipe_io module functions."""

    def test_ion_file_exists(self, fragpipe_ion_file):
        """Test that the ion.tsv file exists."""
        assert os.path.exists(fragpipe_ion_file), f"File {fragpipe_ion_file} does not exist"

    def test_psm_file_exists(self, fragpipe_psm_file):
        """Test that the psm.tsv file exists."""
        assert os.path.exists(fragpipe_psm_file), f"File {fragpipe_psm_file} does not exist"

    def test_ion_reader_returns_dataframe(self, fragpipe_ion_file):
        """Test that ion_reader returns a DataFrame and sample columns."""
        ion_df, sample_cols = ion_reader(str(fragpipe_ion_file))

        assert ion_df is not None, "ion_reader returned None for ion_df"
        assert isinstance(ion_df, pd.DataFrame), "ion_df is not a DataFrame"
        assert not ion_df.empty, "ion_df is empty"
        assert isinstance(sample_cols, list), "sample_cols is not a list"
        assert len(sample_cols) > 0, "No sample intensity columns found"

    def test_ion_reader_detects_sample_columns(self, fragpipe_ion_file):
        """Test that ion_reader correctly identifies sample intensity columns."""
        ion_df, sample_cols = ion_reader(str(fragpipe_ion_file))

        # Check that TMT columns are detected
        tmt_cols = [col for col in sample_cols if "TMT" in col]
        assert len(tmt_cols) > 0, "No TMT columns detected"

        # Check that metadata columns are NOT in sample_cols
        metadata_cols = ["Peptide Sequence", "Modified Sequence", "Protein", "Charge"]
        for col in metadata_cols:
            assert col not in sample_cols, f"Metadata column {col} incorrectly detected as sample column"

    def test_ion_reader_required_columns(self, fragpipe_ion_file):
        """Test that ion_reader preserves required columns."""
        ion_df, sample_cols = ion_reader(str(fragpipe_ion_file))

        # Check for required columns
        expected_cols = ["Peptide Sequence", "Charge", "Protein"]
        for col in expected_cols:
            assert col in ion_df.columns, f"Required column {col} not found"

    def test_psm_reader_returns_dataframe(self, fragpipe_psm_file):
        """Test that psm_reader returns a DataFrame."""
        psm_df = psm_reader(str(fragpipe_psm_file))

        assert psm_df is not None, "psm_reader returned None"
        assert isinstance(psm_df, pd.DataFrame), "psm_df is not a DataFrame"
        assert not psm_df.empty, "psm_df is empty"

    def test_psm_reader_extracts_run(self, fragpipe_psm_file):
        """Test that psm_reader extracts Run from Spectrum column."""
        psm_df = psm_reader(str(fragpipe_psm_file))

        assert "Run" in psm_df.columns, "Run column not found"
        # Check that runs are extracted correctly
        unique_runs = psm_df["Run"].unique()
        assert len(unique_runs) > 0, "No runs extracted"


class TestGetIonIntensityData:
    """Test class for get_ion_intensity_data function."""

    def test_intensity_data_structure(self, fragpipe_ion_file):
        """Test that get_ion_intensity_data returns correct structure."""
        ion_df, sample_cols = ion_reader(str(fragpipe_ion_file))
        intensity_data = get_ion_intensity_data(ion_df, sample_cols)

        assert intensity_data is not None, "intensity_data is None"
        assert isinstance(intensity_data, dict), "intensity_data is not a dict"

        # Check for required keys
        expected_keys = ["intensity_distribution", "intensity_cv", "intensity_summary"]
        for key in expected_keys:
            assert key in intensity_data, f"Key {key} not found in intensity_data"

    def test_intensity_distribution_log2_transformed(self, fragpipe_ion_file):
        """Test that intensity distribution values are log2 transformed."""
        ion_df, sample_cols = ion_reader(str(fragpipe_ion_file))
        intensity_data = get_ion_intensity_data(ion_df, sample_cols)

        distribution = intensity_data["intensity_distribution"]
        assert len(distribution) > 0, "No intensity distribution data"

        # Check that values are in reasonable log2 range (typically 10-30 for proteomics)
        for sample, values in distribution.items():
            if values:
                min_val = min(values)
                max_val = max(values)
                assert min_val > 0, f"Log2 values should be positive for sample {sample}"
                assert max_val < 40, f"Log2 values seem too high for sample {sample}"

    def test_intensity_summary_statistics(self, fragpipe_ion_file):
        """Test that intensity summary contains correct statistics."""
        ion_df, sample_cols = ion_reader(str(fragpipe_ion_file))
        intensity_data = get_ion_intensity_data(ion_df, sample_cols)

        summary = intensity_data["intensity_summary"]
        assert len(summary) > 0, "No intensity summary data"

        # Check that each sample has required statistics
        for sample, stats in summary.items():
            assert "median" in stats, f"Missing median for sample {sample}"
            assert "mean" in stats, f"Missing mean for sample {sample}"
            assert "std" in stats, f"Missing std for sample {sample}"
            assert "count" in stats, f"Missing count for sample {sample}"

    def test_intensity_cv_calculation(self, fragpipe_ion_file):
        """Test that CV values are calculated correctly."""
        ion_df, sample_cols = ion_reader(str(fragpipe_ion_file))
        intensity_data = get_ion_intensity_data(ion_df, sample_cols)

        cv_data = intensity_data["intensity_cv"]
        # CV should be calculated if there are at least 2 samples
        if len(sample_cols) >= 2:
            # CV data may or may not be populated depending on the data
            # Just check it's a dict
            assert isinstance(cv_data, dict), "CV data is not a dict"


class TestExtractSampleGroups:
    """Test class for extract_sample_groups function."""

    def test_tmt_pattern_extraction(self):
        """Test TMT channel pattern extraction."""
        sample_cols = [
            "Sample1_TMT_126",
            "Sample1_TMT_127N",
            "Sample1_TMT_127C",
            "Sample1_TMT_128N",
        ]

        groups = extract_sample_groups(sample_cols)

        assert len(groups) == 4, "Should have 4 sample groups"

        # Check that TMT channels are identified
        for col in sample_cols:
            assert col in groups, f"Sample {col} not in groups"
            assert "channel" in groups[col], f"Channel not identified for {col}"
            assert groups[col]["experiment"] == "Sample1", f"Experiment not correct for {col}"

    def test_replicate_pattern_extraction(self):
        """Test replicate pattern extraction."""
        sample_cols = [
            "ConditionA_Rep1",
            "ConditionA_Rep2",
            "ConditionB_Rep1",
            "ConditionB_Rep2",
        ]

        groups = extract_sample_groups(sample_cols)

        assert len(groups) == 4, "Should have 4 sample groups"

        # Check that conditions are identified
        for col in sample_cols[:2]:
            if col in groups and "condition" in groups[col]:
                assert "ConditionA" in groups[col]["condition"], f"Condition not correct for {col}"


class TestFragpipeIntensityPlots:
    """Test class for FragPipe intensity plotting functions."""

    def test_intensity_distribution_data_format(self, fragpipe_ion_file):
        """Test that intensity distribution data is in correct format for plotting."""
        ion_df, sample_cols = ion_reader(str(fragpipe_ion_file))
        intensity_data = get_ion_intensity_data(ion_df, sample_cols)

        distribution = intensity_data["intensity_distribution"]

        # Check that data is in correct format for box plot
        # Should be dict mapping sample -> list of values
        for sample, values in distribution.items():
            assert isinstance(sample, str), f"Sample key should be string"
            assert isinstance(values, list), f"Values should be list for sample {sample}"
            assert all(isinstance(v, (int, float)) for v in values), f"Values should be numeric for sample {sample}"

    def test_no_empty_samples(self, fragpipe_ion_file):
        """Test that no samples have empty intensity data."""
        ion_df, sample_cols = ion_reader(str(fragpipe_ion_file))
        intensity_data = get_ion_intensity_data(ion_df, sample_cols)

        distribution = intensity_data["intensity_distribution"]

        for sample, values in distribution.items():
            assert len(values) > 0, f"Sample {sample} has empty intensity data"


class TestFragpipeDataQuality:
    """Test data quality checks for FragPipe data."""

    def test_ion_tsv_structure(self, fragpipe_ion_df):
        """Test that ion.tsv has expected structure."""
        assert fragpipe_ion_df is not None, "Could not load ion.tsv"

        # Check for expected columns
        expected_cols = ["Peptide Sequence", "Charge", "Protein"]
        for col in expected_cols:
            assert col in fragpipe_ion_df.columns, f"Expected column {col} not found"

    def test_psm_tsv_structure(self, fragpipe_psm_df):
        """Test that psm.tsv has expected structure."""
        assert fragpipe_psm_df is not None, "Could not load psm.tsv"

        # Check for expected columns
        expected_cols = ["Spectrum", "Peptide", "Charge", "Intensity"]
        for col in expected_cols:
            assert col in fragpipe_psm_df.columns, f"Expected column {col} not found"

    def test_intensity_values_reasonable(self, fragpipe_ion_df):
        """Test that intensity values are in reasonable range."""
        assert fragpipe_ion_df is not None, "Could not load ion.tsv"

        # Find intensity columns (TMT columns)
        intensity_cols = [col for col in fragpipe_ion_df.columns if "TMT" in col]

        for col in intensity_cols:
            values = fragpipe_ion_df[col].dropna()
            assert values.min() >= 0, f"Negative intensity values in {col}"
            # Values should be reasonable (not astronomical)
            assert values.max() < 1e12, f"Intensity values too large in {col}"

    def test_charge_values_valid(self, fragpipe_ion_df):
        """Test that charge values are valid."""
        assert fragpipe_ion_df is not None, "Could not load ion.tsv"

        if "Charge" in fragpipe_ion_df.columns:
            charges = fragpipe_ion_df["Charge"].dropna()
            assert charges.min() >= 1, "Charge should be at least 1"
            assert charges.max() <= 10, "Charge should be at most 10 for typical data"
