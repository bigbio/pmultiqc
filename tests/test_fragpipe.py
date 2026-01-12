"""
Test module for FragPipe module functions.

These tests verify that the FragPipe ion.tsv parsing and intensity plot
generation functions work correctly.

Test data is from PRIDE dataset PXD066146:
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD066146/

The test data files (ion.tsv, psm.tsv) are subsets (first 1000 rows)
of the full dataset to keep test execution fast.
"""

import os

import numpy as np
import pandas as pd

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
        _, sample_cols = ion_reader(str(fragpipe_ion_file))

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
        assert "intensity_distribution" in intensity_data, "Key intensity_distribution not found"

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
            assert isinstance(sample, str), "Sample key should be string"
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


class TestPXD066146Data:
    """Test class specific to PXD066146 dataset characteristics."""

    def test_tmt_channels_detected(self, fragpipe_ion_file):
        """Test that TMT channels from PXD066146 are correctly detected."""
        ion_df, sample_cols = ion_reader(str(fragpipe_ion_file))

        # PXD066146 uses TMT labeling with 33075_TMT_* pattern
        tmt_cols = [col for col in sample_cols if "33075_TMT" in col]
        assert len(tmt_cols) >= 10, f"Expected at least 10 TMT channels, found {len(tmt_cols)}"

        # Check for specific TMT channels
        expected_channels = ["33075_TMT_126", "33075_TMT_127N", "33075_TMT_127C"]
        for channel in expected_channels:
            assert channel in sample_cols, f"Expected TMT channel {channel} not found"

    def test_ion_count_reasonable(self, fragpipe_ion_file):
        """Test that ion.tsv has reasonable number of ions."""
        ion_df, sample_cols = ion_reader(str(fragpipe_ion_file))

        # Test data has 1000 rows (subset of full dataset)
        assert len(ion_df) == 1000, f"Expected 1000 ions in test data, found {len(ion_df)}"

    def test_protein_identifiers_format(self, fragpipe_ion_df):
        """Test that protein identifiers follow UniProt format."""
        assert fragpipe_ion_df is not None, "Could not load ion.tsv"

        proteins = fragpipe_ion_df["Protein"].dropna().unique()

        # Check for UniProt accession pattern (sp|XXXXX|YYYY_SPECIES)
        uniprot_pattern_count = sum(1 for p in proteins if "sp|" in str(p) or "tr|" in str(p))
        assert uniprot_pattern_count > 0, "No UniProt-style protein identifiers found"

    def test_intensity_distribution_values(self, fragpipe_ion_file):
        """Test intensity distribution values for PXD066146 data."""
        ion_df, sample_cols = ion_reader(str(fragpipe_ion_file))
        intensity_data = get_ion_intensity_data(ion_df, sample_cols)

        distribution = intensity_data["intensity_distribution"]

        # Check that values are in typical proteomics range
        # log2(1000) ~ 10, log2(1e9) ~ 30
        for sample, values in distribution.items():
            if values:
                median = np.median(values)
                assert 8 < median < 25, f"Median intensity {median} outside typical range for {sample}"

    def test_modified_sequences_present(self, fragpipe_ion_df):
        """Test that modified sequences are properly formatted."""
        assert fragpipe_ion_df is not None, "Could not load ion.tsv"

        if "Modified Sequence" in fragpipe_ion_df.columns:
            mod_seqs = fragpipe_ion_df["Modified Sequence"].dropna()
            assert len(mod_seqs) > 0, "No modified sequences found"

            # Check for typical modification patterns (e.g., n[43] for N-term acetyl)
            has_modifications = any("[" in str(seq) for seq in mod_seqs.head(100))
            assert has_modifications, "No modification annotations found in sequences"

    def test_psm_spectrum_format(self, fragpipe_psm_df):
        """Test that PSM spectrum identifiers follow FragPipe format."""
        assert fragpipe_psm_df is not None, "Could not load psm.tsv"

        spectra = fragpipe_psm_df["Spectrum"].dropna()
        assert len(spectra) > 0, "No spectra found"

        # FragPipe spectrum format: filename.scan.scan.charge
        sample_spectrum = str(spectra.iloc[0])
        parts = sample_spectrum.split(".")
        assert len(parts) >= 3, f"Spectrum format unexpected: {sample_spectrum}"

    def test_run_extraction_from_spectrum(self, fragpipe_psm_file):
        """Test that run names are correctly extracted from spectrum IDs."""
        psm_df = psm_reader(str(fragpipe_psm_file))

        runs = psm_df["Run"].unique()
        assert len(runs) >= 1, "No runs extracted from spectra"

        # Check that run names don't contain scan numbers
        for run in runs:
            # Run should not be just numbers
            assert not str(run).isdigit(), f"Run name {run} appears to be just a scan number"
