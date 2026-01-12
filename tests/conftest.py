"""
Configuration file for pytest.

This file contains fixtures and configuration for the pytest test suite.
"""

import gzip
import os
from pathlib import Path

import pandas as pd
import pytest

# Define the path to the test data
COMPRESSED_DATA_DIR = Path("tests/resources/maxquant/PXD003133-22min")
FRAGPIPE_DATA_DIR = Path("tests/resources/fragpipe")


@pytest.fixture
def maxquant_test_data_dir():
    """Return the path to the MaxQuant test data directory."""
    return COMPRESSED_DATA_DIR


@pytest.fixture
def proteingroups_file():
    """Return the path to the proteinGroups.txt.gz file."""
    return COMPRESSED_DATA_DIR / "proteinGroups.txt.gz"


@pytest.fixture
def evidence_file():
    """Return the path to the evidence.txt.gz file."""
    return COMPRESSED_DATA_DIR / "evidence.txt.gz"


@pytest.fixture
def summary_file():
    """Return the path to the summary.txt.gz file."""
    return COMPRESSED_DATA_DIR / "summary.txt.gz"


@pytest.fixture
def msms_file():
    """Return the path to the msms.txt.gz file."""
    return COMPRESSED_DATA_DIR / "msms.txt.gz"


@pytest.fixture
def msmsscans_file():
    """Return the path to the msmsScans.txt.gz file."""
    return COMPRESSED_DATA_DIR / "msmsScans.txt.gz"


@pytest.fixture
def parameters_file():
    """Return the path to the parameters.txt.gz file."""
    return COMPRESSED_DATA_DIR / "parameters.txt.gz"


@pytest.fixture
def proteingroups_df():
    """Return a pandas DataFrame of the proteinGroups.txt.gz file."""
    file_path = COMPRESSED_DATA_DIR / "proteinGroups.txt.gz"
    if os.path.exists(file_path):
        with gzip.open(file_path, "rt") as f:
            return pd.read_csv(f, sep="\t", low_memory=False)
    return None


@pytest.fixture
def evidence_df():
    """Return a pandas DataFrame of the evidence.txt.gz file."""
    file_path = COMPRESSED_DATA_DIR / "evidence.txt.gz"
    if os.path.exists(file_path):
        with gzip.open(file_path, "rt") as f:
            return pd.read_csv(f, sep="\t", low_memory=False)
    return None


@pytest.fixture
def summary_df():
    """Return a pandas DataFrame of the summary.txt.gz file."""
    file_path = COMPRESSED_DATA_DIR / "summary.txt.gz"
    if os.path.exists(file_path):
        with gzip.open(file_path, "rt") as f:
            return pd.read_csv(f, sep="\t", low_memory=False)
    return None


# FragPipe fixtures
@pytest.fixture
def fragpipe_test_data_dir():
    """Return the path to the FragPipe test data directory."""
    return FRAGPIPE_DATA_DIR


@pytest.fixture
def fragpipe_ion_file():
    """Return the path to the ion.tsv file."""
    return FRAGPIPE_DATA_DIR / "ion.tsv"


@pytest.fixture
def fragpipe_psm_file():
    """Return the path to the psm.tsv file."""
    return FRAGPIPE_DATA_DIR / "psm.tsv"


@pytest.fixture
def fragpipe_ion_df():
    """Return a pandas DataFrame of the ion.tsv file."""
    file_path = FRAGPIPE_DATA_DIR / "ion.tsv"
    if os.path.exists(file_path):
        return pd.read_csv(file_path, sep="\t", low_memory=False)
    return None


@pytest.fixture
def fragpipe_psm_df():
    """Return a pandas DataFrame of the psm.tsv file."""
    file_path = FRAGPIPE_DATA_DIR / "psm.tsv"
    if os.path.exists(file_path):
        return pd.read_csv(file_path, sep="\t", low_memory=False)
    return None