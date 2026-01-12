import os
import re
import pandas as pd
import numpy as np

from pmultiqc.modules.common.logging import get_logger


# Initialise the module logger via centralized logger
log = get_logger("pmultiqc.modules.fragpipe.fragpipe_io")


REQUIRED_COLS = {
    "psm": [
        "Spectrum", "Peptide", "Modified Peptide", "Charge", "Retention", "Intensity",
        "Delta Mass", "Number of Missed Cleavages", "Is Unique", "Protein", "Hyperscore",
        "Assigned Modifications"
    ],
    "ion": [
        "Peptide Sequence", "Modified Sequence", "Charge", "Protein", "Intensity"
    ]
}


# FragPipe File Paths
def get_fragpipe_files(find_log_files):

    # FragPipe results (https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html)
    # Main report files
    # psm.tsv (from Philosopher, updated by PTM-Shepherd and IonQuant)
    # ion.tsv (from Philosopher, overwritten by IonQuant)
    # peptide.tsv (from Philosopher, overwritten by IonQuant)
    # protein.tsv (from Philosopher, overwritten by IonQuant)
    # combined_ion.tsv (from Philosopher, overwritten by IonQuant)
    # combined_modified_peptide.tsv (from IonQuant)
    # combined_peptide.tsv (from Philosopher, overwritten by IonQuant)
    # combined_protein.tsv (from Philosopher, overwritten by IonQuant)
    # diann-output files (see DIA-NN documentation)

    required_files = ["psm", "ion"]
    req_set = set(required_files)

    fragpipe_files = {req: [] for req in required_files}

    # FragPipe *tsv Data
    for file_info in find_log_files("pmultiqc/tsv", filecontents=False):
        filename = file_info["fn"]
        full_path = os.path.join(file_info["root"], filename)

        for req in req_set:
            # Match exact file names to avoid conflicts
            # e.g., "ion.tsv" should not match "combined_ion.tsv"
            if req == "ion" and filename == "ion.tsv":
                fragpipe_files[req].append(full_path)
            elif req == "psm" and "psm" in filename:
                fragpipe_files[req].append(full_path)

    if any(fragpipe_files.values()):

        for k, v in fragpipe_files.items():
            log.info(f"FragPipe data loaded: {k} ({len(v)} files).")
            log.debug(f"FragPipe data loaded: {k}: {v}")

        return fragpipe_files

    else:
        log.warning("FragPipe data loading failed: Please verify the input folder path.")

        return None


def psm_reader(file_path: str):

    psm_df = pd.read_csv(file_path, sep="\t")

    if "Spectrum" not in psm_df.columns:
        raise ValueError("psm.tsv must contain a 'Spectrum' column")

    required_cols = [c for c in REQUIRED_COLS["psm"] if c in psm_df.columns]

    psm_df["Run"] = psm_df["Spectrum"].astype(str).str.rsplit(".", n=3).str[0]
    required_cols.append("Run")

    return psm_df[required_cols].copy()


def ion_reader(file_path: str):
    """
    Read ion.tsv file from FragPipe/IonQuant output.

    The ion.tsv file contains ion-level quantification data with per-sample
    intensity columns (e.g., TMT channels or LFQ intensities).

    Parameters
    ----------
    file_path : str
        Path to the ion.tsv file.

    Returns
    -------
    tuple
        (ion_df, sample_intensity_cols) where:
        - ion_df: DataFrame with ion-level data including sample intensities
        - sample_intensity_cols: List of column names containing sample intensities
    """
    ion_df = pd.read_csv(file_path, sep="\t", low_memory=False)

    if "Peptide Sequence" not in ion_df.columns:
        raise ValueError("ion.tsv must contain a 'Peptide Sequence' column")

    # Get required columns that exist
    required_cols = [c for c in REQUIRED_COLS["ion"] if c in ion_df.columns]

    # Identify sample intensity columns
    # These are columns that are NOT in the metadata columns
    metadata_cols = {
        "Peptide Sequence", "Modified Sequence", "Prev AA", "Next AA",
        "Peptide Length", "M/Z", "Charge", "Observed Mass", "Probability",
        "Expectation", "Spectral Count", "Intensity", "Assigned Modifications",
        "Observed Modifications", "Protein", "Protein ID", "Entry Name",
        "Gene", "Protein Description", "Mapped Genes", "Mapped Proteins"
    }

    # Sample intensity columns are those not in metadata
    all_cols = set(ion_df.columns)
    sample_intensity_cols = sorted(list(all_cols - metadata_cols))

    # Filter to only numeric columns that look like sample intensities
    # These typically have numeric values and may contain patterns like TMT, LFQ, etc.
    valid_sample_cols = []
    for col in sample_intensity_cols:
        # Check if column contains numeric data
        if ion_df[col].dtype in [np.float64, np.int64, np.float32, np.int32]:
            valid_sample_cols.append(col)
        elif ion_df[col].dtype == object:
            # Try to convert to numeric
            try:
                pd.to_numeric(ion_df[col], errors='raise')
                valid_sample_cols.append(col)
            except (ValueError, TypeError):
                pass

    sample_intensity_cols = valid_sample_cols

    log.info(f"Found {len(sample_intensity_cols)} sample intensity columns in ion.tsv")
    if sample_intensity_cols:
        log.debug(f"Sample intensity columns: {sample_intensity_cols[:5]}...")

    # Select required columns plus sample intensity columns
    select_cols = required_cols + sample_intensity_cols
    select_cols = [c for c in select_cols if c in ion_df.columns]

    return ion_df[select_cols].copy(), sample_intensity_cols


def get_ion_intensity_data(ion_df, sample_cols):
    """
    Extract and process intensity data from ion.tsv for QC plots.

    Parameters
    ----------
    ion_df : pd.DataFrame
        DataFrame from ion_reader containing ion-level data.
    sample_cols : list
        List of column names containing sample intensities.

    Returns
    -------
    dict
        Dictionary containing:
        - 'intensity_distribution': dict mapping sample -> list of log2 intensities
    """
    if not sample_cols:
        log.warning("No sample intensity columns found in ion.tsv")
        return None

    result = {
        'intensity_distribution': {},
    }

    # Calculate intensity distribution (log2 transformed)
    for sample in sample_cols:
        intensities = ion_df[sample].dropna()
        # Filter out zero and negative values before log transformation
        valid_intensities = intensities[intensities > 0]
        if len(valid_intensities) > 0:
            log_intensities = np.log2(valid_intensities)
            result['intensity_distribution'][sample] = log_intensities.tolist()

    return result


def extract_sample_groups(sample_cols):
    """
    Extract sample groups from column names for grouping in plots.

    Tries to identify experimental conditions from sample column names.
    Common patterns include: TMT channels, replicates, conditions.

    Parameters
    ----------
    sample_cols : list
        List of sample column names.

    Returns
    -------
    dict
        Dictionary mapping sample name to group name.
    """
    sample_groups = {}

    # Try to extract sample groups from column names
    if len(sample_cols) >= 2:
        # TMT pattern: e.g., "33075_TMT_126", "33075_TMT_127N"
        tmt_pattern = re.compile(r"(.+?)_TMT_(\d+[NC]?)$")

        for col in sample_cols:
            match = tmt_pattern.match(col)
            if match:
                experiment = match.group(1)
                channel = match.group(2)
                sample_groups[col] = {
                    "experiment": experiment,
                    "channel": f"TMT_{channel}",
                    "group": experiment
                }
            else:
                # Try to extract replicate pattern: e.g., "Sample_Rep1", "Sample_Rep2"
                rep_pattern = re.compile(r"(.+?)_?[Rr]ep(\d+)$")
                match = rep_pattern.match(col)
                if match:
                    condition = match.group(1)
                    replicate = match.group(2)
                    sample_groups[col] = {
                        "condition": condition,
                        "replicate": replicate,
                        "group": condition
                    }
                else:
                    # Default: use column name as is
                    sample_groups[col] = {
                        "group": col
                    }

    return sample_groups

