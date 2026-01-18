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
    ],
    "combined_protein": [
        "Protein", "Protein ID", "Entry Name", "Gene", "Protein Length",
        "Combined Total Peptides", "Combined Spectral Count", "Combined Unique Spectral Count",
        "Combined Total Spectral Count"
    ],
    "combined_peptide": [
        "Peptide", "Peptide Length", "Charges", "Protein", "Protein Start", "Protein End",
        "Combined Spectral Count"
    ],
    "combined_ion": [
        "Peptide Sequence", "Modified Sequence", "Charge", "Protein", "Gene", "Assigned Modifications"
    ]
}

REQUIRED_KEYWORDS = {
    "combined_ion": {
        "Spectral Count": False,
        "Match Type": False,
        "Intensity": False
    },
    "combined_peptide": {
        "Match Type": False,
        "Intensity": False
    },
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
    # fragpipe.workflow (FragPipe configuration parameters)
    # fragpipe-files.fp-manifest (experiment design/manifest)
    # diann-output files (see DIA-NN documentation)

    # Define all file types to look for
    file_types = [
        "psm", "ion", "combined_protein", "combined_peptide", "combined_ion",
        "workflow", "manifest", "fragger_params"
    ]
    fragpipe_files = {ft: [] for ft in file_types}

    # FragPipe *tsv Data
    for file_info in find_log_files("pmultiqc/tsv", filecontents=False):
        filename = file_info["fn"]
        full_path = os.path.join(file_info["root"], filename)

        # Match exact file names to avoid conflicts
        # e.g., "ion.tsv" should not match "combined_ion.tsv"
        if filename == "ion.tsv":
            fragpipe_files["ion"].append(full_path)
        elif filename == "combined_ion.tsv":
            fragpipe_files["combined_ion"].append(full_path)
        elif filename == "combined_protein.tsv":
            fragpipe_files["combined_protein"].append(full_path)
        elif filename == "combined_peptide.tsv":
            fragpipe_files["combined_peptide"].append(full_path)
        elif "psm" in filename:
            fragpipe_files["psm"].append(full_path)

    # FragPipe workflow file (configuration parameters)
    for file_info in find_log_files("pmultiqc/workflow", filecontents=False):
        filename = file_info["fn"]
        full_path = os.path.join(file_info["root"], filename)
        if filename == "fragpipe.workflow":
            fragpipe_files["workflow"].append(full_path)

    # FragPipe manifest file (experiment design)
    for file_info in find_log_files("pmultiqc/fp-manifest", filecontents=False):
        filename = file_info["fn"]
        full_path = os.path.join(file_info["root"], filename)
        fragpipe_files["manifest"].append(full_path)

    # MSFragger params file (search engine parameters)
    for file_info in find_log_files("pmultiqc/fragger_params", filecontents=False):
        filename = file_info["fn"]
        full_path = os.path.join(file_info["root"], filename)
        if filename == "fragger.params":
            fragpipe_files["fragger_params"].append(full_path)

    if any(fragpipe_files.values()):
        for k, v in fragpipe_files.items():
            if v:
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

    try:
        anchor_idx = ion_df.columns.get_loc("Mapped Proteins")
    except KeyError:
        return []

    sample_intensity_cols = list(ion_df.columns[anchor_idx + 1:])

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


def workflow_reader(file_path: str):
    """
    Read fragpipe.workflow file containing FragPipe configuration parameters.

    The workflow file is a key=value format file containing all parameters
    used in the FragPipe analysis.

    Parameters
    ----------
    file_path : str
        Path to the fragpipe.workflow file.

    Returns
    -------
    dict
        Dictionary containing parameter name -> value pairs.
    """
    parameters = {}

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if '=' in line:
                    key, value = line.split('=', 1)
                    parameters[key.strip()] = value.strip()
    except Exception as e:
        log.warning(f"Error reading workflow file {file_path}: {e}")
        return None

    log.info(f"Loaded {len(parameters)} parameters from workflow file")
    return parameters


def fragger_params_reader(file_path: str):
    """
    Read fragger.params file containing MSFragger search engine parameters.

    The fragger.params file is a key=value format file with MSFragger-specific
    search parameters including mass tolerances, enzyme settings, and modifications.

    Parameters
    ----------
    file_path : str
        Path to the fragger.params file.

    Returns
    -------
    dict
        Dictionary containing parameter name -> value pairs.
    """
    parameters = {}

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                # fragger.params uses "key = value" format with spaces
                if '=' in line:
                    key, value = line.split('=', 1)
                    parameters[key.strip()] = value.strip()
    except Exception as e:
        log.warning(f"Error reading fragger.params file {file_path}: {e}")
        return None

    log.info(f"Loaded {len(parameters)} parameters from fragger.params file")
    return parameters


def get_workflow_parameters_table(parameters: dict, fragger_params: dict = None):
    """
    Convert workflow and fragger parameters to table format for display.

    Extracts key parameters relevant for QC reporting similar to MaxQuant.
    Parameters from fragger.params can supplement or provide fallback for
    workflow parameters.

    Parameters
    ----------
    parameters : dict
        Dictionary of workflow parameters (from fragpipe.workflow).
    fragger_params : dict, optional
        Dictionary of MSFragger parameters (from fragger.params).

    Returns
    -------
    dict
        Dictionary formatted for table display with parameter/value structure.
    """
    if not parameters and not fragger_params:
        return None

    # Merge parameters - workflow takes precedence, fragger provides fallback
    merged_params = {}
    if fragger_params:
        merged_params.update(fragger_params)
    if parameters:
        merged_params.update(parameters)

    # Key parameters to display (similar to MaxQuant parameters table)
    # Format: (key_from_workflow, key_from_fragger, display_name)
    key_params = [
        # FragPipe version
        ("fragpipe.version", None, "FragPipe Version"),
        # Search engine settings - workflow keys and fragger.params equivalents
        ("msfragger.search_enzyme_name_1", "search_enzyme_name_1", "Enzyme"),
        ("msfragger.search_enzyme_cut_1", "search_enzyme_cut_1", "Enzyme Cut Site"),
        ("msfragger.allowed_missed_cleavage_1", "allowed_missed_cleavage_1", "Max Missed Cleavages"),
        ("msfragger.precursor_mass_lower", "precursor_mass_lower", "Precursor Mass Tolerance (Lower)"),
        ("msfragger.precursor_mass_upper", "precursor_mass_upper", "Precursor Mass Tolerance (Upper)"),
        ("msfragger.precursor_mass_units", "precursor_mass_units", "Precursor Mass Units"),
        ("msfragger.fragment_mass_tolerance", "fragment_mass_tolerance", "Fragment Mass Tolerance"),
        ("msfragger.fragment_mass_units", "fragment_mass_units", "Fragment Mass Units"),
        # Modifications - from workflow or fragger.params
        ("msfragger.variable_mod_01", "variable_mod_01", "Variable Modification 1"),
        ("msfragger.variable_mod_02", "variable_mod_02", "Variable Modification 2"),
        ("msfragger.variable_mod_03", "variable_mod_03", "Variable Modification 3"),
        # Database
        ("database.db-path", "database_name", "Database Path"),
        # IonQuant settings (workflow only)
        ("ionquant.mbr", None, "Match Between Runs (MBR)"),
        ("ionquant.normalization", None, "Normalization"),
        ("ionquant.requantify", None, "Requantify"),
        # TMT settings (workflow only)
        ("tmtintegrator.channel_num", None, "TMT Channels"),
        ("tmtintegrator.ref_tag", None, "TMT Reference Tag"),
        # FDR
        ("philosopher.filter--prot", None, "Protein FDR"),
        ("philosopher.filter--pep", None, "Peptide FDR"),
        ("philosopher.filter--psm", None, "PSM FDR"),
        # Additional fragger.params specific settings
        (None, "num_threads", "Number of Threads"),
        (None, "decoy_prefix", "Decoy Prefix"),
        (None, "isotope_error", "Isotope Error"),
        (None, "mass_offsets", "Mass Offsets"),
        (None, "precursor_true_tolerance", "Precursor True Tolerance"),
        (None, "precursor_true_units", "Precursor True Units"),
        (None, "calibrate_mass", "Calibrate Mass"),
        (None, "clip_nTerm_M", "Clip N-term Met"),
        (None, "digest_min_length", "Min Peptide Length"),
        (None, "digest_max_length", "Max Peptide Length"),
    ]

    table_data = {}
    row_num = 1

    for workflow_key, fragger_key, display_name in key_params:
        value = None

        # Try workflow key first
        if workflow_key and workflow_key in merged_params:
            value = merged_params[workflow_key]
        # Fallback to fragger key
        elif fragger_key and fragger_key in merged_params:
            value = merged_params[fragger_key]

        # Clean up value for display
        if value and value != "null" and str(value).strip():
            value = str(value).strip()
            # Extract filename from paths
            if display_name == "Database Path" and ("/" in value or "\\" in value):
                value = os.path.basename(value.replace("\\", "/"))
            # Skip empty modification slots
            if "Modification" in display_name and (not value or value == "0.0000 X 0"):
                continue
            table_data[row_num] = {
                "parameter": display_name,
                "value": value
            }
            row_num += 1

    if not table_data:
        return None

    return table_data


def manifest_reader(file_path: str):
    """
    Read FragPipe manifest file containing experiment design information.

    The manifest file (fp-manifest) contains file paths and experimental
    design information like sample names, groups, and data types.

    Parameters
    ----------
    file_path : str
        Path to the manifest file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing experiment design information.
    """
    try:
        # Manifest file is tab-separated with columns:
        # file_path, experiment, bioreplicate, data_type (optional)
        manifest_df = pd.read_csv(file_path, sep='\t', header=None)

        # Assign column names based on number of columns
        if len(manifest_df.columns) >= 4:
            manifest_df.columns = ['file_path', 'experiment', 'bioreplicate', 'data_type'] + \
                                  [f'col_{i}' for i in range(4, len(manifest_df.columns))]
        elif len(manifest_df.columns) == 3:
            manifest_df.columns = ['file_path', 'experiment', 'bioreplicate']
        elif len(manifest_df.columns) == 2:
            manifest_df.columns = ['file_path', 'experiment']
        else:
            manifest_df.columns = ['file_path']

        # Extract filename from path for display
        if 'file_path' in manifest_df.columns:
            manifest_df['file_name'] = manifest_df['file_path'].apply(
                lambda x: os.path.basename(str(x).replace("\\", "/"))
            )

        log.info(f"Loaded manifest with {len(manifest_df)} entries")
        return manifest_df

    except Exception as e:
        log.warning(f"Error reading manifest file {file_path}: {e}")
        return None


def get_experiment_design_table(manifest_df: pd.DataFrame):
    """
    Convert manifest DataFrame to experiment design table format.

    Parameters
    ----------
    manifest_df : pd.DataFrame
        DataFrame from manifest_reader.

    Returns
    -------
    dict
        Dictionary formatted for experiment design table display.
    """
    if manifest_df is None or manifest_df.empty:
        return None

    table_data = {}

    for idx, row in manifest_df.iterrows():
        file_name = row.get('file_name', row.get('file_path', f'File_{idx}'))
        entry = {
            "file_name": file_name,
        }
        if 'experiment' in row:
            entry["experiment"] = row['experiment']

        if 'bioreplicate' in row and pd.notna(row['bioreplicate']) and row['bioreplicate'] != '':
            entry["bioreplicate"] = row['bioreplicate']
        else:
            entry["bioreplicate"] = "-"

        if 'data_type' in row:
            entry["data_type"] = row['data_type']

        table_data[idx + 1] = entry

    return table_data


def combined_protein_reader(file_path: str):
    """
    Read combined_protein.tsv file from FragPipe output.

    The combined_protein.tsv file contains protein-level quantification
    data including MBR (Match Between Runs) information and per-sample
    intensity columns.

    Parameters
    ----------
    file_path : str
        Path to the combined_protein.tsv file.

    Returns
    -------
    tuple
        (protein_df, sample_intensity_cols, mbr_cols) where:
        - protein_df: DataFrame with protein-level data
        - sample_intensity_cols: List of sample intensity column names
        - mbr_cols: Dict mapping sample to MBR-related columns
    """
    try:
        protein_df = pd.read_csv(file_path, sep="\t", low_memory=False)
    except Exception as e:
        log.warning(f"Error reading combined_protein.tsv: {e}")
        return None, [], {}

    if protein_df.empty:
        log.warning("combined_protein.tsv is empty")
        return None, [], {}

    log.info(f"Loaded combined_protein.tsv with {len(protein_df)} proteins and {len(protein_df.columns)} columns")

    # Identify sample intensity columns
    # These typically follow patterns like "Sample MaxLFQ Intensity" or just sample names
    sample_intensity_cols = []

    # Look for MaxLFQ or Intensity columns per sample
    for col in protein_df.columns:
        col_lower = col.lower()
        if 'maxlfq' in col_lower or (
            'intensity' in col_lower and 'combined' not in col_lower
        ):
            # Skip metadata columns
            if not any(skip in col_lower for skip in ['total', 'spectral', 'razor']):
                sample_intensity_cols.append(col)

    # Look for MBR-related columns
    # FragPipe uses columns like "Sample Spectral Count" with MBR info
    # spectral_cols = [col for col in protein_df.columns if 'spectral count' in col.lower()]

    # Identify unique vs total spectral counts which can indicate MBR contribution
    # for col in spectral_cols:
    #     sample_base = col.replace(' Spectral Count', '').replace(' Total Spectral Count', '').strip()
    #     if sample_base not in mbr_cols:
    #         mbr_cols[sample_base] = {}
    #     if 'unique' in col.lower():
    #         mbr_cols[sample_base]['unique'] = col
    #     elif 'total' in col.lower():
    #         mbr_cols[sample_base]['total'] = col
    #     else:
    #         mbr_cols[sample_base]['spectral'] = col

    log.info(f"Found {len(sample_intensity_cols)} sample intensity columns")

    return protein_df, sample_intensity_cols


def get_protein_intensity_distribution(protein_df, sample_cols, contam_affix="CONT"):
    """
    Extract protein intensity distribution data for QC plots.

    Parameters
    ----------
    protein_df : pd.DataFrame
        DataFrame from combined_protein_reader.
    sample_cols : list
        List of sample intensity column names.
    contam_affix : str
        Contaminant identifier prefix.

    Returns
    -------
    tuple
        (sample_distribution, contaminant_distribution) - log2 transformed intensity distributions.
    """
    if protein_df is None or not sample_cols:
        return None, None

    sample_distribution = {}
    contaminant_distribution = {}

    # Identify protein column
    protein_col = None
    for col in ['Protein', 'Protein ID', 'Protein Group']:
        if col in protein_df.columns:
            protein_col = col
            break

    # Separate contaminants
    if protein_col:
        is_contaminant = protein_df[protein_col].str.contains(contam_affix, na=False, case=False)
        sample_df = protein_df[~is_contaminant]
        contam_df = protein_df[is_contaminant]
    else:
        sample_df = protein_df
        contam_df = pd.DataFrame()

    for col in sample_cols:
        if col not in sample_df.columns:
            continue

        if not pd.to_numeric(sample_df[col], errors='coerce').notna().all():
            continue

        # Sample intensities
        intensities = sample_df[col].dropna()
        valid_intensities = intensities[intensities > 0]
        if len(valid_intensities) > 0:
            sample_distribution[col] = np.log2(valid_intensities).tolist()

        # Contaminant intensities
        if not contam_df.empty and col in contam_df.columns:
            cont_intensities = contam_df[col].dropna()
            valid_cont = cont_intensities[cont_intensities > 0]
            if len(valid_cont) > 0:
                contaminant_distribution[col] = np.log2(valid_cont).tolist()

    return sample_distribution, contaminant_distribution


def combined_peptide_reader(file_path: str):
    """
    Read combined_peptide.tsv file from FragPipe output.

    The combined_peptide.tsv file contains peptide-level quantification
    data including MBR information.

    Parameters
    ----------
    file_path : str
        Path to the combined_peptide.tsv file.

    Returns
    -------
    tuple
        (peptide_df, sample_cols, mbr_info) where:
        - peptide_df: DataFrame with peptide-level data
        - sample_cols: List of sample column names
        - mbr_info: Dict with MBR-related statistics
    """
    try:
        peptide_df = pd.read_csv(file_path, sep="\t", low_memory=False)
    except Exception as e:
        log.warning(f"Error reading combined_peptide.tsv: {e}")
        return None, False

    if peptide_df.empty:
        log.warning("combined_peptide.tsv is empty")
        return None, False

    log.info(f"Loaded combined_peptide.tsv with {len(peptide_df)} peptides")

    if validate_columns_existence(
        df_columns=peptide_df.columns,
        data_name="combined_peptide"
    ):
        validate_columns = True
    else:
        validate_columns = False

    return peptide_df, validate_columns


def get_mbr_stats(protein_df, peptide_df, sample_cols):
    """
    Calculate Match-Between-Runs (MBR) statistics.

    Analyzes identification types to determine MBR contribution
    to protein and peptide identification counts.

    Parameters
    ----------
    protein_df : pd.DataFrame
        DataFrame from combined_protein_reader.
    peptide_df : pd.DataFrame
        DataFrame from combined_peptide_reader.
    sample_cols : list
        List of sample column names.

    Returns
    -------
    dict
        Dictionary containing MBR statistics per sample with:
        - 'proteins': dict with msms_only, mbr_only, both counts
        - 'peptides': dict with msms_only, mbr_only, both counts
    """
    mbr_stats = {}

    if protein_df is None and peptide_df is None:
        return mbr_stats

    # For FragPipe, we look at spectral counts to infer MBR contribution
    # Proteins/peptides with spectral count > 0 are MS/MS identified
    # Proteins/peptides that are quantified but have 0 spectral count may be MBR

    for sample in sample_cols:
        mbr_stats[sample] = {
            'proteins': {'msms_only': 0, 'mbr_only': 0, 'both': 0},
            'peptides': {'msms_only': 0, 'mbr_only': 0, 'both': 0}
        }

        # Protein-level MBR stats
        if protein_df is not None:
            spectral_col = None
            intensity_col = sample

            # Find corresponding spectral count column
            for col in protein_df.columns:
                if sample.replace(' MaxLFQ Intensity', '').replace(' Intensity', '') in col:
                    if 'spectral count' in col.lower():
                        spectral_col = col
                        break

            if spectral_col and intensity_col in protein_df.columns:
                has_spectral = protein_df[spectral_col] > 0
                has_intensity = protein_df[intensity_col] > 0

                msms_only = ((has_spectral) & (~has_intensity)).sum()
                mbr_only = ((~has_spectral) & (has_intensity)).sum()
                both = ((has_spectral) & (has_intensity)).sum()

                mbr_stats[sample]['proteins'] = {
                    'msms_only': int(msms_only),
                    'mbr_only': int(mbr_only),
                    'both': int(both)
                }

    return mbr_stats


def combined_ion_reader(file_path: str):
    """
    Read combined_ion.tsv file from FragPipe output.

    The combined_ion.tsv file contains ion-level quantification data
    across all samples, suitable for MS/MS counts per peak analysis.

    Parameters
    ----------
    file_path : str
        Path to the combined_ion.tsv file.

    Returns
    -------
    tuple
        (ion_df, sample_cols) where:
        - ion_df: DataFrame with ion-level data
        - sample_cols: List of sample column names
    """
    try:
        ion_df = pd.read_csv(file_path, sep="\t", low_memory=False)
    except Exception as e:
        log.warning(f"Error reading combined_ion.tsv: {e}")
        return None, False

    if ion_df.empty:
        log.warning("combined_ion.tsv is empty")
        return None, False

    log.info(f"Loaded combined_ion.tsv with {len(ion_df)} ions and {len(ion_df.columns)} columns")

    if validate_columns_existence(
        df_columns=ion_df.columns,
        data_name="combined_ion"
    ):
        validate_columns = True
    else:
        validate_columns = False

    return ion_df, validate_columns


def validate_columns_existence(df_columns, data_name: str):

    col_list = list(df_columns)

    required_keywords = REQUIRED_KEYWORDS[data_name]

    for col in col_list:
        if "Spectral Count" in col:
            required_keywords["Spectral Count"] = True
        if "Match Type" in col:
            required_keywords["Match Type"] = True
        if "Intensity" in col:
            if "MaxLFQ" not in col:
                required_keywords["Intensity"] = True

    all_passed = all(required_keywords.values())

    print(f"Check whether the data {data_name} meets the extraction requirements.")


    for key, found in required_keywords.items():
        status = "exists" if found else "is missing"
        print(f"{key}: {status}")
        
    return all_passed


def get_msms_counts_per_peak(ion_df):
    """
    Calculate MS/MS counts per peak statistics.

    Analyzes how many MS/MS spectra support each ion/peak identification.

    Parameters
    ----------
    ion_df : pd.DataFrame
        DataFrame from combined_ion_reader.
    Returns
    -------
    dict
        Dictionary containing MS/MS count statistics per sample.
    """
    if ion_df is None:
        return {}

    df = ion_df.copy()

    samples = [col.replace(' Match Type', '') for col in df.columns if ' Match Type' in col]

    plot_data = []

    for s in samples:
        spec_col = f"{s} Spectral Count"
        match_col = f"{s} Match Type"
        int_col = f"{s} Intensity"

        sample_df = df[df[int_col] > 0].copy()
        
        msms_dist = sample_df[sample_df[match_col] == 'MS/MS'][spec_col].value_counts().to_dict()

        for count, freq in msms_dist.items():
            plot_data.append({'run': s, 'ms/ms_count': int(count), 'peptide_count': freq})

    res_df = pd.DataFrame(plot_data)
    res_df["ms/ms_count"] = res_df["ms/ms_count"].apply(
        lambda x: ">=3" if x >= 3 else x
    )
    res_df = res_df.groupby(['run', 'ms/ms_count'])['peptide_count'].sum().reset_index()

    res_df["ms/ms_count"] = res_df["ms/ms_count"].astype(str)

    plot_dict = {}
    for raw_file, group in res_df.groupby("run"):
        group["freq"] = group["peptide_count"] / group["peptide_count"].sum() * 100
        plot_dict[raw_file] = dict(zip(group["ms/ms_count"], group["freq"]))

    oversampling = {
        "plot_data": plot_dict,
        "cats": list(res_df["ms/ms_count"].unique())
    }

    return oversampling


def cal_peptide_id_gain(df):
    df = df.copy()

    samples = [col.replace(' Match Type', '') for col in df.columns if ' Match Type' in col]

    peptide_counts = []

    for s in samples:
        match_col = f"{s} Match Type"
        int_col = f"{s} Intensity"

        sample_df = df[df[int_col] > 0].copy()

        match_type_count = sample_df[match_col].value_counts().to_dict()

        ms_count = 0
        mbr_count = 0
        for match_type, count in match_type_count.items():
            if match_type == "MS/MS":
                ms_count = count
            elif match_type == "MBR":
                mbr_count = count
        peptide_counts.append({'run': s, "ms/ms_count": int(ms_count), 'mbr': int(mbr_count)})

    count_df = pd.DataFrame(peptide_counts)

    count_df["MBRgain"] = (count_df["mbr"] / count_df["ms/ms_count"]) * 100

    temp_df = count_df[['run', 'ms/ms_count', 'mbr']].rename(
        columns={'ms/ms_count': 'MS/MS', 'mbr': 'MBR'}
    )
    plot_data = temp_df.set_index('run').to_dict(orient='index')

    mbr_gain = round(count_df["MBRgain"].mean(), 2)
    title_value = f"MBR gain: +{mbr_gain}%" if mbr_gain is not None else ""

    return {
        "plot_data": plot_data,
        "cats": ["MS/MS", "MBR"],
        "title_value": title_value
    }

