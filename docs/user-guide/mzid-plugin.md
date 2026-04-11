# Using pmultiqc with mzIdentML Files

The `mzidentml` plugin processes [mzIdentML](https://www.psidev.info/mzidentml) (`.mzid`) files, the PSI standard format for peptide and protein identification results. It pairs mzIdentML files with their corresponding spectrum files (`.mzML` or `.mgf`) to produce combined identification and spectral quality metrics.

## Supported Input Files

| File Pattern | Description |
|---|---|
| `*.mzid` | mzIdentML identification results (required) |
| `*.mzML` | Spectrum files for MS-level QC (optional but recommended) |
| `*.mgf` | Mascot Generic Format spectrum files (alternative to mzML) |

At minimum, one or more `.mzid` files must be present. When spectrum files are also available in the same directory, MS1-level metrics (TIC, BPC, charge distributions) are included in the report.

## Running the Report

```bash
# Basic usage — directory with .mzid files
multiqc --mzid-plugin /path/to/mzid/results -o ./report

# With both mzid and mzML files in the same directory
multiqc --mzid-plugin /path/to/mzid/results -o ./report

# Disable hover tooltips
multiqc --mzid-plugin /path/to/mzid/results -o ./report --disable-hoverinfo

# Use a custom decoy prefix (default: DECOY_)
multiqc --mzid-plugin /path/to/mzid/results -o ./report --decoy-affix REV_ --affix-type prefix

# Disable protein/peptide tables for large datasets
multiqc --mzid-plugin /path/to/mzid/results -o ./report --disable-table
```

## QC Sections Generated

### Identification Summary Table

A per-file table with:
- Total spectra in the run
- Number of identified PSMs
- MS/MS identification rate (%)
- Total peptides identified
- Total proteins identified

### Charge State Distribution

Bar chart of precursor charge states extracted from mzIdentML `SpectrumIdentificationResult` entries. Parsed via `pmultiqc/modules/common/mzidentml_utils.py:get_mzidentml_charge()`.

### Identification Rate over Retention Time

Histogram of PSM identifications binned by retention time (seconds). Derived from `get_mzid_rt_id()`. Shows whether identifications are uniformly distributed across the gradient or concentrated at specific elution windows.

### Delta Mass (ppm and Da)

Mass accuracy distribution computed from the experimental vs. theoretical mass differences for all identified PSMs. A tight, zero-centered distribution indicates well-calibrated mass accuracy.

### Oversampling

Fraction of precursors that were selected for MS2 more than once per LC run. High oversampling indicates inefficient data acquisition (e.g., too many MS2 per cycle) and is a common indicator of suboptimal DDA settings.

### Number of Peptides per Protein

Distribution of unique peptide counts per protein. Proteins with a single peptide identification (singletons) are highlighted.

### Long-Trend Plots

Cumulative identification counts (PSMs, peptides, proteins) as a function of retention time. Useful for detecting sudden drops in identification rate that could indicate a clogged column or ESI instability.

### MS Spectral Metrics (when mzML or MGF is provided)

When spectrum files are co-located with mzIdentML files, the following additional sections are generated:

- **Total Ion Chromatogram (TIC)** — ion signal over retention time
- **Base Peak Chromatogram (BPC)** — most abundant ion at each time point
- **Precursor Charge Distribution** — charge state frequencies from survey scan headers
- **Peaks per MS2 Spectrum** — fragment ion count histogram
- **Peak Intensity Distribution** — MS2 fragment intensity histogram

### Quantification Table (if present)

If the mzIdentML file contains quantification data (via `DataCollection/AnalysisData` elements), a summary table of protein-level quantification values is generated via `draw_mzid_quant_table()`.

## File Pairing

The mzIdentML module (`pmultiqc/modules/mzidentml/mzidentml.py`) automatically pairs `.mzid` files with spectrum files based on matching file name prefixes. For example, `sample1.mzid` is paired with `sample1.mzML` or `sample1.mgf` if present.

Unpaired mzIdentML files produce identification statistics only, without MS-level chromatographic plots.

## Notes

- Parsing is implemented using the [pyteomics](https://pyteomics.readthedocs.io/) `mzid` and `mgf` readers.
- Very large mzIdentML files (>500 MB) may require additional memory. Consider using `--disable-table` to reduce rendering overhead.
- The mzIdentML standard does not enforce a specific search engine score; pmultiqc reads whatever PSM-level score is present in the `SpectrumIdentificationItem` elements.
