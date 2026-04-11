# Using pmultiqc with MaxQuant Output

The `maxquant` plugin parses the standard text output files produced by [MaxQuant](https://www.maxquant.org) and generates a QC report comparable to the quantms module. MaxQuant output directories typically contain a `txt/` subfolder with all relevant files.

## Supported Input Files

| File | Description |
|---|---|
| `parameters.txt` | Analysis parameters and software settings |
| `proteinGroups.txt` | Protein group identification and quantification |
| `summary.txt` | Per-raw-file MS/MS identification summary |
| `evidence.txt` | Peptide-level evidence with retention times and intensities |
| `msms.txt` | MS/MS scan results with fragment ion information |
| `msmsScans.txt` | Detailed MS/MS scan information |
| `*sdrf.tsv` | SDRF-Proteomics metadata file (optional) |

Place all files in the same directory (or in a `txt/` subfolder). pmultiqc will locate them automatically.

## Running the Report

```bash
# Basic usage
multiqc --maxquant-plugin /path/to/maxquant/results -o ./report

# Include SDRF for experimental design visualization
multiqc --maxquant-plugin /path/to/maxquant/results -o ./report

# Disable large protein/peptide tables for datasets with many samples
multiqc --maxquant-plugin /path/to/maxquant/results -o ./report --disable-table

# Customize the contaminant prefix (MaxQuant default is CON__)
multiqc --maxquant-plugin /path/to/maxquant/results -o ./report --contaminant-affix CON__

# Disable hover tooltips for static reporting
multiqc --maxquant-plugin /path/to/maxquant/results -o ./report --disable-hoverinfo
```

## QC Sections Generated

### Experimental Design

If an SDRF file is present alongside the MaxQuant output, pmultiqc renders an experimental design table showing sample-to-file mapping, biological replicates, fractions, and technical replicates.

### Parameters Summary

Key parameters extracted from `parameters.txt`:

- MaxQuant version
- Enzyme and missed cleavages setting
- Fixed and variable modifications
- Mass tolerance (MS1 and MS2)
- FDR thresholds

### Protein Groups Statistics

From `proteinGroups.txt`:

- **Protein Identification Summary** — total protein groups, unique proteins, and contaminant fraction
- **Potential Contaminants per Group** — bar chart showing the proportion of contaminant signal per sample group
- **Number of Peptides per Protein** — distribution of unique peptide counts per protein group

### MS/MS Identification Statistics

From `summary.txt`:

- **MS/MS Identification Rate** — per-raw-file identification percentage bar chart
- **Spectra Tracking** — MS1 scans, MS2 scans, and identified PSM counts per file

### Evidence-Level Metrics

From `evidence.txt`:

- **Charge State Distribution** — frequency of precursor charge states across all identified peptides
- **Retention Time Distribution** — histogram of peptide identifications over the LC gradient
- **Delta Mass** — mass error distribution in Da and ppm
- **Peptide Intensity Distribution** — log-transformed intensity histogram for LFQ or iBAQ values
- **Missed Cleavages** — proportion of peptides with 0, 1, or 2 missed tryptic cleavages
- **Peptide Length Distribution** — sequence length histogram

### MS/MS Scan Metrics

From `msms.txt` and `msmsScans.txt`:

- **Search Engine Scores** — Andromeda score distribution for identified PSMs
- **Peaks per MS2 Spectrum** — histogram of fragment ion peak counts per spectrum
- **Peak Intensity Distribution** — distribution of MS2 fragment ion intensities
- **Oversampling** — analysis of precursors selected multiple times per LC run

### Heatmap

Sample-by-sample correlation matrix computed from LFQ intensities in `proteinGroups.txt`. Useful for identifying outliers or unexpected batch structure.

## File Detection

The MaxQuant module (`pmultiqc/modules/maxquant/maxquant_io.py`) searches for files using the following precedence:

1. Direct file names in the analysis directory
2. Files in a `txt/` subdirectory
3. SDRF file ending in `sdrf.tsv` in the same location

If none of the required files are found, the module is skipped and no MaxQuant section appears in the report.

## Notes

- MaxQuant LFQ intensities are used for quantification comparisons; iBAQ values are also parsed when present.
- The `msmsScans.txt` file is optional but enables additional scan-level metrics when present.
- For TMT or SILAC experiments run in MaxQuant, the relevant intensity columns are detected automatically from the column headers.
