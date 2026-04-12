# Using pmultiqc with quantms Pipeline Output

The `quantms` plugin is the primary integration point between pmultiqc and the [quantms](https://github.com/bigbio/quantms) pipeline. It parses mzTab, mzML, idXML, and experimental design files produced by quantms **DDA workflows** (LFQ and TMT/iTRAQ).

!!! info "DIA experiments use a different plugin"
    For DIA experiments processed with **[quantmsdiann](https://github.com/bigbio/quantmsdiann)**, use the `--diann-plugin` instead. See the [DIA-NN Plugin Guide](diann-plugin.md). DIA output is `report.tsv`/`report.parquet` — not mzTab.

## Supported Input Files

The quantms module recognizes the following files from a quantms **DDA** results directory:

| File Pattern | Description |
|---|---|
| `experimental_design.tsv` | Experimental design linking raw files to samples |
| `*.mzTab` | DDA identification and quantification results |
| `*msstats*.csv` | MSstats or MSstatsTMT input tables |
| `*.mzML` | Spectrum files for MS1/MS2 signal inspection |
| `*ms_info.tsv` | Precomputed MS quality control information |
| `*.idXML` | OpenMS identification results (optional) |
| `*.yml` | Pipeline parameter file (optional) |

## Running the Report

```bash
# Basic usage — point to the quantms results directory
multiqc --quantms-plugin /path/to/quantms/results -o ./report

# Remove decoy hits before counting (default: enabled)
multiqc --quantms-plugin results/ -o report/ --remove-decoy

# Use a custom decoy prefix
multiqc --quantms-plugin results/ -o report/ --decoy-affix REV_ --affix-type prefix

# Assign conditions from an SDRF factor column
multiqc --quantms-plugin results/ -o report/ --condition "factor value[disease]"

# Use spectral counting instead of feature intensity for LFQ quantification
multiqc --quantms-plugin results/ -o report/ --quantification-method spectral_counting

# Skip idXML parsing for faster execution on large datasets
multiqc --quantms-plugin results/ -o report/ --ignored-idxml

# Disable protein/peptide tables on very large datasets
multiqc --quantms-plugin results/ -o report/ --disable-table
```

## QC Sections Generated

### Experimental Design

An interactive table mapping each raw file to its sample, biological replicate, fraction, and technical replicate as defined in the SDRF/experimental design file.

### Pipeline Performance Overview

A sparkline summary panel showing per-sample scores for:

- **Contaminant Score** — fraction of total signal attributable to common contaminants (e.g., keratins, trypsin)
- **Peptide Intensity** — distribution uniformity of peptide intensities across samples
- **Charge Score** — expected vs. observed charge state composition
- **Missed Cleavages** — proportion of peptides with 0, 1, or 2 missed cleavages
- **ID Rate over RT** — identification rate across the retention time gradient
- **MS2 Oversampling** — fraction of precursors selected for MS2 more than once in the same run
- **Peptide Missing Values** — fraction of missing quantification values per sample

### MS1 Information

Chromatographic and spectral quality metrics extracted from mzML files:

- **Total Ion Chromatograms (TIC)** — ion signal over retention time per sample
- **Base Peak Chromatograms (BPC)** — most abundant ion at each RT point
- **MS1 Peak Count** — number of survey scans per file

### Identification Statistics

- **Spectra Tracking** — per-file counts of MS2 spectra, identified PSMs, and identification rate
- **MS/MS Identification Rate** — bar chart of identification percentage by run and by sample
- **Search Engine Scores** — distribution of PSM-level identification scores
- **Delta Mass (Da / ppm)** — mass accuracy distribution for all identified PSMs
- **Precursor Charge Distribution** — frequency of charge states +1 through +6+
- **Peaks per MS2 Spectrum** — histogram of fragment ion counts per spectrum
- **Peak Intensity Distribution** — histogram of MS2 peak intensities

### Peptide and Protein Results

- **Identifications Summary Table** — protein and peptide counts per sample
- **Number of Peptides per Protein** — distribution of peptides assigned per protein group
- **Peptide Length Distribution** — frequency of identified peptide sequence lengths
- **Peptide Intensity Table** — quantification values across conditions (first 500 peptides)
- **PSM Table** — per-PSM metrics for the first 500 entries
- **Long-Trend Plots** — cumulative identification and quantification trends

### Contaminants

- **Potential Contaminants** — per-sample percentage of signal from contaminant proteins
- **Top N Contaminants** — ranked list of the most abundant contaminant proteins detected

### Heatmap

A sample-by-sample correlation heatmap based on peptide or protein quantification values, useful for detecting batch effects and outlier samples.

## Notes

- The quantms module requires at least one mzTab file to produce identification-level statistics.
- DIA experiments are handled by the dedicated [quantmsdiann pipeline](https://github.com/bigbio/quantmsdiann) and the `--diann-plugin`. See the [DIA-NN Plugin Guide](diann-plugin.md).
- For TMT experiments, MSstatsTMT input files are parsed for normalization and quantification summaries.
