# Using pmultiqc with DIA-NN Output

The `diann` plugin processes results from [DIA-NN](https://github.com/vdemichev/DiaNN) for data-independent acquisition (DIA) proteomics experiments.

!!! warning "DIA analysis uses quantmsdiann, not quantms"
    DIA experiments are processed by the dedicated **[quantmsdiann](https://github.com/bigbio/quantmsdiann)** pipeline — not quantms. The quantms pipeline handles DDA workflows only (LFQ, TMT/iTRAQ). If you ran DIA-NN through quantmsdiann, your results will contain `report.tsv`/`report.parquet` files — not mzTab.

## Supported Input Files

DIA-NN output does **not** use mzTab format. The input files are DIA-NN native reports:

| File Pattern | Description | Source |
|---|---|---|
| `report.tsv` or `report.parquet` | DIA-NN main report with precursor-level results | DIA-NN / quantmsdiann |
| `report.log.txt` or `diannsummary.log` | DIA-NN log file with software version and parameters | DIA-NN / quantmsdiann |
| `*sdrf.tsv` | SDRF-Proteomics metadata file (optional) | User-provided |
| `*ms_info.parquet` | Precomputed mzML statistics (optional) | quantms-utils |

!!! info "No mzTab for DIA"
    Unlike the quantms DDA pipeline which produces mzTab output, DIA-NN and quantmsdiann produce `report.tsv` / `report.parquet` as the primary output. pmultiqc reads these directly — there is no mzTab conversion step for DIA experiments.

## Running the Report

### From quantmsdiann pipeline output

```bash
# Point to the quantmsdiann results directory
multiqc --diann-plugin /path/to/quantmsdiann/results -o ./report
```

### From standalone DIA-NN

```bash
# Point to the directory containing report.tsv
multiqc --diann-plugin /path/to/diann/output -o ./report
```

### Common options

```bash
# Disable hover tooltips for static/publication-ready reports
multiqc --diann-plugin /path/to/results -o ./report --disable-hoverinfo

# Disable large peptide/protein tables (faster for large experiments)
multiqc --diann-plugin /path/to/results -o ./report --disable-table

# Custom contaminant prefix
multiqc --diann-plugin /path/to/results -o ./report --contaminant-affix CONT
```

## DIA-Specific Metrics

DIA-NN report files contain precursor-level data with DIA-specific confidence metrics (Q.Value, PG.Q.Value). The following QC sections are generated:

### DIA-NN Metadata Table

Software version, library type (in-silico predicted vs. empirical), and key search parameters extracted from the log file:

- DIA-NN version
- Spectral library source
- Mass accuracy settings (MS1 and MS2)
- Protein inference method (genes vs. protein groups)

### Experimental Design

When an SDRF file is present, an experimental design table links DIA raw files to samples, conditions, biological replicates, and fractions.

### Identification Statistics

- **Identification Summary Table** — per-sample counts of precursors, peptides, and protein groups passing FDR
- **MS/MS Identification Rate** — fraction of MS2 events matched to library entries, by run
- **Number of Peptides per Protein** — unique peptide distribution per protein group

### Precursor-Level Metrics

- **Charge State Distribution** — precursor charge frequencies (+1 through +6+)
- **Peptide Length Distribution** — sequence length histogram
- **Peak Intensity Distribution** — precursor intensities (log scale)
- **Peaks per MS2 Spectrum** — fragment ion peak count histogram

### MS1 Information (when `ms_info.parquet` is available)

- **Total Ion Chromatograms (TIC)** — full-gradient signal per sample
- **Base Peak Chromatograms (BPC)** — most abundant ion at each RT point
- **MS1 Peak Count** — number of survey scans per run

### Long-Trend Analysis

Cumulative identification plots across the retention time gradient — useful for detecting chromatographic issues.

## DIA vs DDA Output Comparison

| Feature | DIA (quantmsdiann) | DDA (quantms) |
|---------|-------------------|---------------|
| Primary output | `report.tsv` / `report.parquet` | mzTab |
| pmultiqc plugin | `--diann-plugin` | `--quantms-plugin` |
| Quantification level | Precursor (MS2) | Feature (MS1) |
| FDR columns | Q.Value, PG.Q.Value | PSM-level PEP |
| Library | In-silico or empirical | Not applicable |
| Pipeline | quantmsdiann | quantms |

## File Format Notes

- DIA-NN produces both `.tsv` and `.parquet` report formats. The `.parquet` format is preferred for large experiments (>50 samples) as it loads significantly faster.
- When running through quantmsdiann, the report is typically named `diann_report.tsv` or `diann_report.parquet`.
- Fragment-level filtering uses Q.Value and PG.Q.Value thresholds from the DIA-NN log.
- Bruker `.d` raw formats are detected automatically; chromatographic metrics may differ slightly from Thermo `.raw` data.
