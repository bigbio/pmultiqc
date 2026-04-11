# CLI Reference

pmultiqc extends the MultiQC command-line interface with proteomics-specific options. All options are passed as arguments to the `multiqc` command after installing pmultiqc.

## Basic Invocation

```bash
multiqc [ANALYSIS_DIR] -o [OUTPUT_DIR] [OPTIONS]
```

Each plugin is activated by passing its corresponding flag. Only one plugin should be active per run.

```bash
# quantms pipeline output
multiqc --quantms-plugin /path/to/results -o ./report

# MaxQuant output
multiqc --maxquant-plugin /path/to/results -o ./report

# DIA-NN output
multiqc --diann-plugin /path/to/results -o ./report

# mzIdentML files
multiqc --mzid-plugin /path/to/results -o ./report

# ProteoBench result file
multiqc --proteobench-plugin /path/to/results -o ./report

# FragPipe output
multiqc --fragpipe-plugin /path/to/results -o ./report

# nf-core/mhcquant output
multiqc --mhcquant-plugin /path/to/results -o ./report
```

## Plugin Selection Options

| Option | Description | Default |
|---|---|---|
| `--quantms-plugin` | Enable the quantms module | `False` |
| `--maxquant-plugin` | Enable the MaxQuant module | `False` |
| `--diann-plugin` | Enable the DIA-NN module | `False` |
| `--mzid-plugin` | Enable the mzIdentML module | `False` |
| `--proteobench-plugin` | Enable the ProteoBench module | `False` |
| `--fragpipe-plugin` | Enable the FragPipe module | `False` |
| `--mhcquant-plugin` | Enable the mhcquant module | `False` |

## Identification Options

These options control how PSMs and peptides are counted and filtered.

| Option | Description | Default |
|---|---|---|
| `--remove-decoy` / `--no-remove-decoy` | Remove decoy peptides from counts | `True` (remove) |
| `--decoy-affix TEXT` | Prefix or suffix identifying decoy proteins | `DECOY_` |
| `--contaminant-affix TEXT` | Prefix or suffix identifying contaminant proteins | `CONT` |
| `--affix-type [prefix\|suffix]` | Whether the decoy/contaminant affix is a prefix or suffix | `prefix` |

Example: if your search used `REV_` as a decoy suffix:

```bash
multiqc --quantms-plugin results/ -o report/ \
  --decoy-affix REV_ --affix-type suffix
```

## Quantification Options

| Option | Description | Default |
|---|---|---|
| `--quantification-method [feature_intensity\|spectral_counting]` | LFQ quantification method | `feature_intensity` |
| `--condition TEXT` | Column name to use for condition grouping (e.g., SDRF factor column) | — |
| `--keep-raw` | Keep raw filenames in the experimental design output table | `False` |

Example: create conditions from an SDRF factor column:

```bash
multiqc --quantms-plugin results/ -o report/ \
  --condition "factor value[disease]"
```

## Performance Options

| Option | Description | Default |
|---|---|---|
| `--disable-table` | Skip rendering protein and peptide tables (recommended for >100 samples) | `False` |
| `--ignored-idxml` | Skip parsing idXML files (faster for large quantms outputs) | `False` |

## Display Options

| Option | Description | Default |
|---|---|---|
| `--disable-hoverinfo` | Disable interactive hover tooltips on all plots | `False` |

Use `--disable-hoverinfo` when generating static reports for PDF export or archiving.

## Versioning

```bash
# Print the pmultiqc version
multiqc --pmultiqc-version
```

## Standard MultiQC Options

pmultiqc inherits all standard MultiQC options. Commonly used ones in proteomics workflows:

```bash
# Set report title
multiqc --quantms-plugin results/ -o report/ --title "PXD012345 LFQ QC"

# Set report filename
multiqc --quantms-plugin results/ -o report/ --filename my_report.html

# Force overwrite of existing report
multiqc --quantms-plugin results/ -o report/ --force

# Verbose logging
multiqc --quantms-plugin results/ -o report/ --verbose

# Limit search depth to avoid scanning irrelevant subdirectories
multiqc --quantms-plugin results/ -o report/ --dirs-depth 2
```

Refer to the [MultiQC documentation](https://multiqc.info/docs/) for the full list of standard options.
