# Example Reports

Browse interactive QC reports generated from real proteomics experiments.

## LFQ Experiments

| Dataset | Description | Report |
|---------|-------------|--------|
| PXD007683 (LFQ) | Label-free quantification benchmark | [View Report](https://pmultiqc.quantms.org/LFQ_PXD007683/multiqc_report.html){target="_blank"} |

## TMT Experiments

| Dataset | Description | Report |
|---------|-------------|--------|
| PXD007683 (TMT) | TMT isobaric labeling benchmark | [View Report](https://pmultiqc.quantms.org/TMT_PXD007683/multiqc_report.html){target="_blank"} |

## DIA Experiments

| Dataset | Description | Report |
|---------|-------------|--------|
| DIA-NN | DIA-NN analysis results | [View Report](https://pmultiqc.quantms.org/DIANN/multiqc_report.html){target="_blank"} |
| MaxDIA | MaxQuant DIA analysis results | [View Report](https://pmultiqc.quantms.org/MaxDIA/multiqc_report.html){target="_blank"} |

## Special Formats

| Dataset | Description | Report |
|---------|-------------|--------|
| PXD003133 | Large-scale proteomics dataset | [View Report](https://pmultiqc.quantms.org/PXD003133/multiqc_report.html){target="_blank"} |
| PXD051187 (mzIdentML) | mzIdentML format results | [View Report](https://pmultiqc.quantms.org/PXD051187/multiqc_report.html){target="_blank"} |
| ProteoBench | Benchmarking comparison | [View Report](https://pmultiqc.quantms.org/ProteoBench/multiqc_report.html){target="_blank"} |

!!! tip "Generate Your Own"
    Install pmultiqc and run `multiqc --quantms-plugin /path/to/results -o ./report` to generate reports from your own data. See the [Getting Started](../getting-started.md) guide.
