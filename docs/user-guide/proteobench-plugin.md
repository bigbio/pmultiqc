# Using pmultiqc with ProteoBench

The `proteobench` plugin integrates with [ProteoBench](https://proteobench.readthedocs.io), a community benchmarking platform for proteomics data analysis workflows. It visualizes the standardized result file produced by ProteoBench to compare quantification performance across conditions.

## Supported Input Files

| File | Description |
|---|---|
| `result_performance.csv` (or `.tsv`, `.txt`) | ProteoBench result file with precursor-level quantification |

The file must contain columns produced by ProteoBench's `result_performance` module. Key expected columns include precursor ion identifiers, charge states, condition-grouped intensities (columns containing `abundance_` or `_Condition_`), and computed fold-change statistics.

## Running the Report

```bash
# Basic usage — point to the directory containing the ProteoBench result file
multiqc --proteobench-plugin /path/to/proteobench/results -o ./report

# Only one result file is allowed per run; multiple files will raise an error
multiqc --proteobench-plugin /path/to/proteobench/results -o ./report --disable-hoverinfo
```

## What ProteoBench Measures

ProteoBench standardizes benchmarking by quantifying a two-condition mixture (Condition A and Condition B) with known fold changes. The result file captures precursor-level quantification across replicates of both conditions.

### Precursor Ion Charge Distribution

A bar chart showing the distribution of precursor charge states across all detected precursors. This provides a first-pass view of sample complexity and instrument performance.

### Log2 Mean Intensity Distribution

Line graph and bar chart showing:
- Distribution of mean log2-transformed precursor intensities for Condition A and Condition B
- Number of missing (NA) intensity values per condition

Missing values indicate precursors detected in one condition but not the other, reflecting sensitivity and reproducibility of the quantification workflow.

### Per-Run Intensity Distribution

Log2 intensity distributions per individual raw file (columns matching `abundance_` or `_Condition_`). This reveals run-to-run variability and potential batch effects across replicates.

### Intensity Count per File

Bar chart of the number of quantified precursors in each raw file. Significant variation across replicates suggests data quality issues in specific runs.

### Log2 Intensity Standard Deviation

Distribution of intra-condition standard deviations (log scale). Lower standard deviation indicates higher quantification reproducibility within each condition.

### Coefficient of Variation (CV)

Distribution of CVs for precursor intensities across replicates. Median CV is a standard metric for comparing quantification robustness across methods and pipelines.

### Log2 Fold Change (A vs. B)

Distribution of observed log2 fold changes between Condition A and Condition B. For standard ProteoBench benchmarks (e.g., LFQ benchmark dataset), the expected fold changes are known; this plot allows direct comparison of observed vs. expected ratios.

### Epsilon (Deviation from Expected Fold Change)

Difference between observed and expected log2 fold change per precursor. Epsilon close to 0 indicates accurate quantification; systematic offsets suggest normalization or calibration issues.

### Log2FC vs. Log2 Mean Intensity (MA-Style Plot)

Scatter plot of log2 fold change against mean log2 intensity across conditions. Intensity-dependent fold change bias (common in low-abundance precursors) is visible as a curved or tilted trend in this plot.

## Interpreting Results

| Metric | Better Performance |
|---|---|
| Missing value count | Lower |
| CV (intra-condition) | Lower |
| Epsilon | Closer to 0 |
| Log2FC deviation from expected | Smaller spread |
| Precursors per file | Higher and consistent |

## Notes

- The ProteoBench plugin expects exactly one result file in the target directory. If multiple files are detected, pmultiqc raises a `ValueError` and exits.
- The module is implemented in `pmultiqc/modules/proteobench/proteobench.py` with utility functions in `proteobench_utils.py`.
- ProteoBench result files from any workflow (MaxQuant, DIA-NN, FragPipe, etc.) are supported as long as the column naming convention is consistent with ProteoBench output format.
