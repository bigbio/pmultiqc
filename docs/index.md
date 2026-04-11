# pmultiqc

**pmultiqc** is a [MultiQC](https://multiqc.info/) plugin for comprehensive quality control reporting of proteomics data. It generates interactive HTML reports with visualizations and metrics to help you assess the quality of your mass spectrometry-based proteomics experiments.

## Key Features

- Works with multiple proteomics data formats and analysis pipelines
- Generates interactive HTML reports with rich visualizations
- Provides comprehensive QC metrics for MS data
- Supports different quantification methods (LFQ, TMT, DIA)
- Integrates seamlessly with the [MultiQC](https://multiqc.info/) framework

## Installation

### From PyPI (recommended)

```bash
pip install pmultiqc
```

### From Source

```bash
git clone https://github.com/bigbio/pmultiqc.git
cd pmultiqc
pip install .
```

## Basic Usage

pmultiqc is a plugin for [MultiQC](https://multiqc.info/). After installation, use it via the MultiQC command-line interface:

```bash
multiqc {analysis_dir} -o {output_dir}
```

Where:

- `{analysis_dir}` is the directory containing your proteomics data files
- `{output_dir}` is the directory where you want to save the report

## Quick Examples

=== "quantms Pipeline"

    ```bash
    multiqc --quantms-plugin /path/to/quantms/results -o ./report
    ```

=== "MaxQuant"

    ```bash
    multiqc --maxquant-plugin /path/to/maxquant/results -o ./report
    ```

=== "DIA-NN"

    ```bash
    multiqc --diann-plugin /path/to/diann/results -o ./report
    ```

=== "mzIdentML"

    ```bash
    multiqc --mzid-plugin /path/to/mzid/files -o ./report
    ```

=== "ProteoBench"

    ```bash
    multiqc --proteobench-plugin /path/to/proteobench/files -o ./report
    ```

## Common Options

```bash
# Remove decoy peptides and set condition column
multiqc --quantms-plugin /path/to/results \
    --remove-decoy \
    --condition factor \
    -o ./report

# Disable interactive tooltips for large datasets
multiqc --quantms-plugin /path/to/results \
    --disable-hoverinfo \
    -o ./report
```

## Part of the quantms Ecosystem

pmultiqc is part of the [quantms ecosystem](https://quantms.org) for quantitative proteomics analysis. It works seamlessly with the [quantms pipeline](https://docs.quantms.org), [mokume](https://mokume.quantms.org) library, and [qpx](https://qpx.quantms.org) format tools.

## Citation

> Yue QX, Dai C, Kamatchinathan S, et al. **pmultiqc: An open-source, lightweight, and metadata-oriented QC reporting library for MS proteomics.** *bioRxiv* (2025). [doi:10.1101/2025.11.02.685980](https://doi.org/10.1101/2025.11.02.685980)
