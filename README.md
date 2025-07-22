# pmultiqc

[![Python application](https://github.com/bigbio/pmultiqc/actions/workflows/python-app.yml/badge.svg?branch=main)](https://github.com/bigbio/pmultiqc/actions/workflows/python-app.yml)
[![Upload Python Package](https://github.com/bigbio/pmultiqc/actions/workflows/python-publish.yml/badge.svg)](https://github.com/bigbio/pmultiqc/actions/workflows/python-publish.yml)
![PyPI - Version](https://img.shields.io/pypi/v/pmultiqc?style=flat)
![PyPI - Downloads](https://img.shields.io/pypi/dm/pmultiqc)
![Pepy Total Downloads](https://img.shields.io/pepy/dt/pmultiqc)
![GitHub Repo stars](https://img.shields.io/github/stars/bigbio/pmultiqc)

## What is pmultiqc?

pmultiqc is a MultiQC plugin for comprehensive quality control reporting of proteomics data. It generates interactive HTML reports with visualizations and metrics to help you assess the quality of your mass spectrometry-based proteomics experiments.

### Key Features

- Works with multiple proteomics data formats and analysis pipelines
- Generates interactive HTML reports with visualizations
- Provides comprehensive QC metrics for MS data
- Supports different quantification methods (LFQ, TMT, DIA)
- Integrates with the MultiQC framework

## Supported Data Sources

pmultiqc supports the following data sources:

1. **[quantms pipeline](https://github.com/nf-core/quantms)** output files:
   - `experimental_design.tsv`: Experimental design file
   - `*.mzTab`: Results of the identification
   - `*msstats*.csv`: MSstats/MSstatsTMT input files
   - `*.mzML`: Spectra files
   - `*ms_info.tsv`: MS quality control information
   - `*.idXML`: Identification results
   - `*.yml`: Pipeline parameters (optional)
   - `diann_report.tsv` or `diann_report.parquet`: DIA-NN main report (DIA analysis only)

2. **[MaxQuant](https://www.maxquant.org)** result files:
   - `parameters.txt`: Analysis parameters
   - `proteinGroups.txt`: Protein identification results
   - `summary.txt`: Summary statistics
   - `evidence.txt`: Peptide evidence
   - `msms.txt`: MS/MS scan information
   - `msmsScans.txt`: MS/MS scan details

3. **[DIA-NN](https://aptila.bio)** result files:
   - `*ms_info.parquet`: mzML statistics after Raw-to-mzML conversion (using **[quantms-utils](https://github.com/bigbio/quantms-utils)**)
   - `report.tsv` or `report.parquet`: DIA-NN main report

4. **[ProteoBench](https://proteobench.readthedocs.io)** file:
   - `result_performance.csv`: ProteoBench result file

5. **mzIdentML** files:
   - `*.mzid`: Identification results
   - `*.mzML` or `*.mgf`: Corresponding spectra files

## Installation

### Install from PyPI

```bash
# To install the stable release from PyPI:
pip install pmultiqc
```

### Install from Source (Without PyPI)

```bash
# Fork the repository on GitHub

# Clone the repository
git clone https://github.com/your-username/pmultiqc.git
cd pmultiqc

# Install the package locally
pip install .

# Now you can run pmultiqc on your own dataset
```

## Usage

pmultiqc is used as a plugin for MultiQC. After installation, you can run it using the MultiQC command-line interface.

### Basic Usage

```bash
multiqc {analysis_dir} -o {output_dir}
```

Where:
- `{analysis_dir}` is the directory containing your proteomics data files
- `{output_dir}` is the directory where you want to save the report

### Examples

#### For quantms pipeline results

```bash
# Basic usage
multiqc /path/to/quantms/results -o ./report

# With specific options
multiqc /path/to/quantms/results -o ./report --remove_decoy --condition factor
```

#### For MaxQuant results

```bash
multiqc --parse_maxquant /path/to/maxquant/results -o ./report
```

#### For DIA-NN results

```bash
multiqc /path/to/diann/results -o ./report
```

#### For ProteoBench files

```bash
multiqc --parse_proteobench /path/to/proteobench/files -o ./report
```

#### For mzIdentML files

```bash
multiqc --mzid_plugin /path/to/mzid/files -o ./report
```


### Command-line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--raw` | Keep filenames in experimental design output as raw | `False` |
| `--condition` | Create conditions from provided columns | - |
| `--remove_decoy` | Remove decoy peptides when counting | `True` |
| `--decoy_affix` | Pre- or suffix of decoy proteins in their accession | `DECOY_` |
| `--contaminant_affix` | The contaminant prefix or suffix | `CONT` |
| `--affix_type` | Location of the decoy marker (prefix or suffix) | `prefix` |
| `--disable_plugin` | Disable pmultiqc plugin | `False` |
| `--quantification_method` | Quantification method for LFQ experiment | `feature_intensity` |
| `--disable_table` | Disable protein/peptide table plots for large datasets | `False` |
| `--ignored_idxml` | Ignore idXML files for faster processing | `False` |
| `--parse_maxquant` | Generate reports based on MaxQuant results | `False` |
| `--parse_proteobench` | Generate reports based on ProteoBench result | `False` |
| `--mzid_plugin` | Generate reports based on mzIdentML files | `False` |

## QC Metrics and Visualizations

pmultiqc generates a comprehensive report with multiple sections:

### General Report

- **Experimental Design**: Overview of the dataset structure
- **Pipeline Performance Overview**: Key metrics including:
  - Contaminants Score
  - Peptide Intensity
  - Charge Score
  - Missed Cleavages
  - ID rate over RT
  - MS2 OverSampling
  - Peptide Missing Value
- **Summary Table**: Spectra counts, identification rates, peptide and protein counts
- **MS1 Information**: Quality metrics at MS1 level
- **Pipeline Results Statistics**: Overall identification results
- **Number of Peptides per Protein**: Distribution of peptide counts per protein

### Results Tables

- **Peptide Table**: First 500 peptides in the dataset
- **PSM Table**: First 500 PSMs (Peptide-Spectrum Matches)

### Identification Statistics

- **Spectra Tracking**: Summary of identification results by file
- **Search Engine Scores**: Distribution of search engine scores
- **Precursor Charges Distribution**: Distribution of precursor ion charges
- **Number of Peaks per MS/MS Spectrum**: Peak count distribution
- **Peak Intensity Distribution**: MS2 peak intensity distribution
- **Oversampling Distribution**: Analysis of MS2 oversampling
- **Delta Mass**: Mass accuracy distribution
- **Peptide/Protein Quantification Tables**: Quantitative levels across conditions

## Example Reports

You can find example reports on the [docs page](https://bigbio.github.io/pmultiqc).

## Development

To contribute to pmultiqc:

1. Fork the repository
2. Clone your fork: `git clone https://github.com/YOUR-USERNAME/pmultiqc`
3. Create a feature branch: `git checkout -b new-feature`
4. Make your changes
5. Install in development mode: `pip install -e .`
6. Test your changes: `cd tests && multiqc resources/LFQ -o ./`
7. Commit your changes: `git commit -am 'Add new feature'`
8. Push to the branch: `git push origin new-feature`
9. Submit a pull request

## License

This project is licensed under the terms of the LICENSE file included in the repository.

## Citation

If you use pmultiqc in your research, please cite:

```
pmultiqc: A MultiQC plugin for proteomics quality control
https://github.com/bigbio/pmultiqc
```