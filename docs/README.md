<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="images/pmultiqc_logo_darkbg.svg">
    <source media="(prefers-color-scheme: light)" srcset="images/pmultiqc_logo.svg">
    <img src="images/pmultiqc_logo.svg" width="40%" alt="pmultiqc Logo"/>
  </picture>
</p>

[![Python application](https://github.com/bigbio/pmultiqc/actions/workflows/python-app.yml/badge.svg?branch=main)](https://github.com/bigbio/pmultiqc/actions/workflows/python-app.yml)
[![Upload Python Package](https://github.com/bigbio/pmultiqc/actions/workflows/python-publish.yml/badge.svg)](https://github.com/bigbio/pmultiqc/actions/workflows/python-publish.yml)
![PyPI - Version](https://img.shields.io/pypi/v/pmultiqc?style=flat)
![PyPI - Downloads](https://img.shields.io/pypi/dm/pmultiqc)
![Pepy Total Downloads](https://img.shields.io/pepy/dt/pmultiqc)
![GitHub Repo stars](https://img.shields.io/github/stars/bigbio/pmultiqc)

## 🚀 What is pmultiqc?

**pmultiqc** is a [MultiQC](https://multiqc.info/) plugin for comprehensive quality control reporting of proteomics data. It generates interactive HTML reports with visualizations and metrics to help you assess the quality of your mass spectrometry-based proteomics experiments.

### ✨ Key Features

- 📊 Works with multiple proteomics data formats and analysis pipelines
- 💻 Generates interactive HTML reports with visualizations
- 📈 Provides comprehensive QC metrics for MS data
- 🔬 Supports different quantification methods (LFQ, TMT, DIA)
- 🧩 Integrates seamlessly with the MultiQC framework
- 🌐 Web service for easy access and PRIDE integration

## 🌐 Public Services

You can use pmultiqc through our public web services:

| Service | URL | Status | Description |
|---------|-----|--------|-------------|
| **EBI PRIDE Service** | [https://www.ebi.ac.uk/pride/services/pmultiqc/](https://www.ebi.ac.uk/pride/services/pmultiqc/) | ![PRIDE pmultiqc service](https://img.shields.io/endpoint?url=https%3A%2F%2Fwww.ebi.ac.uk%2Fpride%2Fservices%2Fpmultiqc%2Fhealth-check)| Official EBI service with PRIDE integration |
| **FU Berlin University Service** | [https://pmultiqc.bsc.fu-berlin.de](https://pmultiqc.bsc.fu-berlin.de) | ![FU Berlin pmultiqc service](https://img.shields.io/endpoint?url=https%3A%2F%2Fpmultiqc.bsc.fu-berlin.de%2Fhealth-check) | pmultiqc service at Freie Universität Berlin |
| **Tübingen University Service** | [https://abi-services.cs.uni-tuebingen.de/pmultiqc/](https://abi-services.cs.uni-tuebingen.de/pmultiqc/) | ![Tübingen pmultiqc service](https://img.shields.io/endpoint?url=https%3A%2F%2Fabi-services.cs.uni-tuebingen.de%2Fpmultiqc%2Fhealth-check) | pmultiqc service at Tübingen University |

### 🎯 Service Features

- **📁 File Upload**: Upload ZIP files with your proteomics data
- **🌐 PRIDE Integration**: Process datasets directly from PRIDE database
- **📋 Job Tracking**: Monitor processing status and download results
- **🔒 Security**: CAPTCHA verification for large file uploads
- **⚡ Multiple Processing**: Handle multiple search engine files from PRIDE
- **📊 Real-time Output**: View console output during processing

## 📋 Supported Data Sources

pmultiqc supports the following data sources:

### 1. [quantms pipeline](https://github.com/nf-core/quantms) output files:
   - `experimental_design.tsv`: Experimental design file
   - `*.mzTab`: Results of the identification
   - `*msstats*.csv`: [MSstats](https://www.bioconductor.org/packages/release/bioc/html/MSstats.html)/MSstatsTMT input files
   - `*.mzML`: Spectra files
   - `*ms_info.tsv`: MS quality control information
   - `*.idXML`: Identification results
   - `*.yml`: Pipeline parameters (optional)
   - `diann_report.tsv` or `diann_report.parquet`: [DIA-NN](https://github.com/vdemichev/DiaNN) main report (DIA analysis only)

### 2. [MaxQuant](https://www.maxquant.org) result files:
   - `parameters.txt`: Analysis parameters
   - `proteinGroups.txt`: Protein identification results
   - `summary.txt`: Summary statistics
   - `evidence.txt`: Peptide evidence
   - `msms.txt`: MS/MS scan information
   - `msmsScans.txt`: MS/MS scan details
   - `*sdrf.tsv`: SDRF-Proteomics (optional)

### 3. **[DIA-NN](https://aptila.bio)** result files:
   - `report.tsv` or `report.parquet`: DIA-NN main report
   - `*sdrf.tsv`: SDRF-Proteomics (optional)
   - `*ms_info.parquet`: mzML statistics after RAW-to-mzML conversion (using **[quantms-utils](https://github.com/bigbio/quantms-utils)**) (optional)

### 4. **[ProteoBench](https://proteobench.readthedocs.io)** file:
   - `result_performance.csv`: ProteoBench result file

### 5. [mzIdentML](https://www.psidev.info/mzidentml) files:
   - `*.mzid`: Identification results
   - `*.mzML` or `*.mgf`: Corresponding spectra files

### 6. [FragPipe](https://fragpipe.nesvilab.org) files:
   - `psm.tsv`: FDR-filtered PSMs
   - `ion.tsv`: FDR-filtered ions
   - `combined_ion.tsv`: FDR-filtered ions
   - `combined_peptide.tsv`: FDR-filtered peptides
   - `combined_protein.tsv`: FDR-filtered proteins

## 💾 Installation

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

## 🛠️ Usage

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
multiqc --quantms-plugin /path/to/quantms/results -o ./report

# With specific options
multiqc --quantms-plugin /path/to/quantms/results -o ./report --remove-decoy --condition factor
```

#### For MaxQuant results

```bash
multiqc --maxquant-plugin /path/to/maxquant/results -o ./report
```

#### For DIA-NN results

```bash
# Discover report inside a results folder
multiqc --diann-plugin /path/to/diann/results -o ./report

# Directly pass a DIA-NN report file (TSV or Parquet)
# Note: MultiQC requires an analysis directory argument; use '.' as a placeholder
multiqc --diann-plugin --diann-report /path/to/report.tsv . -o ./report
multiqc --diann-plugin --diann-report /path/to/report.parquet . -o ./report
```

#### For ProteoBench files

```bash
multiqc --proteobench-plugin /path/to/proteobench/files -o ./report
```

#### For mzIdentML files

```bash
multiqc --mzid-plugin /path/to/mzid/files -o ./report
```

#### For FragPipe files

```bash
multiqc --fragpipe-plugin /path/to/fragpipe/files -o ./report
```

### Command-line Options

| Option | Description | Default |
|---|---|---|
| `--keep-raw` | Keep filenames in experimental design output as raw | `False` |
| `--condition` | Create conditions from provided columns | - |
| `--remove-decoy` | Remove decoy peptides when counting | `True` |
| `--decoy-affix` | Pre- or suffix of decoy proteins in their accession | `DECOY_` |
| `--contaminant-affix` | The contaminant prefix or suffix | `CONT` |
| `--affix-type` | Location of the decoy marker (prefix or suffix) | `prefix` |
| `--disable-plugin` | Disable pmultiqc plugin | `False` |
| `--quantification-method` | Quantification method for LFQ experiment | `feature_intensity` |
| `--disable-table` | Disable protein/peptide table plots for large datasets | `False` |
| `--ignored-idxml` | Ignore idXML files for faster processing | `False` |
| `--quantms-plugin` | Generate reports based on Quantms results | `False` |
| `--diann-plugin` | Generate reports based on DIANN results | `False` |
| `--diann-report` | Path to DIA-NN main report (.tsv or .parquet). When provided with `--diann-plugin`, you can use `.` as the analysis directory placeholder. | - |
| `--maxquant-plugin` | Generate reports based on MaxQuant results | `False` |
| `--proteobench-plugin` | Generate reports based on ProteoBench result | `False` |
| `--mzid-plugin` | Generate reports based on mzIdentML files | `False` |
| `--fragpipe-plugin` | Generate reports based on FragPipe files | `False` |
| `--disable-hoverinfo` | Disable interactive hover tooltips in the plots | `False` |

## 📊 QC Metrics and Visualizations

pmultiqc generates a comprehensive report with multiple sections:

### 📑 General Report

- **Experimental Design**: Overview of the dataset structure
- **Pipeline Performance Overview**: Key metrics including:
  - 🔍 Contaminants Score
  - 📊 Peptide Intensity
  - ⚡ Charge Score
  - ✂️ Missed Cleavages
  - 📈 ID rate over RT
  - 🔄 MS2 OverSampling
  - 🧩 Peptide Missing Value
- **Summary Table**: Spectra counts, identification rates, peptide and protein counts
- **MS1 Information**: Quality metrics at MS1 level
- **Pipeline Results Statistics**: Overall identification results
- **Number of Peptides per Protein**: Distribution of peptide counts per protein

### 📚 Results Tables

- **Peptide Table**: First 500 peptides in the dataset
- **PSM Table**: First 500 PSMs (Peptide-Spectrum Matches)

### 📈 Identification Statistics

- **Spectra Tracking**: Summary of identification results by file
- **Search Engine Scores**: Distribution of search engine scores
- **Precursor Charges Distribution**: Distribution of precursor ion charges
- **Number of Peaks per MS/MS Spectrum**: Peak count distribution
- **Peak Intensity Distribution**: MS2 peak intensity distribution
- **Oversampling Distribution**: Analysis of MS2 oversampling
- **Delta Mass**: Mass accuracy distribution
- **Peptide/Protein Quantification Tables**: Quantitative levels across conditions

## 🔍 Example Reports

You can find example reports on the [docs page](https://bigbio.github.io/pmultiqc). Here are the direct links to different example reports:

| Example Type | Description | Link | Dataset Download |
|---|---|---|---|
| LFQ | Label-free quantification | [LFQ Example](https://pmultiqc.quantms.org/LFQ_PXD007683/multiqc_report.html) ([disable_hoverinfo](https://pmultiqc.quantms.org/LFQ_PXD007683_disable_hoverinfo/multiqc_report.html)) | [LFQ_PXD007683.zip](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/pmultiqc/example-projects/LFQ_PXD007683.zip) |
| TMT | Tandem mass tag | [TMT Example](https://pmultiqc.quantms.org/TMT_PXD007683/multiqc_report.html) ([disable_hoverinfo](https://pmultiqc.quantms.org/TMT_PXD007683_disable_hoverinfo/multiqc_report.html)) | [TMT_PXD007683.zip](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/pmultiqc/example-projects/TMT_PXD007683.zip) |
| quantms DIA | Data-independent acquisition | [quantms DIA Example](https://pmultiqc.quantms.org/dia/multiqc_report.html) ([disable_hoverinfo](https://pmultiqc.quantms.org/dia_disable_hoverinfo/multiqc_report.html)) | [dia.zip](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/pmultiqc/dia/dia.zip) |
| DIA-NN | Data-independent acquisition | [DIA-NN Example](https://pmultiqc.quantms.org/DIANN/multiqc_report.html) ([disable_hoverinfo](https://pmultiqc.quantms.org/DIANN_disable_hoverinfo/multiqc_report.html)) | [PXD063291.zip](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/pmultiqc/example-projects/PXD063291.zip) |
| MaxQuant | MaxQuant results | [MaxQuant Example](https://pmultiqc.quantms.org/PXD003133/multiqc_report.html) ([disable_hoverinfo](https://pmultiqc.quantms.org/PXD003133_disable_hoverinfo/multiqc_report.html)) | [txt_20min.zip](https://ftp.pride.ebi.ac.uk/pride/data/archive/2015/11/PXD003133/txt_20min.zip) |
| MaxQuant DIA | MaxQuant DIA results | [MaxQuant DIA Example](https://pmultiqc.quantms.org/MaxDIA/multiqc_report.html) ([disable_hoverinfo](https://pmultiqc.quantms.org/MaxDIA_disable_hoverinfo/multiqc_report.html)) | [MaxDIA_txt.zip](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/pmultiqc/maxquant/MaxDIA_txt.zip) |
| ProteoBench | ProteoBench results | [ProteoBench Example](https://pmultiqc.quantms.org/ProteoBench/multiqc_report.html) ([disable_hoverinfo](https://pmultiqc.quantms.org/ProteoBench_disable_hoverinfo/multiqc_report.html)) | [ProteoBench data](https://proteobench.cubimed.rub.de/datasets/d01e87b997b84c985868204b1ed26749902fd7f9/d01e87b997b84c985868204b1ed26749902fd7f9_data.zip) |
| mzIdentML with mzML | mzIdentML with mzML files | [mzIdentML with mzML Example](https://pmultiqc.quantms.org/PXD053068/multiqc_report.html) ([disable_hoverinfo](https://pmultiqc.quantms.org/PXD053068_disable_hoverinfo/multiqc_report.html)) | [PXD053068 folder](https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/05/PXD053068/) |
| mzIdentML with MGF | mzIdentML with MGF files | [mzIdentML with MGF Example](https://pmultiqc.quantms.org/PXD054720/multiqc_report.html) ([disable_hoverinfo](https://pmultiqc.quantms.org/PXD054720_disable_hoverinfo/multiqc_report.html)) | [PXD054720 folder](https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/08/PXD054720/) |
| FragPipe | FragPipe results | [FragPipe Example](https://pmultiqc.quantms.org/PXD062399/multiqc_report.html) ([disable_hoverinfo](https://pmultiqc.quantms.org/PXD062399_disable_hoverinfo/multiqc_report.html)) | [PXD062399.zip](https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/pmultiqc/example-projects/PXD062399.zip) |

## 👥 Contributing

To contribute to pmultiqc:

1. 🍴 Fork the repository
2. 📥 Clone your fork: `git clone https://github.com/YOUR-USERNAME/pmultiqc`
3. 🌿 Create a feature branch: `git checkout -b new-feature`
4. ✏️ Make your changes
5. 🔧 Install in development mode: `pip install -e .`
6. 🧪 Test your changes: `cd tests && multiqc resources/LFQ -o ./`
7. 💾 Commit your changes: `git commit -am 'Add new feature'`
8. 📤 Push to the branch: `git push origin new-feature`
9. 📩 Submit a pull request

## 📖 Documentation

For more detailed information, visit the [pmultiqc GitHub repository](https://github.com/bigbio/pmultiqc) or check the [documentation site](https://bigbio.github.io/pmultiqc).

## 📜 License

This project is licensed under the terms of the [LICENSE](https://github.com/bigbio/pmultiqc/blob/main/LICENSE) file included in the repository.

## 📝 Citation

If you use pmultiqc in your research, please cite our preprint:

**pmultiqc: An open-source, lightweight, and metadata-oriented QC reporting library for MS proteomics**  
Qi-Xuan Yue, Chengxin Dai, Selvakumar Kamatchinathan, Chakradhar Bandla, Henry Webel, Asier Larrea, Wout Bittremieux, Julian Uszkoreit, Tom David Muller, Jinqiu Xiao, Juergen Cox, Philip Ewels, Vadim Demichev, Oliver Kohlbacher, Timo Sachsenberg, Chris Bielow, Mingze Bai, Yasset Perez-Riverol  
*bioRxiv* 2025. DOI: [10.1101/2025.11.02.685980](https://doi.org/10.1101/2025.11.02.685980)

Preprint: https://www.biorxiv.org/content/early/2025/11/03/2025.11.02.685980

### BibTeX

```bibtex
@article{Yue2025.11.02.685980,
  author = {Yue, Qi-Xuan and Dai, Chengxin and Kamatchinathan, Selvakumar and Bandla, Chakradhar and Webel, Henry and Larrea, Asier and Bittremieux, Wout and Uszkoreit, Julian and Muller, Tom David and Xiao, Jinqiu and Cox, Juergen and Ewels, Philip and Demichev, Vadim and Kohlbacher, Oliver and Sachsenberg, Timo and Bielow, Chris and Bai, Mingze and Perez-Riverol, Yasset},
  title = {pmultiqc: An open-source, lightweight, and metadata-oriented QC reporting library for MS proteomics},
  elocation-id = {2025.11.02.685980},
  year = {2025},
  doi = {10.1101/2025.11.02.685980},
  publisher = {Cold Spring Harbor Laboratory},
  URL = {https://www.biorxiv.org/content/early/2025/11/03/2025.11.02.685980},
  eprint = {https://www.biorxiv.org/content/early/2025/11/03/2025.11.02.685980.full.pdf},
  journal = {bioRxiv}
}
```

## 🔗 Related Tools

- [MultiQC](https://multiqc.info/) - Base framework for pmultiqc
- [nf-core/quantms](https://nf-co.re/quantms) - Proteomics quantification pipeline
- [OpenMS](https://www.openms.de/) - Open-source tools for mass spectrometry
- [ProteomeXchange](http://www.proteomexchange.org/) - Repository for proteomics data

## ❓ Need Help?

If you have questions or need assistance:
- [Open an issue](https://github.com/bigbio/pmultiqc/issues) on GitHub
- Check [existing issues](https://github.com/bigbio/pmultiqc/issues?q=is%3Aissue) for solutions
