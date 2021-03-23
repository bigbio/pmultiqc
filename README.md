# pmultiqc
[![Python application](https://github.com/bigbio/pmultiqc/actions/workflows/python-app.yml/badge.svg?branch=main)](https://github.com/bigbio/pmultiqc/actions/workflows/python-app.yml)

A library for proteomics QC report based on MultiQC framework. The library generates a QC report for the [proteomicsLFQ pipeline](https://github.com/nf-core/proteomicslfq). The library read the input of the proteomicsLFQ pipeline, with the following structure:


- consensus_ids       : Identification results from ConsesusId tool in OpenMS
- dbs                 : Database used for the peptide/protein identification step.
- ids                 : Identification results from each search engine.
- logs                : Log files for each independent step
- pipeline_info       : Pipeline info.
- proteomics_lfq      : Final results of the pipeline
  - out.consensusXML  : Feature map output of OpenMS including non-id features.
  - out.mzTab         : mzTab with results of the identification
  - out_msstats.csv   : Input of MSstats software
  - out_triqler.tsv   : Input of Triqler software
- raw_ids             : Identification results from search + percolator

An example report can be found in [multiqc_report.html](http://bigbio.xyz/pmultiqc/multiqc_report.html)

Most of the metrics are compute based on the `out.mzTab` and the `consensus_ids` which contains the filtered peptides and protein identifications.

## Metrics

### General report

- First we show the experimental design of the dataset project http://bigbio.xyz/pmultiqc/multiqc_report.html#proteomicslfq_exp_design . This is a translation from the SDRF proteomics standard to OpenMS standard configuration.
- **Summary Table**: shows the number of spectra, % of identified spectra, total peptide count, total identified proteins (including protein groups - if two proteins are identified by the same peptide the two proteins are count) http://bigbio.xyz/pmultiqc/multiqc_report.html#proteomicslfq_summary_table

- **Number of peptides per Protein**: Includes an histogram with the number of peptides per proteins http://bigbio.xyz/pmultiqc/multiqc_report.html#num_of_pep_per_prot

### Results tables

Two tables are shown to the user with the first [500 peptides](http://bigbio.xyz/pmultiqc/multiqc_report.html#quant_result) in the mzTab and the first [500 PSMs](http://bigbio.xyz/pmultiqc/multiqc_report.html#psm). This tables enable to show some of the most relevant peptide and PSMs in the experiment.

### Identification Statistics

A table called [Spectra Tracking](http://bigbio.xyz/pmultiqc/multiqc_report.html#spectra_tracking) summarize the Identification results by mzML file. The table capture the following numbers:

- MS1_num: Number of MS1 in the mzML
- MS2_num: Number of MS2 in the mzML
- MSGF: Number of Peptides identified using the MSGF+ search engine
- Comet: Number of Peptides identified using the Comet search engine
- Final result of Spectra: Final number of PSMs reported in the mzTab?
- Final result of Peptides: Final number of Peptides identified in the mzTab






