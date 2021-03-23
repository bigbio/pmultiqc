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


Most of the metrics are compute based on the `out.mzTab` and the `consensus_ids` which contains the filtered peptides and protein identifications.

## Metrics




