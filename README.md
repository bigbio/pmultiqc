# pmultiqc
[![Python application](https://github.com/bigbio/pmultiqc/actions/workflows/python-app.yml/badge.svg?branch=main)](https://github.com/bigbio/pmultiqc/actions/workflows/python-app.yml)
[![Upload Python Package](https://github.com/bigbio/pmultiqc/actions/workflows/python-publish.yml/badge.svg)](https://github.com/bigbio/pmultiqc/actions/workflows/python-publish.yml)

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

## Usage
```multiqc --exp_design/sdrf {expdesign_file/sdrf file} --mzMLs {mzMLs file dir} --raw_ids {raw identification dir} {proteomicslfq result dir} -o {output dir}```

example: ```multiqc --exp_design ./UPS1/experimental_design.tsv --mzMLs ./UPS1/shared-peptides-star-align-stricter-pep-protein-FDR/mzMLs --raw_ids ./UPS1/shared-peptides-star-align-stricter-pep-protein-FDR/raw_ids ./UPS1/shared-peptides-star-align-stricter-pep-protein-FDR/proteomics_lfq -o ./shared-peptides-star-align-stricter-pep-protein-FDR-statistics ```

### parameters
- --exp_design: The experimental design file path, the most entries can be derived from the sdrf file
- --sdrf: Sample and Data Relationship Format file path
- --raw: Keep filenames in experimental design output as raw when exp_design file is provided
- --condition: Create conditions from provided (e.g., factor) columns when exp_design file is provided
- --quant_method: quantification method (e.g lfq or tmt. default lfq)
- --mzMLs: mzMLs file directory
- --raw_ids: raw identification file dir
- --remove_decoy: Whether to remove the decoy peptides when counting


An example report can be found in [multiqc_report.html](http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html)

Most of the metrics are compute based on the `out.mzTab` and the `consensus_ids` which contains the filtered peptides and protein identifications.

## Metrics

### General report

- First we show the experimental design of the dataset project http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#proteomicslfq_exp_design . This is a translation from the SDRF proteomics standard to OpenMS standard configuration.
- **Summary Table**: shows the number of spectra, % of identified spectra, total peptide count, total identified proteins (including protein groups - if two proteins are identified by the same peptide the two proteins are count) http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#proteomicslfq_summary_table

- **Number of peptides per Protein**: Includes an histogram with the number of peptides per proteins http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#num_of_pep_per_prot

### Results tables

Two tables are shown to the user with the first [500 peptides](http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#quant_result) in the mzTab and the first [500 PSMs](http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#psm). This tables enable to show some of the most relevant peptide and PSMs in the experiment.

### Identification Statistics

A table called [Spectra Tracking](http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#spectra_tracking) summarize the Identification results by mzML file. The table capture the following numbers:

- MS1_num: Number of MS1 in the mzML
- MS2_num: Number of MS2 in the mzML
- MSGF: Number of Peptides identified using the MSGF+ search engine
- Comet: Number of Peptides identified using the Comet search engine
- Final result of Spectra: Final number of PSMs reported in the mzTab?
- Final result of Peptides: Final number of Peptides identified in the mzTab

### Peak Intensity Distribution

The [Peak Intensity Distribution](http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#Peak_Intensity_Distribution) aims to show the Peak instensity in the MS2 spectra for all the experiment but also for the identified spectra. The plot split the intesity in chunks of 0-10, 10-100, 100-300, ... 6k-10k, >10k.

This is a histogram representing the ion intensity vs. the frequency for all MS2 spectra in a whole given experiment. It is possible to filter the information for all, identified and unidentified spectra. This plot can give a general estimation of the noise level of the spectra. Generally, one should expect to have a high number of low intensity noise peaks with a low number of high intensity signal peaks. A disproportionate number of high signal peaks may indicate heavy spectrum pre-filtering or potential experimental problems. In the case of data reuse this plot can be useful in identifying the requirement for pre-processing of the spectra prior to any downstream analysis. The quality of the identifications is not linked to this data as most search engines perform internal spectrum pre-processing before matching the spectra. Thus, the spectra reported are not necessarily pre-processed since the search engine may have applied the pre-processing step internally. This pre-processing is not necessarily reported in the experimental metadata.




