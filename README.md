# pmultiqc
[![Python application](https://github.com/bigbio/pmultiqc/actions/workflows/python-app.yml/badge.svg?branch=main)](https://github.com/bigbio/pmultiqc/actions/workflows/python-app.yml)
[![Upload Python Package](https://github.com/bigbio/pmultiqc/actions/workflows/python-publish.yml/badge.svg)](https://github.com/bigbio/pmultiqc/actions/workflows/python-publish.yml)

A library for proteomics QC report based on MultiQC framework. The library generates a QC report for the [quantms pipeline](https://github.com/nf-core/quantms). The library read the input of the quantms pipeline by specified analysis dir, with the following structure:

- analysis_dir        : Final results of the pipeline
  - exp_design        : experimental design in two-table format
  - out.mzTab         : mzTab with results of the identification
  - *.mzML            : mzML spectra files
  - *.idXML           : Identification results from search + percolator
  - *.yml             : summary software information and parameters of quantms pipeline (optional)
  - diann_report.tsv  : DIA-NN main report file. Only for DIA analysis.

## Usage
```multiqc {analysis_dir} -o {output dir}```

example: ```multiqc resources/LFQ -o ./```

### parameters
- --exp_design: The experimental design file path, the most entries can be derived from the sdrf file
- --sdrf: Sample and Data Relationship Format file path
- --raw: Keep filenames in experimental design output as raw when exp_design file is provided
- --condition: Create conditions from provided (e.g., factor) columns when exp_design file is provided
- --remove_decoy: Whether to remove the decoy peptides when counting
- --decoy_affix: Pre- or suffix of decoy proteins in their accession
- --affix_type: Location of the decoy marker string in the fasta accession. Before (prefix) or after (suffix)
- --disable_plugin: disable pmultiqc plugin


An example report can be found in [multiqc_report.html](http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html)

Most of the metrics are compute based on the `out.mzTab` and the `raw_ids` which contains the peptides and protein identifications.

## Metrics

### General report

- First we show the experimental design of the dataset project http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#proteomicslfq_exp_design . This is a translation from the SDRF proteomics standard to OpenMS standard configuration.
- **Pipeline performance overview**: show the quantms performance overview including.
    - Contaminants Score: as fraction of summed intensity with 0 = sample full of contaminants; 1 = no contaminants
    - Peptide Intensity: linear scale of the median intensity reaching the threshold, i.e. reaching 2^21 of 2^23 gives score 0.25.
    - Charge Score: deviation of the charge 2 proportion from a representative Raw file (median). For typtic digests, peptides of charge 2 (one N-terminal and one at tryptic C-terminal R or K residue) should be dominant. Ionization issues (voltage?), in-source fragmentation, missed cleavages and buffer irregularities can cause a shift (see Bittremieux 2017, DOI: 10.1002/mas.21544 ).
    - Missed Cleavages: the fraction (0% - 100%) of fully cleaved peptides per Raw file
    - Missed Cleavages Var: each raw file is scored for its deviation from the ‘average’ digestion state of the        current study.
    - ID rate over RT: the number of identifications (here, after FDR filtering) over time. Scored using ‘Uniform’ scoring function.
    - MS2 OverSampling: The percentage of non-oversampled 3D-peaks.
    - Pep Missing Value: Linear scale of the fraction of missing peptides.

- **Summary Table**: shows the number of spectra, % of identified spectra, total peptide count, total identified proteins (including protein groups - if two proteins are identified by the same peptide the two proteins are count) http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#proteomicslfq_summary_table

- **Pipeline Results Statistics**: shows quantms pipeline final results, total peptide identified, total identified proteins et al (The data comes from mzTab and the experimental design file).
  
- **Number of peptides per Protein**: Includes an histogram with the number of peptides per proteins http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#num_of_pep_per_prot

### Results tables

Two tables are shown to the user with the first [500 peptides](http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#quant_result) in the mzTab and the first [500 PSMs](http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#psm). This tables enable to show some of the most relevant peptide and PSMs in the experiment.

### Identification Statistics

A table called [Spectra Tracking](http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#spectra_tracking) summarize the Identification results by mzML file. The table capture the following numbers:

- MS1_num: Number of MS1 in the mzML
- MS2_num: Number of MS2 in the mzML
- MSGF: Number of Peptides identified using the MSGF+ search engine
- Comet: Number of Peptides identified using the Comet search engine
- Final result of Spectra: Final number of PSMs reported in the mzTab
- Final result of Peptides: Final number of Peptides identified in the mzTab

### Precursor Charges Distribution
The [Precursor Charges Distribution](http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#Distribution_of_precursor_charges) aims to show the distribution of the precursor ion charges for a given whole experiment, but also for the identified spectra and unidentified spectra. This information can be used to identify potential ionization problems including many 1+ charges from an ESI ionization source or an unexpected distribution of charges. MALDI experiments are expected to contain almost exclusively 1+ charged ions. An unexpected charge distribution may furthermore be caused by specific search engine parameter settings such as limiting the search to specific ion charges.

### Number of Peaks per MS/MS spectrum
The [Number of Peaks per MS/MS spectrum](https://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#Number_of_Peaks_per_MS_MS_spectrum) aims to show the number of peaks per MS/MS spectrum in a given experiment. Too few peaks can identify poor fragmentation or a detector fault, as opposed to a large number of peaks representing very noisy spectra. This chart is extensively dependent on the pre-processing steps performed to the spectra (centroiding, deconvolution, peak picking approach, etc).

### Peak Intensity Distribution

The [Peak Intensity Distribution](http://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#Peak_Intensity_Distribution) aims to show the Peak instensity in the MS2 spectra for all the experiment but also for the identified spectra. The plot split the intesity in chunks of 0-10, 10-100, 100-300, ... 6k-10k, >10k.

This is a histogram representing the ion intensity vs. the frequency for all MS2 spectra in a whole given experiment. It is possible to filter the information for all, identified and unidentified spectra. This plot can give a general estimation of the noise level of the spectra. Generally, one should expect to have a high number of low intensity noise peaks with a low number of high intensity signal peaks. A disproportionate number of high signal peaks may indicate heavy spectrum pre-filtering or potential experimental problems. In the case of data reuse this plot can be useful in identifying the requirement for pre-processing of the spectra prior to any downstream analysis. The quality of the identifications is not linked to this data as most search engines perform internal spectrum pre-processing before matching the spectra. Thus, the spectra reported are not necessarily pre-processed since the search engine may have applied the pre-processing step internally. This pre-processing is not necessarily reported in the experimental metadata.

### Oversampling Distribution
The [Oversampling Distribution] aims to show the OverSampling information. An oversampled 3D-peak is defined as a peak whose peptide ion (same sequence and same charge state) was identified by at least two distinct MS2 spectra in the same Raw file. For high complexity samples, oversampling of individual 3D-peaks automatically leads to undersampling or even omission of other 3D-peaks, reducing the number of identified peptides. Oversampling occurs in low-complexity samples or long LC gradients, as well as undersized dynamic exclusion windows for data independent acquisitions.

### Delta Mass

The [Delta Mass](https://bigbio.xyz/pmultiqc/shared-peptides-star-align-stricter-pep-protein-FDR/multiqc_report.html#delta_mass-1) aims to show the Peak instensity in the MS2 spectra for all the experiment but also for the identified spectra. The plot split the intesity in chunks of 0-10, 10-100, 100-300, ... 6k-10k, >10k. Mass deltas close to zero reflect more accurate identifications and also that the reporting of the amino acid modifications and charges have been done accurately. This plot can highlight systematic bias if not centered on zero. Other distributions can reflect modifications not being reported properly. Also it is easy to see the different between the target and the decoys identifications.

```Note: Because DIA-NN has much difference in output file !!! So some metrics are difficult to calculate```

```Note: If you want to disable this plugin and use the multiqc function, please set disable_plugin```