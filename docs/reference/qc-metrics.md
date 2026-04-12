# QC Metrics Reference

This page lists all quality control metrics computed by pmultiqc, organized by data level. Each metric includes a description, the source file(s) it is derived from, and the modules that produce it.

## PSM-Level Metrics

| Metric | Description | Source | Modules |
|---|---|---|---|
| **PSM Count** | Total number of peptide-spectrum matches passing FDR | mzTab, evidence.txt, report.tsv | quantms, maxquant, diann |
| **MS/MS Identification Rate** | Percentage of MS2 spectra matched to a peptide | summary.txt, mzTab | quantms, maxquant |
| **Search Engine Score** | Distribution of PSM-level confidence scores (e.g., Andromeda score, OpenMS score) | idXML, msms.txt | quantms, maxquant |
| **Delta Mass (Da)** | Difference between experimental and theoretical precursor mass in daltons | mzTab, evidence.txt, mzid | quantms, maxquant, mzidentml |
| **Delta Mass (ppm)** | Mass error expressed in parts per million; calibrated vs. uncalibrated | mzTab, evidence.txt | quantms, maxquant |
| **Precursor Charge State** | Frequency of +1, +2, +3, +4, +5, +6+ charge states among identified PSMs | mzTab, evidence.txt, report.tsv, mzid | all |
| **Retention Time Distribution** | Histogram of PSM identification times across the LC gradient | mzTab, evidence.txt | quantms, maxquant, mzidentml |
| **Oversampling** | Fraction of precursors fragmented more than once in a single run | evidence.txt, mzML | quantms, maxquant, mzidentml |

## Peptide-Level Metrics

| Metric | Description | Source | Modules |
|---|---|---|---|
| **Peptide Count** | Total unique peptide sequences identified per sample | mzTab, evidence.txt, report.tsv | all |
| **Missed Cleavages** | Distribution of peptides with 0, 1, or 2 missed tryptic cleavage sites | mzTab, evidence.txt | quantms, maxquant |
| **Peptide Length Distribution** | Histogram of peptide sequence lengths (number of amino acid residues) | mzTab, evidence.txt, report.tsv | quantms, maxquant, diann |
| **Peptide Intensity** | Log2-transformed peptide-level signal intensity per sample | mzTab, evidence.txt | quantms, maxquant |
| **Missing Values (Peptide)** | Fraction of peptides with no quantification value in a given sample | mzTab | quantms |
| **Modifications** | Summary of variable modifications identified (phospho, oxidation, etc.) | mzTab, evidence.txt | quantms, maxquant |
| **Peptide Table** | First 500 peptides with identification scores and quantification values | mzTab | quantms |

## Protein-Level Metrics

| Metric | Description | Source | Modules |
|---|---|---|---|
| **Protein Group Count** | Total unique protein groups identified per sample | mzTab, proteinGroups.txt, report.tsv | all |
| **Peptides per Protein** | Distribution of unique peptide counts assigned to each protein group | mzTab, proteinGroups.txt | quantms, maxquant, diann, mzidentml |
| **Protein Identification Summary Table** | Per-sample table with PSM, peptide, and protein counts | mzTab, proteinGroups.txt | quantms, maxquant |
| **Contaminant Fraction** | Percentage of total protein signal attributable to common contaminants | proteinGroups.txt, mzTab | quantms, maxquant, fragpipe |
| **Top N Contaminants** | Ranked list of the most abundant contaminant proteins detected | mzTab, proteinGroups.txt | quantms, maxquant |
| **Protein Quantification Table** | LFQ, iBAQ, or TMT intensities per protein group across conditions | mzTab, proteinGroups.txt | quantms, maxquant |

## MS1-Level Metrics

| Metric | Description | Source | Modules |
|---|---|---|---|
| **Total Ion Chromatogram (TIC)** | Sum of all ion intensities at each retention time point | mzML | quantms, diann, mzidentml |
| **Base Peak Chromatogram (BPC)** | Intensity of the most abundant ion at each retention time point | mzML | quantms, diann, mzidentml |
| **MS1 Scan Count** | Total number of survey (MS1) scans per raw file | mzML, ms_info.tsv | quantms, diann |
| **MS2 Scan Count** | Total number of fragmentation (MS2) scans per raw file | mzML, ms_info.tsv | quantms, diann |

## MS2-Level Metrics

| Metric | Description | Source | Modules |
|---|---|---|---|
| **Peaks per MS2 Spectrum** | Histogram of fragment ion peak counts per spectrum | mzML, mgf | quantms, maxquant, mzidentml, diann |
| **Peak Intensity Distribution** | Distribution of MS2 fragment ion intensities (log scale) | mzML, mgf | quantms, maxquant, mzidentml, diann |

## Experimental Design and Sample-Level Metrics

| Metric | Description | Source | Modules |
|---|---|---|---|
| **Experimental Design Table** | Sample-to-file mapping with biological replicate, fraction, and technical replicate | experimental_design.tsv, sdrf.tsv | quantms, maxquant, diann |
| **Sample Correlation Heatmap** | Pairwise sample correlation matrix based on peptide/protein intensities | mzTab, proteinGroups.txt | quantms, maxquant |
| **ID Rate over RT** | Identification rate binned by retention time; normalized per run | mzML + mzTab | quantms, mzidentml |

## ProteoBench-Specific Metrics

| Metric | Description | Source |
|---|---|---|
| **Log2 Mean Intensity** | Mean log2 intensity per precursor, split by Condition A and B | result_performance.csv |
| **Coefficient of Variation (CV)** | Intra-condition CV of precursor intensities | result_performance.csv |
| **Log2 Fold Change (A vs. B)** | Observed fold change between two conditions | result_performance.csv |
| **Epsilon** | Deviation of observed from expected log2 fold change | result_performance.csv |
| **Missing Value Count** | Number of precursors without intensity in each condition | result_performance.csv |
| **Intensity per File** | Per-run precursor counts and intensity distributions | result_performance.csv |

## mhcquant-Specific Metrics

| Metric | Description | Source |
|---|---|---|
| **m/z Histogram** | Distribution of precursor m/z values | multiqc_histogram_mz.txt |
| **RT Histogram** | Distribution of peptide retention times | multiqc_histogram_rt.txt |
| **Ion Mobility Histogram** | Distribution of ion mobility values (Bruker timsTOF) | multiqc_histogram_im.txt |
| **Score Histogram** | Distribution of PSM identification scores | multiqc_histogram_scores.txt |
| **Xcorr Scores** | Sequest cross-correlation score distribution | multiqc_scores_xcorr.txt |
| **Mass Error** | Mass accuracy histogram for identified MHC peptides | multiqc_mass_error.txt |
| **Peptide Intensity** | Quantification signal distribution for MHC-bound peptides | multiqc_peptide_intensity.txt |
| **Percolator Output** | Post-processing score distribution from Percolator | percolator_plot.txt |
| **Peptide Length Distribution** | Sequence length histogram (MHC-I: 8–11 aa; MHC-II: 13–25 aa) | multiqc_length_dist.txt |
