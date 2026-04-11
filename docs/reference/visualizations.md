# Visualizations Reference

pmultiqc generates interactive HTML plots using the MultiQC plotting framework (Plotly-backed). This page describes each plot type used across modules, what it shows, and which sections it appears in.

All plots are rendered in the MultiQC HTML report. Interactive features (hover tooltips, zoom, pan, legend toggle) can be disabled with `--disable-hoverinfo` for static output.

## Bar Charts

Bar charts (`multiqc.plots.bargraph`) are used for count and percentage comparisons across samples.

| Plot | Section | Description |
|---|---|---|
| **MS/MS Identification Rate** | Identification | Percentage of MS2 spectra matched to a peptide, per run and per sample. |
| **Precursor Charge Distribution** | Identification | Frequency of precursor charge states (+1 to +6+) across identified PSMs. |
| **Missed Cleavages** | Identification | Proportion of peptides with 0, 1, or 2 missed tryptic cleavage sites per sample. |
| **Contaminants per Group** | Contaminants | Percentage of total signal from contaminant proteins per sample group. |
| **Top N Contaminants** | Contaminants | Ranked bar chart of the most abundant contaminant proteins. |
| **Precursor Ion Charge (ProteoBench)** | ProteoBench | Charge state distribution from ProteoBench result file. |
| **Missing Values per Condition** | ProteoBench | Count of NA intensity values in Condition A and B. |
| **Intensity Count per File** | ProteoBench | Number of quantified precursors per raw file. |

## Line Graphs / Chromatograms

Line graphs (`multiqc.plots.linegraph`) display continuous signal over a numeric axis such as retention time or intensity.

| Plot | Section | Description |
|---|---|---|
| **Total Ion Chromatogram (TIC)** | MS1 Information | Sum of all ion intensities at each retention time point; one trace per sample. |
| **Base Peak Chromatogram (BPC)** | MS1 Information | Most abundant ion intensity at each retention time point; one trace per sample. |
| **Log2 Mean Intensity Distribution** | ProteoBench | Distribution of mean log2 intensities for Condition A and B as a density-like line. |
| **Log2 Intensity per File** | ProteoBench | Per-run log2 intensity distribution as overlaid line plots. |
| **Log2 Intensity Std Distribution** | ProteoBench | Distribution of intra-condition standard deviations across precursors. |
| **CV Distribution** | ProteoBench | Distribution of precursor coefficient of variation (CV) values. |
| **ID Rate over RT** | Identification | Identification rate per retention time bin; shows where in the gradient IDs are most dense. |
| **Long-Trend Plots** | Identification | Cumulative PSM, peptide, and protein identifications across the RT gradient. |

## Box Plots

Box plots (`multiqc.plots.box`) show the spread of a distribution per sample, providing median, quartiles, and outliers.

| Plot | Section | Description |
|---|---|---|
| **Peptide Intensity Box Plot** | Quantification | Log2 peptide intensities per sample; useful for assessing normalization. |
| **Search Engine Score Distribution** | Identification | PSM-level score distribution per sample; one box per file. |
| **mhcquant Score Distributions** | mhcquant | Score distributions (Xcorr, Percolator) per sample for MHC peptides. |
| **mhcquant Mass Error** | mhcquant | Mass accuracy distribution per sample. |
| **mhcquant Peptide Intensities** | mhcquant | Quantification signal per sample for MHC-bound peptides. |

## Histograms / Distribution Plots

Histograms are implemented as line graphs with binned x-axes or as dedicated histogram elements.

| Plot | Section | Description |
|---|---|---|
| **Delta Mass (ppm and Da)** | Identification | Mass error distribution; a tight peak at 0 indicates good calibration. |
| **Peaks per MS2 Spectrum** | MS Spectra | Frequency distribution of fragment ion peak counts per spectrum. |
| **Peak Intensity Distribution** | MS Spectra | Log2 distribution of MS2 fragment peak intensities. |
| **Peptide Length Distribution** | Identification | Frequency of peptide sequence lengths (number of residues). |
| **mhcquant m/z Histogram** | mhcquant | Precursor m/z distribution for immunopeptidomics data. |
| **mhcquant RT Histogram** | mhcquant | Retention time distribution of identified MHC peptides. |
| **mhcquant Ion Mobility Histogram** | mhcquant | Ion mobility distribution for timsTOF data (when available). |

## Scatter Plots

Scatter plots (`multiqc.plots.scatter`) display relationships between two continuous variables.

| Plot | Section | Description |
|---|---|---|
| **Log2FC vs. Log2 Mean Intensity** | ProteoBench | MA-style plot of fold change vs. mean intensity; reveals intensity-dependent bias. |
| **Epsilon (Deviation) Plot** | ProteoBench | Observed minus expected fold change per precursor; centered at 0 for unbiased quantification. |
| **Log2 A vs. B** | ProteoBench | Pairwise scatter of log2 intensities for Condition A and Condition B. |

## Heatmaps

Heatmaps (`multiqc.plots.heatmap`) display pairwise relationships in a matrix format.

| Plot | Section | Description |
|---|---|---|
| **Sample Correlation Heatmap** | Heatmap | Pairwise Pearson or Spearman correlation of protein/peptide intensities between samples. Color scale from 0 (no correlation) to 1 (perfect correlation). Useful for detecting outlier samples and batch effects. |

## Tables

Interactive sortable/filterable tables (`multiqc.plots.table`) present structured per-sample or per-feature data.

| Plot | Section | Description |
|---|---|---|
| **Experimental Design Table** | Experiment | Sample-to-file mapping with replicate and fraction annotations. |
| **Identification Summary Table** | Identification | Per-sample PSM, peptide, and protein counts with identification rates. |
| **Peptide Table** | Quantification | First 500 peptides with scores and quantification values per condition. |
| **PSM Table** | Quantification | First 500 PSMs with scores, mass errors, and charge states. |
| **MaxQuant Parameters Table** | MaxQuant | Key search parameters extracted from parameters.txt. |
| **DIA-NN Metadata Table** | DIA-NN | Software version, library type, and key run parameters. |
| **mzIdentML Quantification Table** | mzIdentML | Protein-level quantification values when present in the mzid file. |

## Pipeline Performance Overview (Sparklines)

Available in the `quantms` module only. A condensed visual panel with one sparkline per QC dimension per sample:

| Score | What It Summarizes |
|---|---|
| Contaminant Score | Fraction of total intensity from contaminants |
| Peptide Intensity Score | Uniformity of peptide signal across samples |
| Charge Score | Expected vs. observed precursor charge composition |
| Missed Cleavage Score | Proportion of fully cleaved peptides |
| ID Rate / RT Score | Identification coverage across the gradient |
| Oversampling Score | Efficiency of MS2 precursor selection |
| Missing Value Score | Completeness of quantification matrix |
