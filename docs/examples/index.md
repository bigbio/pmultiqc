# Example Reports

Browse interactive QC reports generated from real proteomics experiments reanalyzed with the quantms pipeline. Each report includes identification rates, quantification metrics, and quality control visualizations.

## LFQ Experiments

| Dataset | Description | Interactive Report | Lightweight Report |
|---------|-------------|--------------------|--------------------|
| PXD007683 | Label-free quantification benchmark | [View](https://pmultiqc.quantms.org/LFQ_PXD007683/multiqc_report.html){target="_blank"} | [View](https://pmultiqc.quantms.org/LFQ_PXD007683_disable_hoverinfo/multiqc_report.html){target="_blank"} |
| PXD010899 | Large-scale LFQ proteomics | [View](https://pmultiqc.quantms.org/PXD010899/multiqc_report.html){target="_blank"} | [View](https://pmultiqc.quantms.org/PXD010899_disable_hoverinfo/multiqc_report.html){target="_blank"} |

## TMT Experiments

| Dataset | Description | Interactive Report | Lightweight Report |
|---------|-------------|--------------------|--------------------|
| PXD007683 (TMT) | TMT isobaric labeling benchmark | [View](https://pmultiqc.quantms.org/TMT_PXD007683/multiqc_report.html){target="_blank"} | [View](https://pmultiqc.quantms.org/TMT_PXD007683_disable_hoverinfo/multiqc_report.html){target="_blank"} |
| PXD003133 | Large-scale TMT proteomics | [View](https://pmultiqc.quantms.org/PXD003133/multiqc_report.html){target="_blank"} | [View](https://pmultiqc.quantms.org/PXD003133_disable_hoverinfo/multiqc_report.html){target="_blank"} |
| PXD053068 | TMT multi-plex experiment | [View](https://pmultiqc.quantms.org/PXD053068/multiqc_report.html){target="_blank"} | [View](https://pmultiqc.quantms.org/PXD053068_disable_hoverinfo/multiqc_report.html){target="_blank"} |
| PXD053464 | TMT proteomics dataset | [View](https://pmultiqc.quantms.org/PXD053464/multiqc_report.html){target="_blank"} | [View](https://pmultiqc.quantms.org/PXD053464_disable_hoverinfo/multiqc_report.html){target="_blank"} |
| PXD054720 | TMT quantification dataset | [View](https://pmultiqc.quantms.org/PXD054720/multiqc_report.html){target="_blank"} | [View](https://pmultiqc.quantms.org/PXD054720_disable_hoverinfo/multiqc_report.html){target="_blank"} |

## DIA Experiments

| Dataset | Description | Interactive Report | Lightweight Report |
|---------|-------------|--------------------|--------------------|
| DIA (quantms) | DIA analysis via quantms pipeline | [View](https://pmultiqc.quantms.org/dia/multiqc_report.html){target="_blank"} | [View](https://pmultiqc.quantms.org/dia_disable_hoverinfo/multiqc_report.html){target="_blank"} |
| DIA-NN | DIA-NN native output | [View](https://pmultiqc.quantms.org/DIANN/multiqc_report.html){target="_blank"} | [View](https://pmultiqc.quantms.org/DIANN_disable_hoverinfo/multiqc_report.html){target="_blank"} |
| MaxDIA | MaxQuant DIA analysis | [View](https://pmultiqc.quantms.org/MaxDIA/multiqc_report.html){target="_blank"} | [View](https://pmultiqc.quantms.org/MaxDIA_disable_hoverinfo/multiqc_report.html){target="_blank"} |

## Additional Datasets

| Dataset | Description | Interactive Report | Lightweight Report |
|---------|-------------|--------------------|--------------------|
| PXD062383 | Proteomics QC dataset | [View](https://pmultiqc.quantms.org/PXD062383/multiqc_report.html){target="_blank"} | [View](https://pmultiqc.quantms.org/PXD062383_disable_hoverinfo/multiqc_report.html){target="_blank"} |
| PXD062399 | Proteomics QC dataset | [View](https://pmultiqc.quantms.org/PXD062399/multiqc_report.html){target="_blank"} | [View](https://pmultiqc.quantms.org/PXD062399_disable_hoverinfo/multiqc_report.html){target="_blank"} |
| PXD066146 | Proteomics QC dataset | [View](https://pmultiqc.quantms.org/PXD066146/multiqc_report.html){target="_blank"} | [View](https://pmultiqc.quantms.org/PXD066146_disable_hoverinfo/multiqc_report.html){target="_blank"} |

## Special Formats

| Dataset | Description | Interactive Report | Lightweight Report |
|---------|-------------|--------------------|--------------------|
| ProteoBench | Cross-tool benchmarking comparison | [View](https://pmultiqc.quantms.org/ProteoBench/multiqc_report.html){target="_blank"} | [View](https://pmultiqc.quantms.org/ProteoBench_disable_hoverinfo/multiqc_report.html){target="_blank"} |
| mhcquant | MHC immunopeptidomics | [View](https://pmultiqc.quantms.org/mhcquant/multiqc_report.html){target="_blank"} | [View](https://pmultiqc.quantms.org/mhcquant_disable_hoverinfo/multiqc_report.html){target="_blank"} |

!!! info "Interactive vs Lightweight Reports"
    **Interactive** reports include hover tooltips on all plots for detailed data exploration. **Lightweight** reports disable hover info for faster loading and smaller file size — useful for large experiments or embedding in publications.

!!! tip "Generate Your Own"
    ```bash
    pip install pmultiqc
    multiqc /path/to/results --module pmultiqc -o ./report
    ```
    See the [User Guide](../user-guide/quantms-plugin.md) for detailed instructions per workflow type.

## Citation

If you use pmultiqc for your analysis, please cite:

> Yue QX, Dai C, Kamatchinathan S, Bandla C, Webel H, Larrea A, Bittremieux W, Uszkoreit J, Müller TD, Xiao J, Cox J, Yu F, Ewels P, Demichev V, Kohlbacher O, Sachsenberg T, Bielow C, Bai M, Perez-Riverol Y. **pmultiqc: An open-source, lightweight, and metadata-oriented QC reporting library for MS proteomics.** *Mol Cell Proteomics.* 2026;101530. [DOI: 10.1016/j.mcpro.2026.101530](https://doi.org/10.1016/j.mcpro.2026.101530)
