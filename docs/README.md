# pmultiqc Examples

This directory contains examples for the pmultiqc project.

## Plugins:

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
| `--mzid_plugin` | Generate reports based on mzIdentML files | `False` |

## Example Data
The example data information can be found in the `config.json`.

## Example reports:

- [LFQ Example](LFQ_PXD007683/multiqc_report.html)
- [TMT Example](TMT_PXD007683/multiqc_report.html)
- [DIA Example](dia/multiqc_report.html)
- [MaxQuant Example](PXD003133/multiqc_report.html)
- [mzIdentML with mzML Example](PXD051187/multiqc_report.html)
- [mzIdentML with MGF Example](PXD054720/multiqc_report.html)

## License

This project is licensed under the terms of the LICENSE file included in the repository.

## Citation

If you use pmultiqc in your research, please cite:

pmultiqc: A MultiQC plugin for proteomics quality control
https://github.com/bigbio/pmultiqc