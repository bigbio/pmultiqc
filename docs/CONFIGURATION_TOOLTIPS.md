# Configuration: disable_plot_tooltips

## Overview

This feature adds a global configuration option `disable_plot_tooltips` that allows users to disable tooltips across all interactive plots in pmultiqc.

## Configuration

Add to your `multiqc_config.yaml`:

```yaml
disable_plot_tooltips: true
```

**Default**: `false` (tooltips enabled)

## Affected Plots

This configuration affects all plot types including:

- **Box plots**: Intensity distributions, mass error distributions, intensity standard deviation
- **Bar graphs**: Peptide counts, protein groups, modifications, charge distributions
- **Line graphs**: Retention time plots, peak width over RT, normalization factor over RT, IDs over RT
- **Scatter plots**: PCA plots
- **Heatmaps**: Pipeline performance overview

## Implementation

The feature is implemented using:

1. **Configuration variable**: Set in `pmultiqc/main.py` via `config.disable_plot_tooltips`
2. **Helper module**: `pmultiqc/modules/common/tooltip_config.py` provides:
   - `should_disable_tooltips()`: Check if tooltips should be disabled
   - `apply_tooltip_config(plot_config)`: Apply configuration to plot dictionaries

3. **Updated files**:
   - `pmultiqc/main.py`: Configuration initialization
   - `pmultiqc/modules/common/tooltip_config.py`: Helper functions
   - `pmultiqc/modules/maxquant/maxquant_plots.py`: MaxQuant plots
   - `pmultiqc/modules/proteobench/proteobench_utils.py`: ProteoBench plots  
   - `pmultiqc/modules/common/plots/dia.py`: DIA-specific plots
   - `pmultiqc/modules/common/plots/id.py`: Identification plots
   - `docs/README.md`: Documentation
   - `README.md`: Documentation

## Usage Example

```bash
# Create config file
echo "disable_plot_tooltips: true" > multiqc_config.yaml

# Run pmultiqc
multiqc /path/to/data --config multiqc_config.yaml -o ./report
```

## Context

Some users find the interactive tooltips (especially in box-and-whisker plots) to be overwhelming. This configuration option allows them to disable tooltips when desired while keeping them enabled by default for users who find them helpful.
