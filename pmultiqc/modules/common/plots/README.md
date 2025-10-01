# PMultiQC Plots Package

This package contains organized plotting functions for PMultiQC reports. Functions are separated by logical categories for better maintainability and code organization.

## Package Structure

```
plots/
├── __init__.py              # Package initialization and exports
├── contaminants.py          # Contaminant-related plots
├── identification.py        # Identification statistics plots
├── charge_modifications.py  # Charge state and modification plots
├── mass_analysis.py         # Mass error and missed cleavages plots
├── retention_time.py        # Retention time and oversampling plots
├── tables.py               # Table plotting functions
├── ms_analysis.py          # MS1 and MS2 analysis plots
├── general.py              # General utilities and heatmaps
└── README.md               # This file
```

## Migration Guide

### Old Import Pattern
```python
from pmultiqc.modules.common import common_plots

# Usage
common_plots.draw_potential_contaminants(...)
common_plots.HEATMAP_COLOR_LIST
```

### New Import Pattern

#### Option 1: Import from package (Recommended)
```python
from pmultiqc.modules.common.plots import (
    draw_potential_contaminants,
    draw_top_n_contaminants,
    HEATMAP_COLOR_LIST
)

# Usage
draw_potential_contaminants(...)
```

#### Option 2: Import specific modules
```python
from pmultiqc.modules.common.plots.contaminants import (
    draw_potential_contaminants,
    draw_top_n_contaminants
)
from pmultiqc.modules.common.plots.general import HEATMAP_COLOR_LIST

# Usage
draw_potential_contaminants(...)
```

#### Option 3: Import entire modules
```python
from pmultiqc.modules.common.plots import contaminants
from pmultiqc.modules.common.plots import general

# Usage
contaminants.draw_potential_contaminants(...)
general.HEATMAP_COLOR_LIST
```

## Module Categories

### Contaminants (`contaminants.py`)
- `draw_potential_contaminants()` - Potential contaminants per group
- `draw_top_n_contaminants()` - Top N contaminants per raw file
- `remove_subtitle()` - Utility function to remove plot subtitles

### Identification (`identification.py`)
- `draw_ms_ms_identified()` - MS/MS identification rate per raw file
- `draw_quantms_identification()` - Quantms identification plots
- `draw_num_pep_per_protein()` - Number of peptides per protein

### Charge & Modifications (`charge_modifications.py`)
- `draw_charge_state()` - Charge state distribution per file
- `draw_modifications()` - Modifications per raw file
- `draw_precursor_charge_distribution()` - Precursor charge distribution

### Mass Analysis (`mass_analysis.py`)
- `draw_delta_mass_da_ppm()` - Delta mass analysis (Da or ppm)
- `draw_msms_missed_cleavages()` - Missed cleavages per raw file

### Retention Time (`retention_time.py`)
- `draw_ids_rt_count()` - IDs over retention time
- `draw_oversampling()` - MS/MS oversampling plots

### Tables (`tables.py`)
- `draw_peptides_table()` - Peptides quantification table
- `draw_protein_table()` - Protein quantification table
- `draw_exp_design()` - Experimental design table
- `draw_summary_protein_ident_table()` - Summary protein identification table
- `draw_quantms_identi_num()` - Quantms identification numbers table

### MS Analysis (`ms_analysis.py`)
- `draw_ms_information()` - MS1 information plots (TIC, BPC, peaks, stats)
- `draw_peaks_per_ms2()` - Number of peaks per MS/MS spectrum
- `draw_peak_intensity_distribution()` - Peak intensity distribution

### General (`general.py`)
- `draw_heatmap()` - Heatmap plots
- `HEATMAP_COLOR_LIST` - Color configuration for heatmaps

## Benefits of New Structure

1. **Better Organization**: Functions are grouped by logical categories
2. **Easier Maintenance**: Smaller, focused files are easier to maintain
3. **Clearer Dependencies**: Import only what you need
4. **Better Documentation**: Each module can have focused documentation
5. **Easier Testing**: Smaller modules are easier to unit test
6. **Reduced Complexity**: Large monolithic file is broken into manageable pieces

## Backward Compatibility

The old `common_plots.py` file can be kept as a compatibility layer that imports from the new package structure, allowing for gradual migration of existing code.
