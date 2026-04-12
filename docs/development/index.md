# Development

Welcome to the pmultiqc development guide. This section covers everything you need to contribute to pmultiqc — from setting up your environment to adding new QC metrics and modules.

## Quick Links

| Topic | Description |
|-------|-------------|
| [Contributing](contributing.md) | Setup, code style, PR workflow |
| [Architecture](#architecture) | Module structure and data flow |
| [Adding Metrics](#adding-a-new-metric) | Step-by-step guide |
| [Adding Modules](#adding-a-new-module) | Support a new pipeline |

## Architecture

pmultiqc follows a modular architecture where each proteomics pipeline has its own module:

```
pmultiqc/
  cli.py                  ← CLI options (Click flags)
  main.py                 ← MultiQC plugin hook
  modules/
    core/                 ← Dispatcher, section utilities
    base.py               ← BasePMultiqcModule ABC
    common/
      common_utils.py     ← Shared parsing utilities
      dia_utils.py        ← DIA-specific parsing
      mzidentml_utils.py  ← mzIdentML parsing
      plots/
        id.py             ← Identification plots (PSM scores, peptide counts)
        ms.py             ← MS-level plots (mass error, RT, charge)
        general.py        ← Cross-cutting plots (heatmap, PCA, exp design)
    quantms/              ← quantms pipeline module
    maxquant/             ← MaxQuant module
    diann/                ← DIA-NN standalone module
    mzidentml/            ← mzIdentML module
    proteobench/          ← ProteoBench benchmarking module
    fragpipe/             ← FragPipe module
    mhcquant/             ← MHC immunopeptidomics module
```

### Data Flow

```
Input files (mzTab, evidence.txt, report.tsv, ...)
    │
    ▼
Module.get_data()          ← File discovery + parsing
    │
    ▼
Module.draw_plots()        ← Compute metrics, generate plots
    │
    ▼
common/plots/id.py         ← Identification plot functions
common/plots/ms.py         ← MS-level plot functions
common/plots/general.py    ← General plot functions
    │
    ▼
MultiQC HTML report        ← Assembled by MultiQC framework
```

### Module Interface

Every module extends `BasePMultiqcModule`:

```python
class MyToolModule(BasePMultiqcModule):
    def get_data(self):
        """Discover and parse input files.
        Returns True if data was found."""
        ...

    def draw_plots(self):
        """Create QC visualizations.
        Calls plot functions from common/plots/."""
        ...
```

## Adding a New Metric

1. **Implement the calculation** in the appropriate utility file:
   - `common/common_utils.py` for shared metrics
   - Module-specific `*_utils.py` for tool-specific metrics

2. **Create the plot function** in `common/plots/`:

```python
# common/plots/id.py
def draw_my_metric(sub_section, data_dict):
    pconfig = {
        "id": "my_metric_id",
        "title": "My Metric",
        "ylab": "Count",
        "xlab": "Sample",
    }
    html = bargraph.plot(data=data_dict, pconfig=pconfig)
    html = plot_html_check(html)
    add_sub_section(
        sub_section=sub_section,
        plot=html,
        order=10,
        description="Description of what this metric shows.",
        helptext="Detailed help text for users.",
    )
```

3. **Call the function** from the module's `draw_plots()` method.

4. **Add a test** in `tests/` that verifies the metric with minimal input data.

5. **Update documentation** — add the metric to `docs/reference/qc-metrics.md`.

## Adding a New Module

To support a new proteomics pipeline or file format:

1. **Create the module directory**:

```
pmultiqc/modules/your_tool/
  __init__.py
  your_tool.py
```

2. **Implement the module**:

```python
# your_tool.py
from pmultiqc.modules.base import BasePMultiqcModule

class YourToolModule(BasePMultiqcModule):
    def get_data(self):
        # Find and parse your tool's output files
        files = self.find_files("your_pattern*.tsv")
        if not files:
            return False
        # Parse files into internal data structures
        self.parse_files(files)
        return True

    def draw_plots(self):
        # Generate QC plots using common plot functions
        from pmultiqc.modules.common.plots import id, ms
        id.draw_psm_table(self.sub_sections, self.psm_data)
        ms.draw_mass_error(self.sub_sections, self.ms_data)
```

3. **Register the module** in `pmultiqc/modules/core/__init__.py`.

4. **Add a CLI flag** in `pmultiqc/cli.py`:

```python
your_tool_plugin = click.option(
    "--your-tool-plugin",
    "your_tool_plugin",
    is_flag=True,
    help="Enable YourTool plugin",
)
```

5. **Register the entry point** in `pyproject.toml`:

```toml
[project.entry-points."multiqc.cli_options.v1"]
your_tool_plugin = "pmultiqc.cli:your_tool_plugin"
```

6. **Add test data** in `tests/resources/your_tool/` and write tests.

## Running Tests

```bash
# Full test suite
pytest tests/

# Test with specific data
cd tests
multiqc resources/LFQ -o ./test_output --quantms-plugin
multiqc resources/maxquant -o ./test_output --maxquant-plugin
```

## Code Style

```bash
# Format code
black pmultiqc/ --line-length 99

# Sort imports
isort pmultiqc/

# Run linting
ruff check pmultiqc/
```

## Submitting a Pull Request

1. Fork and create a feature branch: `git checkout -b feature/my-metric`
2. Make changes, format with `black` and `isort`
3. Run tests: `pytest tests/`
4. Push and open a PR against `dev`
5. Fill in the PR template

### PR Checklist

- [ ] Code formatted with black/isort
- [ ] Tests added or updated
- [ ] Documentation updated (qc-metrics.md or user guide)
- [ ] PR description explains what changed and why
