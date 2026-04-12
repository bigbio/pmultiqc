# Contributing to pmultiqc

pmultiqc is an open-source MultiQC plugin for proteomics QC reporting. Contributions are welcome — bug fixes, new metrics, new module support, and documentation improvements.

## Setting Up the Development Environment

```bash
# Clone the repository
git clone https://github.com/bigbio/pmultiqc.git
cd pmultiqc

# Create a conda environment (recommended)
conda env create -f environment.yml
conda activate pmultiqc

# Install in editable mode with development dependencies
pip install -e .
pip install pytest black isort
```

Requirements: Python 3.10–3.13, MultiQC 1.29–1.33.

## Running Tests

```bash
# Run the full test suite
pytest tests/

# Test with sample LFQ data
cd tests && multiqc resources/LFQ -o ./test_output --quantms-plugin

# Test with MaxQuant sample data
cd tests && multiqc resources/maxquant -o ./test_output --maxquant-plugin
```

Tests are in `tests/`. Test resources (small representative datasets) are in `tests/resources/`.

## Code Style

- Format with `black` (line length 99): `black pmultiqc/`
- Sort imports with `isort`: `isort pmultiqc/`
- Follow PEP 8; use type hints where practical
- Write docstrings for all public functions and classes

## Adding a New Metric to an Existing Module

1. **Implement the calculation** in the appropriate utility file (e.g., `pmultiqc/modules/common/common_utils.py` for shared utilities, or the module-specific `*_utils.py`).

2. **Add a plot function** in the corresponding plots file under `pmultiqc/modules/common/plots/` (`id.py` for identification metrics, `ms.py` for spectral metrics, `general.py` for cross-cutting plots).

   Plot functions should follow the pattern:
   ```python
   def draw_my_metric(sub_section, data_dict):
       pconfig = {"id": "my_metric_id", "title": "My Metric", ...}
       html = bargraph.plot(data=data_dict, pconfig=pconfig)
       html = plot_html_check(html)
       add_sub_section(sub_section=sub_section, plot=html, order=N, description="...", helptext="...")
   ```

3. **Call the function** from the module's `draw_plots()` method.

4. **Add a test** that verifies the metric is computed correctly with a minimal input fixture.

5. **Update `docs/reference/qc-metrics.md`** with the new metric entry.

## Adding a New Module

New modules (supporting a new pipeline or file format) follow the `BasePMultiqcModule` pattern:

1. Create a directory under `pmultiqc/modules/your_tool/` with `__init__.py` and `your_tool.py`.

2. Subclass `BasePMultiqcModule` and implement:
   - `get_data()` — file discovery and parsing; returns `True` if data was found
   - `draw_plots()` — calls plot functions and populates `sub_sections`

3. Register the module in `pmultiqc/modules/core/__init__.py` in the `PMultiQC` dispatcher.

4. Add a CLI flag in `pmultiqc/cli.py`:
   ```python
   your_tool_plugin = click.option(
       "--your-tool-plugin", "your_tool_plugin", is_flag=True, help="Enable YourTool plugin"
   )
   ```
   Register it in `pyproject.toml` under `[project.entry-points."multiqc.cli_options.v1"]`.

5. Add test data in `tests/resources/your_tool/` and a test function in `tests/`.

## Submitting a Pull Request

1. Create a feature branch: `git checkout -b feature/my-metric`
2. Make changes, run `black` and `isort`, run `pytest`
3. Push: `git push origin feature/my-metric`
4. Open a PR against `main` on GitHub
5. Fill in the PR template: what changed, why, how to test it

PR checklist:
- [ ] Code formatted with black/isort
- [ ] Tests added or updated and passing
- [ ] Documentation updated (qc-metrics.md or relevant user guide)
- [ ] PR description explains the change

## Reporting Issues

Use the GitHub issue tracker with the appropriate template:

- **Bug Report** — crashes, wrong metric values, unexpected behavior
- **Metric Request** — new proteomics QC metric (encouraged)
- **Feature Request** — new visualization type, file format support

For urgent issues: ypriverol@gmail.com

## Architecture Overview

```
pmultiqc/
  cli.py                  # Click CLI options registered as MultiQC entry points
  main.py                 # Plugin initialization hook
  modules/
    core/                 # PMultiQC dispatcher, section group utilities
    base.py               # BasePMultiqcModule ABC
    common/
      common_utils.py     # Shared parsing and computation utilities
      dia_utils.py        # DIA-specific utilities (DIA-NN report parsing)
      mzidentml_utils.py  # mzIdentML-specific utilities
      plots/
        id.py             # Identification-level plot functions
        ms.py             # MS-level plot functions
        general.py        # Cross-cutting plot functions (heatmap, exp design)
    quantms/              # quantms pipeline module
    maxquant/             # MaxQuant module
    diann/                # DIA-NN standalone module
    mzidentml/            # mzIdentML module
    proteobench/          # ProteoBench module
    fragpipe/             # FragPipe module
    mhcquant/             # nf-core/mhcquant module
```
