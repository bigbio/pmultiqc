#!/usr/bin/env python
"""
MultiQC command line options - we tie into the MultiQC
core here and add some new command line parameters.

See the Click documentation for more command line flag types:
http://click.pocoo.org/5/
"""

import click
from multiqc.utils import config


def print_version(ctx, params, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo("pmultiqc, version " + config.pmultiqc_version)
    ctx.exit()


pmultiqc_version = click.option(
    "--pmultiqc_version", is_flag=True, callback=print_version, expose_value=False, is_eager=True
)
raw = click.option(
    "--raw", "raw", help="Keep filenames in experimental design output as raw.", default=False
)
condition = click.option(
    "--condition", "condition", help="Create conditions from provided (e.g., factor) columns."
)
remove_decoy = click.option("--remove_decoy", "remove_decoy", default=True)
decoy_affix = click.option(
    "--decoy_affix",
    "decoy_affix",
    default="DECOY_",
    help="The decoy prefix or suffix used or to be used (default: DECOY_)",
)
contaminant_affix = click.option(
    "--contaminant_affix",
    "contaminant_affix",
    default="CONT",
    help="The contaminant prefix or suffix used or to be used (default: CONT_)",
)
quantification_method = click.option(
    "--quantification_method",
    "quantification_method",
    default="feature_intensity",
    help="The quantification method for LFQ experiment (default: feature_intensity)",
    type=click.Choice(["feature_intensity", "spectral_counting"]),
)
disable_table = click.option(
    "--disable_table",
    "disable_table",
    is_flag=True,
    help="disable protein/peptide table plots for large dataset",
)
ignored_idxml = click.option(
    "--ignored_idxml", "ignored_idxml", is_flag=True, help="ignored idxml files for faster running"
)
affix_type = click.option(
    "--affix_type", "affix_type", default="prefix", help="Prefix (default) or suffix"
)
disable_plugin = click.option(
    "--disable_plugin",
    "disable_plugin",
    is_flag=True,
    help="Disable the pmultiqc plugin on this run",
)
mzid_plugin = click.option("--mzid_plugin", "mzid_plugin", is_flag=True, help="Extract mzIdentML")
parse_maxquant = click.option(
    "--parse_maxquant", "parse_maxquant", is_flag=True, help="Parse MaxQuant results"
)
parse_proteobench = click.option(
    "--parse_proteobench", "parse_proteobench", is_flag=True, help="Parse ProteoBench results"
)
