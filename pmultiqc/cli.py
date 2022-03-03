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


exp_design = click.option('--exp_design', 'exp_design')
pmultiqc_version = click.option('--pmultiqc_version', is_flag=True, callback=print_version, expose_value=False,
                                is_eager=True)
sdrf = click.option('--sdrf', 'sdrf', help='The path of sdrf file')
raw = click.option('--raw', 'raw', help='Keep filenames in experimental design output as raw.', default=False)
condition = click.option('--condition', "condition",
                         help='Create conditions from provided (e.g., factor) columns.')
mzMLs = click.option('--mzMLs', 'mzMLs', default="mzMLs/", help='The directory of mzML files')
raw_ids = click.option('--raw_ids', 'raw_ids', default="raw_ids/", help="The directory of idXML files")
remove_decoy = click.option('--remove_decoy', 'remove_decoy', default=True)
decoy_affix = click.option('--decoy_affix', 'decoy_affix', default='DECOY_',
                           help='The decoy prefix or suffix used or to be used (default: DECOY_)')
affix_type = click.option('--affix_type', 'affix_type', default='prefix',
                          help='Prefix (default) or suffix')
disable_plugin = click.option('--disable_plugin', 'disable_plugin',
                              is_flag=True,
                              help="Disable the pmultiqc plugin on this run"
                              )
