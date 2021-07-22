#!/usr/bin/env python
"""
MultiQC command line options - we tie into the MultiQC
core here and add some new command line parameters.

See the Click documentation for more command line flag types:
http://click.pocoo.org/5/
"""

import click

exp_design = click.option('--exp_design', 'exp_design')
sdrf = click.option('--sdrf', 'sdrf')
raw = click.option('--raw', 'raw', help='Keep filenames in experimental design output as raw.', default=False)
condition = click.option('--condition', "condition",
                         help='Create conditions from provided (e.g., factor) columns.')
quant_method = click.option('--quant_method', 'quant_method', help='lfq or tmt', default='lfq')
mzMLs = click.option('--mzMLs', 'mzMLs', required=True)
raw_ids = click.option('--raw_ids', 'raw_ids', required=True)
remove_decoy = click.option('--remove_decoy', 'remove_decoy', default=True)
decoy_affix = click.option('--decoy_affix', 'decoy_affix', default='DECOY_',
                           help='The decoy prefix or suffix used or to be used (default: DECOY_)')
affix_type = click.option('--affix_type', 'affix_type', default='prefix',
                          help='Prefix (default) or suffix')
disable_plugin = click.option('--disable-example-plugin', 'disable_plugin',
                              is_flag=True,
                              help="Disable the Example MultiQC plugin on this run"
                              )
