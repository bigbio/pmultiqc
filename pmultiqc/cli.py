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
mzMLs = click.option('--mzMLs', 'mzMLs', required=True)
raw_ids = click.option('--raw_ids', 'raw_ids', required=True)
remove_decoy = click.option('--remove_decoy', 'remove_decoy', default=False)

disable_plugin = click.option('--disable-example-plugin', 'disable_plugin',
                              is_flag=True,
                              help="Disable the Example MultiQC plugin on this run"
                              )
