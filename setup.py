#!/usr/bin/env python
"""
ProteomicsLFQ plugin for MultiQC, showing the pipeline result with MultiQC framework.

For more information about MultiQC, see http://multiqc.info
"""

from setuptools import setup, find_packages

version = '0.0.2'

setup(
    name='pmultiqc',
    version=version,
    author='Chengxin Dai',
    author_email='S200502020@cqupt.edu.cn',
    description="proteomics label free analysis pipeline quality control report",
    long_description=__doc__,
    keywords='Proteomics',
    url='https://github.com/bigbio/pmultiqc/',
    download_url='https://github.com/bigbio/pmultiqc/',
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'multiqc',
        'pandas',
        'pyteomics',
        'sdrf_pipelines'
    ],
    entry_points={
        'multiqc.modules.v1': [
            'proteomicslfq = pmultiqc.modules.proteomicslfq:MultiqcModule'
        ],
        'multiqc.cli_options.v1': [
            'disable_plugin = pmultiqc.cli:disable_plugin',
            'exp_design = pmultiqc.cli:exp_design',
            'sdrf = pmultiqc.cli:sdrf',
            'raw = pmultiqc.cli:raw',
            'condition = pmultiqc.cli:condition',
            'mzMLs = pmultiqc.cli:mzMLs',
            'raw_ids = pmultiqc.cli:raw_ids',
            'remove_decoy = pmultiqc.cli:remove_decoy'
        ],
        'multiqc.hooks.v1': [
            'execution_start = pmultiqc.custom_code:example_plugin_execution_start'
        ]
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)
