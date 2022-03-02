#!/usr/bin/env python
"""
pmultiqc plugin for MultiQC, showing the pipeline result with MultiQC framework.

For more information about MultiQC, see http://multiqc.info
"""

from setuptools import setup, find_packages

version = '0.0.10'


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='pmultiqc',
    version=version,
    description='Python package for quality control of proteomics datasets, based on multiqc package',
    author='Chengxin Dai, Yasset Perez-Riverol',
    author_email='S200502020@cqupt.edu.cn, ypriverol@gmail.com',
    long_description=readme(),
    long_description_content_type='text/markdown',
    keywords='Proteomics, Label-free, quality control, MultiQC',
    url='https://github.com/bigbio/pmultiqc/',
    download_url='https://github.com/bigbio/pmultiqc/',
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'lxml',
        'multiqc',
        'pandas',
        'pyteomics',
        'sdrf-pipelines >= 0.0.18',
        'numpy',
        'pyopenms'
    ],
    entry_points={
        'multiqc.modules.v1': [
            'quantms = pmultiqc.modules.quantms:QuantMSModule'
        ],
        'multiqc.cli_options.v1': [
            'disable_plugin = pmultiqc.cli:disable_plugin',
            'exp_design = pmultiqc.cli:exp_design',
            'pmultiqc_version = pmultiqc.cli:pmultiqc_version',
            'sdrf = pmultiqc.cli:sdrf',
            'raw = pmultiqc.cli:raw',
            'condition = pmultiqc.cli:condition',
            'mzMLs = pmultiqc.cli:mzMLs',
            'raw_ids = pmultiqc.cli:raw_ids',
            'remove_decoy = pmultiqc.cli:remove_decoy',
            'decoy_affix = pmultiqc.cli:decoy_affix',
            'affix_type = pmultiqc.cli:affix_type'
        ],
        'multiqc.hooks.v1': [
            'execution_start = pmultiqc.main:pmultiqc_plugin_execution_start'
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
