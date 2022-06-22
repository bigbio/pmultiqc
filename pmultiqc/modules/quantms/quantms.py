#!/usr/bin/env python

""" MultiQC pmultiqc plugin module """

from __future__ import absolute_import
from collections import OrderedDict
import itertools
from datetime import datetime
from operator import itemgetter
import logging
from sdrf_pipelines.openms.openms import OpenMS, UnimodDatabase
from multiqc import config
from multiqc.plots import table, bargraph, linegraph, heatmap
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.utils.mqc_colour import mqc_colour_scale
import pandas as pd
import math
from functools import reduce
from collections import Counter
import re
from pyteomics import mztab
from pyopenms import IdXMLFile, MzMLFile, MSExperiment, OpenMSBuildInfo
import os
import sqlite3
import numpy as np
import math
import copy
import json

from .histogram import Histogram
from . import sparklines

# Initialise the main MultiQC logger
log = logging.getLogger(__name__)

if config.output_dir:
    if os.path.exists(config.output_dir):
        con = sqlite3.connect(os.path.join(config.output_dir, 'quantms.db'))
    else:
        os.makedirs(config.output_dir)
        con = sqlite3.connect(os.path.join(config.output_dir, 'quantms.db'))
else:
    con = sqlite3.connect('./quantms.db')

cur = con.cursor()
cur.execute("drop table if exists PROTQUANT")
con.commit()

cur.execute("drop table if exists PEPQUANT")
con.commit()


log.warning("pyopenms has: "+str(OpenMSBuildInfo().getOpenMPMaxNumThreads())+" threads.")

class QuantMSModule(BaseMultiqcModule):
    def __init__(self):

        # Initialise the parent module Class object
        super(QuantMSModule, self).__init__(
            name='pmultiqc',
            target="pmultiqc",
            anchor='pmultiqc',
            href='https://github.com/bigbio/pmultiqc',
            info=" is a multiQC module to show the pipeline performance of mass spectrometry based quantification pipelines such as <a href='https://nf-co.re/quantms'>nf-core/quantms</a>."
        )

        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_plugin', True):
            return None

        self.enable_exp = False
        self.enable_sdrf = False
        msstats_input_valid = False
        # TODO what if multiple are found??
        for f in self.find_log_files("quantms/exp_design", filecontents=False):
            self.exp_design = os.path.join(f["root"], f['fn'])
            self.enable_exp = True
        if self.enable_exp == False:
            for f in self.find_log_files("quantms/sdrf", filecontents=False):
                self.sdrf = os.path.join(f["root"], f['fn'])
                OpenMS().openms_convert(self.sdrf, config.kwargs['raw'],
                                    False, True, False, config.kwargs['condition'])
                # experimental_design.tsv is the default output name
                self.exp_design = os.path.join(f["root"], f['experimental_design.tsv'])
                self.enable_sdrf = True

        # TODO in theory the module can work without the design. We just need to remove the according sections!
        if self.enable_sdrf == False and self.enable_exp == False:
            raise AttributeError("Neither exp_design nor sdrf can be found! Please provide or correct your multiqc_config.yml.")

        # TODO what if multiple are found??
        for msstats_input in self.find_log_files("quantms/msstats", filecontents=False):
            self.msstats_input_path = os.path.join(msstats_input["root"], msstats_input["fn"])
            msstats_input_valid = True
        if msstats_input_valid == False:
            raise AttributeError("MSstats input file not found! Please provide or correct your multiqc_config.yml.")

        self.PSM_table = dict()
        self.mzml_peptide_map = dict()
        self.identified_spectrum = dict()
        self.pep_quant_table = dict()
        self.mzml_table = OrderedDict()
        self.total_ms2_spectra = 0
        self.Total_ms2_Spectral_Identified = 0
        self.Total_Peptide_Count = 0
        self.Total_Protein_Identified = 0
        self.Total_Protein_Quantified = 0
        self.out_csv_data = dict()
        self.cal_num_table_data = dict()
        self.mL_spec_ident_final = dict()
        self.delta_mass = dict()
        self.delta_mass_percent = dict()
        self.heatmap_con_score = dict()
        self.heatmap_pep_intensity = {}
        self.heatmap_charge_score = dict()
        self.MissedCleavages_heatmap_score = dict()
        self.MissedCleavagesVar_score = dict()
        self.ID_RT_score = dict()
        self.HeatmapOverSamplingScore = dict()
        self.HeatmapPepMissingScore = dict()
        self.exp_design_table = dict()
        self.mzml_info = dict()
        self.mzml_info['charge_distribution'] = dict()
        self.mzml_info['peaks_per_ms2'] = dict()
        self.mzml_info['peak_distribution'] = dict()
        self.oversampling = dict()

        # parse input data
        # draw the experimental design
        self.draw_exp_design()
        self.pep_table_exists = False            
        self.enable_dia = False
        for f in self.find_log_files("quantms/mztab", filecontents=False):
            self.out_mzTab_path = os.path.join(f["root"], f['fn'])
            self.parse_out_mzTab()

        for report in self.find_log_files("quantms/diann_report", filecontents=False):
            self.diann_report_path = os.path.join(report["root"], report["fn"])
            self.enable_dia = True

        self.mzML_paths = []
        for mzML_file in self.find_log_files("quantms/mzML", filecontents=False):
            self.mzML_paths.append(os.path.join(mzML_file["root"], mzML_file['fn']))

        mt = self.parse_mzml()
        self.idx_paths = []
        for idx_file in self.find_log_files("quantms/idXML", filecontents=False):
            self.idx_paths.append(os.path.join(idx_file["root"], idx_file['fn']))

        if self.enable_dia:
            self.parse_diann_report()
            self.draw_summary_protein_ident_table()
            self.draw_quantms_identi_num()
            self.draw_num_pep_per_protein()
            self.draw_precursor_charge_distribution()
            self.draw_peaks_per_ms2()
            self.draw_peak_intensity_distribution()
        else:
            self.parse_idxml(mt)
            self.CalHeatMapScore()
            self.draw_heatmap()
            self.draw_summary_protein_ident_table()
            self.draw_quantms_identi_num()
            self.draw_num_pep_per_protein()
            self.draw_mzml_ms()
            self.draw_precursor_charge_distribution()
            self.draw_peaks_per_ms2()
            self.draw_peak_intensity_distribution()
            self.draw_oversampling()
            self.draw_delta_mass()

        self.parse_msstats_input()
        self.css = {
            'assets/css/quantms.css':
                os.path.join(os.path.dirname(__file__), 'assets', 'css', 'quantms.css')
        }
        self.js = {
            'assets/js/quantms.js':
                os.path.join(os.path.dirname(__file__), 'assets', 'js', 'quantms.js'),
            'assets/js/highcharts.js':
                os.path.join(os.path.dirname(__file__), 'assets', 'js', 'highcharts.js'),
            'assets/js/axios.min.js':
                os.path.join(os.path.dirname(__file__), 'assets', 'js', 'axios.min.js'),
            'assets/js/sql-optimized.js':
                os.path.join(os.path.dirname(__file__), 'assets', 'js', 'sql-optimized.js')
        }

    def draw_heatmap(self):
        HeatMapScore = []
        if self.pep_table_exists:
            xnames = ['Contaminants', 'Peptide Intensity', 'Charge', 'Missed Cleavages', 'Missed Cleavages Var',
                        'ID rate over RT', 'MS2 OverSampling', 'Pep Missing Values']
            ynames = []
            for k, v in self.heatmap_con_score.items():
                ynames.append(k)
                HeatMapScore.append([v, self.heatmap_pep_intensity[k], self.heatmap_charge_score[k],
                                    self.MissedCleavages_heatmap_score[k],
                                    self.MissedCleavagesVar_score[k], self.ID_RT_score[k],
                                    self.HeatmapOverSamplingScore[k], self.HeatmapPepMissingScore[k]])
        else:
            xnames = ['Charge', 'Missed Cleavages', 'Missed Cleavages Var', 'ID rate over RT', 'MS2 OverSampling', 
                        'Pep Missing Values']
            ynames = []
            for k, v in self.heatmap_charge_score.items():
                ynames.append(k)
                HeatMapScore.append([self.heatmap_charge_score[k], self.MissedCleavages_heatmap_score[k],
                                    self.MissedCleavagesVar_score[k], self.ID_RT_score[k],
                                    self.HeatmapOverSamplingScore[k], self.HeatmapPepMissingScore[k]])

        pconfig = {
            'title': 'Performance Overview',  # Plot title - should be in format "Module Name: Plot Title"
            'xTitle': 'metrics',  # X-axis title
            'yTitle': 'RawName',  # Y-axis title
            'square': False
        }

        hm_html = heatmap.plot(HeatMapScore, xnames, ynames, pconfig)
        # Add a report section with the heatmap plot
        self.add_section(
            name="HeatMap",
            anchor="quantms_heatmap",
            description='This heatmap shows a performance overview of the pipeline',
            helptext='''
                    This plot shows the pipeline performance overview. Some metrics are calculated.
                    
                    * Heatmap score[Contaminants]: as fraction of summed intensity with 0 = sample full of contaminants; 
                        1 = no contaminants
                    * Heatmap score[Pep Intensity (>23.0)]: Linear scale of the median intensity reaching the threshold, 
                        i.e. reaching 2^21 of 2^23 gives score 0.25.
                    * Heatmap score[Charge]: Deviation of the charge 2 proportion from a representative 
                        Raw file (median). For typtic digests, peptides of charge 2 (one N-terminal and one at 
                        tryptic C-terminal R or K residue) should be dominant. Ionization issues (voltage?), 
                        in-source fragmentation, missed cleavages and buffer irregularities can cause a shift 
                        (see Bittremieux 2017, DOI: 10.1002/mas.21544 ).
                    * Heatmap score [MC]: the fraction (0% - 100%) of fully cleaved peptides per Raw file
                    * Heatmap score [MC Var]: each Raw file is scored for its deviation from the ‘average’ digestion 
                        state of the current study.
                    * Heatmap score [ID rate over RT]: Judge column occupancy over retention time. 
                        Ideally, the LC gradient is chosen such that the number of identifications 
                        (here, after FDR filtering) is uniform over time, to ensure consistent instrument duty cycles. 
                        Sharp peaks and uneven distribution of identifications over time indicate potential for LC gradient 
                        optimization.Scored using ‘Uniform’ scoring function. i.e. constant receives good score, extreme shapes are bad.
                    * Heatmap score [MS2 Oversampling]: The percentage of non-oversampled 3D-peaks. An oversampled 
                        3D-peak is defined as a peak whose peptide ion (same sequence and same charge state) was 
                        identified by at least two distinct MS2 spectra in the same Raw file. For high complexity samples, 
                        oversampling of individual 3D-peaks automatically leads to undersampling or even omission of other 3D-peaks, 
                        reducing the number of identified peptides.
                    * Heatmap score [Pep Missing]: Linear scale of the fraction of missing peptides.
                    
                    ''',
            plot=hm_html
        )

    def draw_exp_design(self):
        # Currently this only supports the OpenMS two-table format (default in quantms pipeline)
        # One table format would actually be even easier. You can just use pandas.read_tsv
        self.sample_df, self.file_df = read_openms_design(self.exp_design)

        for file in np.unique(self.file_df["Spectra_Filepath"].tolist()):
            stand_file = os.path.basename(file)
            file_index = self.file_df[self.file_df["Spectra_Filepath"] == file]
            self.exp_design_table[stand_file] = {'Fraction_Group': file_index["Fraction_Group"].tolist()[0]}
            self.exp_design_table[stand_file]['Fraction'] = file_index["Fraction"].tolist()[0]
            self.exp_design_table[stand_file]['Label'] = '|'.join(file_index["Label"])
            sample = file_index["Sample"].tolist()
            self.exp_design_table[stand_file]['Sample'] = '|'.join(sample)
            self.exp_design_table[stand_file]['MSstats_Condition'] = ','.join(
                [row['MSstats_Condition'] for _, row in self.sample_df.iterrows()
                if row['Sample'] in sample])
            self.exp_design_table[stand_file]['MSstats_BioReplicate'] = '|'.join(
                [row['MSstats_BioReplicate'] for _, row in self.sample_df.iterrows()
                if row['Sample'] in sample])

        # Create table plot
        pconfig = {
            'id': 'Experimental_Design',  # ID used for the table
            'table_title': 'Experimental Design',  # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'raw_data_fn': 'multiqc_Experimental_Design_table',  # File basename to use for raw data file
            'sortRows': False,  # Whether to sort rows alphabetically
            'only_defined_headers': False,  # Only show columns that are defined in the headers config
            'col1_header': 'Spectra File',
            'format': '{:,.0f}',  # The header used for the first column
        }
        headers = OrderedDict()
        set3_scale = mqc_colour_scale(name="Set3")
        maxnr = len(self.file_df.index)
        set3_colors = set3_scale.get_colours(name="Set3")
        colors = dict( (str(i+1), set3_colors[i % len(set3_colors)]) for i in range(maxnr) )

        headers['Fraction_Group'] = {
            'description': 'Fraction_Group',
            'bgcols': colors,
        }
        headers['Fraction'] = {
            'description': 'Fraction identifier',
            'bgcols': colors,
        }
        headers['Label'] = {
            'description': 'Label',
            'bgcols': colors,
        }
        headers['Sample'] = {
            'description': 'Sample Name',
            'bgcols': colors,
        }
        headers['MSstats_Condition'] = {
            'description': 'MSstats Condition',
            'bgcols': colors,
        }
        headers['MSstats_BioReplicate'] = {
            'description': 'MSstats BioReplicate',
            'bgcols': colors,
        }
        table_html = table.plot(self.exp_design_table, headers, pconfig)

        # Add a report section with the line plot
        self.add_section(
            name="Experimental Design",
            anchor="exp_design",
            description='This table shows the design of the experiment. I.e., which files and channels correspond to which sample/condition/fraction.',
            helptext='''
            You can see details about it in 
            https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/classOpenMS_1_1ExperimentalDesign.html
            ''',
            plot=table_html
        )

    def draw_summary_protein_ident_table(self):
        headers = OrderedDict()
        if self.enable_dia:
            summary_table = {
                self.total_ms2_spectra: {"#Peptides Quantified": self.Total_Peptide_Count}}
            summary_table[self.total_ms2_spectra]['#Proteins Quantified'] = self.Total_Protein_Quantified    
        else:
            summary_table = {
                self.total_ms2_spectra: {"#Identified MS2 Spectra": self.Total_ms2_Spectral_Identified}}
            coverage = self.Total_ms2_Spectral_Identified / self.total_ms2_spectra * 100
            summary_table[self.total_ms2_spectra]['%Identified MS2 Spectra'] = coverage
            summary_table[self.total_ms2_spectra]['#Peptides Identified'] = self.Total_Peptide_Count
            summary_table[self.total_ms2_spectra]['#Proteins Identified'] = self.Total_Protein_Identified
            summary_table[self.total_ms2_spectra]['#Proteins Quantified'] = self.Total_Protein_Quantified
            headers['#Identified MS2 Spectra'] = {
                'description': 'Total number of MS/MS spectra identified',
            }
            headers['%Identified MS2 Spectra'] = {
                'description': 'Percentage of Identified MS/MS Spectra',
                'format': '{:,.2f}',
                "suffix": "%"
            }
        # Create table plot
        pconfig = {
            'id': 'identification summary table',  # ID used for the table
            'table_title': 'Summary Table',  # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'raw_data_fn': 'multiqc_summary_table_table',  # File basename to use for raw data file
            'sortRows': False,  # Whether to sort rows alphabetically
            'only_defined_headers': False,  # Only show columns that are defined in the headers config
            'col1_header': '#MS2 Spectra',
            'format': '{:,.0f}',  # The header used for the first column
            "scale": "Set1"
        }
        table_html = table.plot(summary_table, headers, pconfig)

        # Add a report section with the table
        self.add_section(
            name="Summary Table",
            anchor="quantms_summary_table",
            description='This table shows the quantms pipeline summary statistics',
            helptext='''
                    This table shows the quantms pipeline summary statistics
                    ''',
            plot=table_html
        )

    def draw_quantms_identi_num(self):
        # Create table plot
        pconfig = {
            'id': 'result statistics',  # ID used for the table
            'table_title': 'pipeline result statistics',  # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'raw_data_fn': 'multiqc_result statistics_table',  # File basename to use for raw data file
            'sortRows': True,  # Whether to sort rows alphabetically
            'only_defined_headers': False,  # Only show columns that are defined in the headers config
            'col1_header': 'Spectra File',
            'format': '{:,.0f}'  # The header used for the first column
        }
        headers = OrderedDict()

        ## TODO BIG! We desperately need rules for naming columns!

        headers['Sample Name'] = {
            'title': 'Sample Name',
            'description': 'Sample identifier',
            'color': "#ffffff"
        }
        headers['Condition'] = {
                'title': 'Condition',
                'description': 'Possible Study Variables',
                'color': "#ffffff"
        }
        headers['Fraction'] = {
            'title': 'Fraction',
            'description': 'Fraction identifier',
            'color': "#ffffff"
        }
        headers['peptide_num'] = {
            'title': '#Peptide IDs',
            'description': 'The number of identified PSMs in the pipeline',
            'color': "#ffffff"
        }
        headers['unique_peptide_num'] = {
            'title': '#Unambiguous Peptide IDs',
            'description': 'The number of unique peptides in the pipeline. Those that match only one protein in the provided database.',
            'color': "#ffffff"
        }
        headers['modified_peptide_num'] = {
            'title': '#Modified Peptide IDs',
            'description': 'Number of modified identified peptides in the pipeline',
            'color': "#ffffff"
        }
        headers['protein_num'] = {
            'title': '#Protein (group) IDs',
            'description': 'The number of identified protein(group)s in the pipeline',
            'color': "#ffffff"
        }
        table_html = table.plot(self.cal_num_table_data, headers, pconfig)

        # Add a report section with the line plot
        self.add_section(
            name="Pipeline Result Statistics",
            anchor="Pipeline Result Statistics",
            description='This plot shows the quantms pipeline final result',
            helptext='''
            This plot shows the quantms pipeline final result.
            Including Sample Name、Possible Study Variables、identified the number of peptide in the pipeline、
            and identified the number of modified peptide in the pipeline, eg. All data in this table are obtained 
            from the out_msstats file. You can also remove the decoy with the `remove_decoy` parameter.
            ''',
            plot=table_html
        )

    # draw number of peptides per proteins
    def draw_num_pep_per_protein(self):
        # Create table plot
        pconfig = {
            'id': 'number_of_peptides_per_proteins',  # ID used for the table
            'cpswitch': False,
            'title': 'Number of Peptides identified per Protein',
            'xlab': 'Number of Peptides',
            'tt_percentages': False, 
            'tt_decimals': 2,
            'data_labels': [
                {'name': 'count', 'ylab': 'Frequency'},
                {'name': 'percentage', 'ylab': 'Percentage(%)'}
            ]
        }
        headers = OrderedDict()
        headers['Frequency'] = {
            'name': 'Frequency',
            'description': 'number of peptides per proteins'
        }
        bar_html = bargraph.plot([self.pep_plot.dict['data']['frequency'], self.pep_plot.dict['data']['percentage']], headers, pconfig)
        # Add a report section with the line plot
        self.add_section(
            name="Number of Peptides Per Protein",
            anchor="num_of_pep_per_prot",
            description='This plot shows the number of peptides per proteins '
                        'in quantms pipeline final result',
            helptext='''
                        This statistic is extracted from the out_msstats file. Proteins supported by more peptide 
                        identifications can constitute more confident results.
                    ''',
            plot=bar_html
        )

    def draw_mzml_ms(self):

        pconfig = {
            'id': 'mzML_tracking ',  # ID used for the table
            'table_title': 'Pipeline spectrum tracking',  # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'raw_data_fn': 'multiqc_spectrum_tracking_table',  # File basename to use for raw data file
            'sortRows': False,  # Whether to sort rows alphabetically
            'only_defined_headers': True,  # Only show columns that are defined in the headers config
            'col1_header': 'Spectra File',
            'format': '{:,.0f}'  # The header used for the first column
        }

        headers = OrderedDict()
        
        headers['MS1_Num'] = {
            'title': '#MS1 Spectra',
            'description': 'Number of MS1 spectra',
            'color': "#ffffff"
        }
        headers['MS2_Num'] = {
            'title': '#MS2 Spectra',
            'description': 'Number of MS2 spectra',
            'color': "#ffffff"
        }

        if any(['MSGF' in v for k,v in self.mzml_table.items()]):
            headers['MSGF'] = {
                'description': 'Number of spectra identified by MSGF search engine',
                'color': "#ffffff"
            }
        if any(['Comet' in v for k,v in self.mzml_table.items()]):
            headers['Comet'] = {
                'description': 'Number of spectra identified by Comet search engine',
                'color': "#ffffff"
            }
        headers['num_quant_psms'] = {
            'title': '#PSMs from quant. peptides',
            'description': 'Number of reliable PSMs from peptides IDs used in quantification',
            'color': "#ffffff"
        }
        headers['num_quant_peps'] = {
            'title': '#Peptides quantified',
            'description': 'Number of quantified peptides that passed final protein and peptide FDR thresholds.',
            'color': "#ffffff"
        }
        table_html = table.plot(self.mzml_table, headers, pconfig)

        # Add a report section with the line plot
        self.add_section(
            name="Spectra Tracking",
            anchor="spectra_tracking",
            description='This plot shows the tracking of the number of spectra along the quantms pipeline',
            helptext='''
                    This table shows the changes in the number of spectra corresponding to each input file 
                    during the pipeline operation. And the number of peptides finally identified and quantified is obtained from 
                    the PSM table in the mzTab file. You can also remove decoys with the `remove_decoy` parameter.:

                    * MS1_Num: The number of MS1 spectra extracted from mzMLs
                    * MS2_Num: The number of MS2 spectra extracted from mzMLs
                    * MSGF: The Number of spectra identified by MSGF search engine
                    * Comet: The Number of spectra identified by Comet search engine
                    * #PSMs from quant. peptides: extracted from PSM table in mzTab file
                    * #Peptides quantified: extracted from PSM table in mzTab file
                    ''',
            plot=table_html
        )


    def draw_peak_intensity_distribution(self):
        # Create table plot
        pconfig = {
            'id': 'peak_intensity_distribution',  # ID used for the table
            'cpswitch': False,
            'title': 'Peak Intensity Distribution',
            # 'xlab': 'Peak Intensity'
        }
        cats = self.mzml_peak_distribution_plot.dict['cats']
        bar_html = bargraph.plot(self.mzml_info['peak_distribution'], cats, pconfig)

        # Add a report section with the line plot
        self.add_section(
            name="Peak Intensity Distribution",
            anchor="Peak Intensity Distribution",
            description='''This is a histogram representing the ion intensity vs.
                        the frequency for all MS2 spectra in a whole given experiment.
                        It is possible to filter the information for all, identified and unidentified spectra.
                        This plot can give a general estimation of the noise level of the spectra.
                        ''',
            helptext='''
                    Generally, one should expect to have a high number of low intensity noise peaks with a low number 
                    of high intensity signal peaks. 
                    A disproportionate number of high signal peaks may indicate heavy spectrum pre-filtering or 
                    potential experimental problems. In the case of data reuse this plot can be useful in 
                    identifying the requirement for pre-processing of the spectra prior to any downstream analysis. 
                    The quality of the identifications is not linked to this data as most search engines perform internal 
                    spectrum pre-processing before matching the spectra. Thus, the spectra reported are not 
                    necessarily pre-processed since the search engine may have applied the pre-processing step 
                    internally. This pre-processing is not necessarily reported in the experimental metadata.
                    ''',
            plot=bar_html
        )

    def draw_precursor_charge_distribution(self):
        # Create table plot
        pconfig = {
            'id': 'Precursor_Ion_Charge_Distribution',  # ID used for the table
            'cpswitch': False,
            'title': 'Precursor Ion Charge Distribution',
        }
        cats = self.mzml_charge_plot.dict['cats']
        bar_html = bargraph.plot(self.mzml_info['charge_distribution'], cats, pconfig)

        # Add a report section with the line plot
        self.add_section(
            name="Distribution of precursor charges",
            anchor="Distribution of precursor charges",
            description='''This is a bar chart representing the distribution of the precursor ion charges 
                        for a given whole experiment. 
                        ''',
            helptext='''This information can be used to identify potential ionization problems 
                        including many 1+ charges from an ESI ionization source or an unexpected 
                        distribution of charges. MALDI experiments are expected to contain almost exclusively 1+ 
                        charged ions. An unexpected charge distribution may furthermore be caused by specific search 
                        engine parameter settings such as limiting the search to specific ion charges.
                    ''',
            plot=bar_html
        )

    def draw_peaks_per_ms2(self):
        # Create table plot
        pconfig = {
            'id': 'peaks_per_ms2',  # ID used for the table
            'cpswitch': False,
            'title': 'Number of Peaks per MS/MS spectrum',
        }
        cats = self.mzml_peaks_ms2_plot.dict['cats']
        bar_html = bargraph.plot(self.mzml_info['peaks_per_ms2'], cats, pconfig)

        self.add_section(
            name="Number of Peaks per MS/MS spectrum",
            anchor="Number of Peaks per MS/MS spectrum",
            description='''This chart represents a histogram containing the number of peaks per MS/MS spectrum 
            in a given experiment. This chart assumes centroid data. Too few peaks can identify poor fragmentation 
            or a detector fault, as opposed to a large number of peaks representing very noisy spectra. 
            This chart is extensively dependent on the pre-processing steps performed to the spectra 
            (centroiding, deconvolution, peak picking approach, etc).
                        ''',
            helptext='''
                    ''',
            plot=bar_html
        )

    def draw_oversampling(self):

        # Create bar plot
        pconfig = {
            'id': 'Oversampling_Distribution',
            'cpswitch': False,
            'title': 'MS2 counts per 3D-peak',
            'scale': "set3"
        }
        cats = self.oversampling_plot.dict['cats']
        bar_html = bargraph.plot(self.oversampling, cats, pconfig)

        # Add a report section with the bar plot
        self.add_section(
            name="Oversampling Distribution",
            anchor="Oversampling (MS/MS counts per 3D-peak)",
            description='''An oversampled 3D-peak is defined as a peak whose peptide ion 
                (same sequence and same charge state) was identified by at least two distinct MS2 spectra 
                in the same Raw file. 
                                ''',
            helptext='''For high complexity samples, oversampling of individual 3D-peaks automatically leads to 
                undersampling or even omission of other 3D-peaks, reducing the number of identified peptides. 
                Oversampling occurs in low-complexity samples or long LC gradients, as well as undersized dynamic 
                exclusion windows for data independent acquisitions.
                
                * Heatmap score [EVD: MS2 Oversampling]: The percentage of non-oversampled 3D-peaks.
                            ''',
            plot=bar_html
        )

    def draw_delta_mass(self):

        self.delta_mass_percent['target'] = dict(zip(self.delta_mass['target'].keys(),
                                                    list(map(lambda v: v / len(self.delta_mass['target']),
                                                            self.delta_mass['target'].values()))))
        if "decoy" in self.delta_mass.keys():
            self.delta_mass_percent['decoy'] = dict(zip(self.delta_mass['decoy'].keys(),
                                                    list(map(lambda v: v / len(self.delta_mass['decoy']),
                                                            self.delta_mass['decoy'].values()))))
            lineconfig = {
                # Building the plot
                'id': 'delta_mass',  # HTML ID used for plot
                "tt_label": "<b>{point.x} Mass delta relative frequency</b>: {point.y}",
                # Plot configuration
                'title': 'Delta m/z',  # Plot title - should be in format "Module Name: Plot Title"
                'xlab': 'Experimental m/z - Theoretical m/z',  # X axis label
                'ylab': 'Relative Frequency',  # Y axis label
                'colors': {'target': '#b2df8a', 'decoy': '#DC143C'},
                'xmax': max(list(self.delta_mass['decoy'].keys()) +
                            (list(self.delta_mass['target'].keys()))) + 0.01,
                'xmin': min(list(self.delta_mass['decoy'].keys()) +
                            (list(self.delta_mass['target'].keys()))) - 0.01,
                "data_labels": [{"name": "Counts", "ylab": "Count", "tt_label": "<b>{point.x} Mass delta counts</b>: {point.y}"},
                                {"name": "Relative Frequency", "ylab": "Relative Frequency", "tt_label": "<b>{point.x} Mass delta relative frequency</b>: {point.y}"}],
            }
        else:
            lineconfig = {
                # Building the plot
                'id': 'delta_mass',  # HTML ID used for plot
                "tt_label": "<b>{point.x} Mass delta relative frequency</b>: {point.y}",
                # Plot configuration
                'title': 'Delta m/z',  # Plot title - should be in format "Module Name: Plot Title"
                'xlab': 'Experimental m/z - Theoretical m/z',  # X axis label
                'ylab': 'Relative Frequency',  # Y axis label
                'colors': {'target': '#b2df8a'},
                'xmax': max(list(self.delta_mass['target'].keys())) + 0.01,
                'xmin': min((list(self.delta_mass['target'].keys()))) - 0.01,
                "data_labels": [{"name": "Counts", "ylab": "Count", "tt_label": "<b>{point.x} Mass delta counts</b>: {point.y}"},
                                {"name": "Relative Frequency", "ylab": "Relative Frequency", "tt_label": "<b>{point.x} Mass delta relative frequency</b>: {point.y}"}],
            }
        line_html = linegraph.plot([self.delta_mass, self.delta_mass_percent], lineconfig)

        self.add_section(
            name="Delta Mass",
            anchor="delta_mass",
            description='''This chart represents the distribution of the relative frequency of experimental 
                            precursor ion mass (m/z) - theoretical precursor ion mass (m/z). 
                            ''',
            helptext='''
                    Mass deltas close to zero reflect more accurate identifications and also 
                    that the reporting of the amino acid modifications and charges have been done accurately. 
                    This plot can highlight systematic bias if not centered on zero. 
                    Other distributions can reflect modifications not being reported properly. 
                    Also it is easy to see the different between the target and the decoys identifications.
                    ''',
            plot=line_html
        )

    def CalHeatMapScore(self):
        log.warning("Calculating Heatmap Scores...")
        mztab_data = mztab.MzTab(self.out_mzTab_path)
        psm = mztab_data.spectrum_match_table
        meta_data = dict(mztab_data.metadata)
        if self.pep_table_exists:
            pep_table = mztab_data.peptide_table
            pep_table.loc[:, 'stand_spectra_ref'] = pep_table.apply(
                lambda x: os.path.basename(meta_data[x.spectra_ref.split(':')[0] + '-location']), axis=1)
            study_variables = list(filter(lambda x: re.match(r'peptide_abundance_study_variable.*?', x) is not None,
                                        pep_table.columns.tolist()))

            for name, group in pep_table.groupby("stand_spectra_ref"):
                self.heatmap_con_score[name] = 1 - np.sum(np.sum(
                    group[group['accession'].str.contains(config.kwargs["contaminant_affix"])][study_variables])) / \
                                            np.sum(np.sum(group[study_variables]))
                if config.kwargs['remove_decoy']:
                    T = sum(group[(group['opt_global_cv_MS:1002217_decoy_peptide'] == 0)][study_variables].values. \
                            tolist(), [])
                else:
                    T = sum(group[study_variables].values.tolist(), [])
                pep_median = np.median([j for j in T if not math.isnan(j) is True])
                self.heatmap_pep_intensity[name] = np.minimum(1.0, pep_median / (2 ** 23))  # Threshold

        #  HeatMapMissedCleavages
        global_peps = set(psm['opt_global_cv_MS:1000889_peptidoform_sequence'])
        global_peps_count = len(global_peps)
        if config.kwargs['remove_decoy'] and 'opt_global_cv_MS:1002217_decoy_peptide' in psm.columns:
            psm = psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] == 0]
        psm.loc[:, 'stand_spectra_ref'] = psm.apply(
            lambda x: os.path.basename(meta_data[x.spectra_ref.split(':')[0] + '-location']), axis=1)
        psm.loc[:, 'missed_cleavages'] = psm.apply(lambda x: self.cal_MissCleavages(x['sequence']), axis=1)

        # Calculate the ID RT Score
        for name, group in psm.groupby('stand_spectra_ref'):
            sc = group['missed_cleavages'].value_counts()
            self.MissedCleavages_heatmap_score[name] = sc[0] / sc[:].sum()

            x = group['retention_time'] / np.sum(group['retention_time'])
            n = len(group['retention_time'])
            y = np.sum(x) / n
            worst = ((1 - y) ** 0.5) * 1 / n + (y ** 0.5) * (n - 1) / n
            sc = np.sum(np.abs(x - y) ** 0.5) / n
            if worst == 0:
                self.ID_RT_score[name] = 1
            else:
                self.ID_RT_score[name] = (worst - sc) / worst

            #  For HeatMapOverSamplingScore
            self.HeatmapOverSamplingScore[name] = self.oversampling[name]['1'] / np.sum(list(self.oversampling[name].values()))

            # For HeatMapPepMissingScore
            idFraction = len(
                set(group['opt_global_cv_MS:1000889_peptidoform_sequence']).intersection(
                    global_peps)) / global_peps_count
            self.HeatmapPepMissingScore[name] = np.minimum(1.0, idFraction)

        median = np.median(list(self.MissedCleavages_heatmap_score.values()))
        self.MissedCleavagesVar_score = dict(zip(self.MissedCleavages_heatmap_score.keys(),
                                                list(map(lambda v: 1 - np.abs(v - median),
                                                        self.MissedCleavages_heatmap_score.values()))))
        log.warning("Done calculating Heatmap Scores.")

    # if missed.cleavages is not given, it is assumed that trypsin was used for digestion
    @staticmethod
    def cal_MissCleavages(sequence):
        miss_cleavages = len(sequence[:-1]) - len(sequence[:-1].replace('K', '').replace('R', ''))
        return miss_cleavages

    @staticmethod
    def dis_decoy(ProteinName):
        if config.kwargs['decoy_affix'] not in ProteinName:
            return "TARGET"
        elif ProteinName.split(';') == 1:
            return "DECOY"
        else:
            if config.kwargs['affix_type'] == 'prefix':
                if list(filter(lambda x: lambda x: not x.startswith(config.kwargs['decoy_affix']),
                                ProteinName.split(';'))):
                    return "TARGET"
                return "DECOY"
            else:
                if list(filter(lambda x: not x.endswith(config.kwargs['decoy_affix']),
                                ProteinName.split(';'))):
                    return "TARGET"
                return "DECOY"

    def parse_mzml(self):

        self.mzml_peak_distribution_plot = Histogram('Peak Intensity', plot_category='range', breaks=[
                                                        0, 10, 100, 300, 500, 700, 900, 1000, 3000, 6000, 10000])

        self.mzml_charge_plot = Histogram('Precursor Charge', plot_category='frequency', breaks=[
                                            i for i in range(1, 8)])

        self.mzml_peaks_ms2_plot = Histogram('#Peaks per MS/MS spectrum', plot_category='range', breaks=[
                                            i for i in range(0, 1001, 100)])

        # New instances are used for dictionary construction.
        self.mzml_peak_distribution_plot_1 = copy.deepcopy(self.mzml_peak_distribution_plot)
        self.mzml_charge_plot_1 = copy.deepcopy(self.mzml_charge_plot)
        self.mzml_peaks_ms2_plot_1 = copy.deepcopy(self.mzml_peaks_ms2_plot)

        mzml_table = {}
        heatmap_charge = {}
        exp = MSExperiment()
        for m in self.mzML_paths:
            ms1_number = 0
            ms2_number = 0
            log.warning("{}: Parsing mzML file {}...".format(datetime.now().strftime("%H:%M:%S"),m))
            MzMLFile().load(m, exp)
            log.warning("{}: Done parsing mzML file {}...".format(datetime.now().strftime("%H:%M:%S"),m))
            m = os.path.basename(m)
            log.warning("{}: Aggregating mzML file {}...".format(datetime.now().strftime("%H:%M:%S"),m))
            charge_2 = 0
            for i in exp:
                if i.getMSLevel() == 1:
                    ms1_number += 1
                elif i.getMSLevel() == 2:
                    ms2_number += 1
                    charge_state = i.getPrecursors()[0].getCharge()
                    peak_per_ms2 = len(i.get_peaks()[0])
                    base_peak_intensity = i.getMetaValue("base peak intensity")
                    if charge_state == 2:
                        charge_2 += 1

                    if self.enable_dia:
                        self.mzml_charge_plot.addValue(charge_state)
                        self.mzml_peak_distribution_plot.addValue(base_peak_intensity)
                        self.mzml_peaks_ms2_plot.addValue(peak_per_ms2)
                        continue

                    if i.getNativeID() in self.identified_spectrum[m]:
                        self.mzml_charge_plot.addValue(charge_state)
                        self.mzml_peak_distribution_plot.addValue(base_peak_intensity)
                        self.mzml_peaks_ms2_plot.addValue(peak_per_ms2)
                    else:
                        self.mzml_charge_plot_1.addValue(charge_state)
                        self.mzml_peak_distribution_plot_1.addValue(base_peak_intensity)
                        self.mzml_peaks_ms2_plot_1.addValue(peak_per_ms2)

            heatmap_charge[m] = charge_2 / ms2_number
            self.total_ms2_spectra = self.total_ms2_spectra + ms2_number
            mzml_table[m] = {'MS1_Num': ms1_number}
            mzml_table[m]['MS2_Num'] = ms2_number
            log.warning("{}: Done aggregating mzML file {}...".format(datetime.now().strftime("%H:%M:%S"),m))
        
        self.mzml_charge_plot.to_dict()
        self.mzml_peaks_ms2_plot.to_dict()
        self.mzml_peak_distribution_plot.to_dict()
        # Construct compound dictionaries to apply to drawing functions.
        if self.enable_dia:
            self.mzml_info['charge_distribution'] = {
                'Whole Experiment': self.mzml_charge_plot.dict['data']
            }
            self.mzml_info['peaks_per_ms2'] = {
                'Whole Experiment': self.mzml_peaks_ms2_plot.dict['data']
            }
            self.mzml_info['peak_distribution'] = {
                'Whole Experiment': self.mzml_peak_distribution_plot.dict['data']
            }
        else:
            self.mzml_charge_plot_1.to_dict()
            self.mzml_peaks_ms2_plot_1.to_dict()
            self.mzml_peak_distribution_plot_1.to_dict()
            self.mzml_info['charge_distribution'] = {
                'identified_spectra': self.mzml_charge_plot.dict['data'],
                'unidentified_spectra': self.mzml_charge_plot_1.dict['data']
            }
            self.mzml_info['peaks_per_ms2'] = {
                'identified_spectra': self.mzml_peaks_ms2_plot.dict['data'],
                'unidentified_spectra': self.mzml_peaks_ms2_plot_1.dict['data']
            }
            self.mzml_info['peak_distribution'] = {
                'identified_spectra': self.mzml_peak_distribution_plot.dict['data'],
                'unidentified_spectra': self.mzml_peak_distribution_plot_1.dict['data']
            }

        median = np.median(list(heatmap_charge.values()))
        self.heatmap_charge_score = dict(zip(heatmap_charge.keys(),
                                            list(map(lambda v: 1 - np.abs(v - median),
                                                    heatmap_charge.values()))))
        
        return mzml_table
                                
    def parse_idxml(self, mzml_table):
        for raw_id in self.idx_paths:
            log.warning("Parsing {}...".format(raw_id))
            protein_ids = []
            peptide_ids = []
            IdXMLFile().load(raw_id, protein_ids, peptide_ids)
            raw_id = os.path.basename(raw_id)
            ## TODO I would use the QC functionality of pyopenms. Should be much faster.
            if config.kwargs['remove_decoy']:
                identified_num = len(set([i.getMetaValue("spectrum_reference") for i in peptide_ids
                                                    if i.getHits()[0].getMetaValue("target_decoy") == 'target']))
            else:
                identified_num = len(peptide_ids)

            ## TODO make clear if this is before or after first step of filtering
            mzML_name = os.path.basename(protein_ids[0].getMetaValue("spectra_data")[0].decode("UTF-8"))
            search_engine = protein_ids[0].getSearchEngine()
            if search_engine == "MSGF+" or "msgf" in raw_id:
                mzml_table[mzML_name]['MSGF'] = identified_num
            elif search_engine == "Comet" or "comet" in raw_id:
                mzml_table[mzML_name]['Comet'] = identified_num
            else:
                mzml_table[mzML_name][search_engine] = identified_num

            mzml_table[mzML_name]['num_quant_psms'] = self.mL_spec_ident_final[mzML_name]
            mzml_table[mzML_name]['num_quant_peps'] = len(self.mzml_peptide_map[mzML_name])

        # mzml file sorted based on experimental file
        for mzML_name in self.exp_design_table.keys():
            self.mzml_table[mzML_name] = mzml_table[mzML_name]

    def parse_out_mzTab(self):
        ## TODO why is this a warning??
        log.warning("{}: Parsing mzTab file {}...".format(datetime.now().strftime("%H:%M:%S"),self.out_mzTab_path))
        mztab_data = mztab.MzTab(self.out_mzTab_path)
        log.warning("{}: Done parsing mzTab file {}.".format(datetime.now().strftime("%H:%M:%S"),self.out_mzTab_path))
        log.warning("{}: Aggregating mzTab file {}...".format(datetime.now().strftime("%H:%M:%S"),self.out_mzTab_path))
        pep_table = mztab_data.peptide_table
        meta_data = dict(mztab_data.metadata)

        self.delta_mass['target'] = dict()
        self.delta_mass['decoy'] = dict()

        # PSM table data
        psm = mztab_data.spectrum_match_table
        prot = mztab_data.protein_table
        self.prot_search_score = dict()
        # TODO what is "stand_spectra_ref". What does it mean? What is it used for?
        psm['stand_spectra_ref'] = psm.apply(
            lambda x: os.path.basename(meta_data[x.spectra_ref.split(':')[0] + '-location']), axis=1)

        prot_abundance_cols = list(filter(lambda x: re.match(r'protein_abundance_assay.*?', x) is not None,
                                    prot.columns.tolist()))
        opt_cols = list(filter(lambda x: x.startswith("opt_"),
                                    prot.columns.tolist()))
        score_cols = list(filter(lambda x: x.startswith("best_search_engine_score"),
                                    prot.columns.tolist()))
        # TODO in theory we do not need accession since the index is the accession
        fixed_cols = ['accession', 'description', 'taxid', 'species', 'database', 'database_version', 'search_engine', 'ambiguity_members', 'modifications', 'protein_coverage']
                                                                    
        prot = prot[fixed_cols + score_cols + prot_abundance_cols + opt_cols]

        # We only need the overall protein (group) scores and abundances. Currently we do not care about details of single proteins (length, description,...)
        prot = prot[prot['opt_global_result_type'] != 'protein_details']

        if config.kwargs['remove_decoy']:
            psm = psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] != 1]
            # TODO do we really want to remove groups that contain a single decoy? I would say ALL members need to be decoy.
            prot = prot[~prot['accession'].str.contains(config.kwargs['decoy_affix'])]

        prot["protein_group"] = prot["ambiguity_members"].apply(lambda x: x.replace(",", ";"))
        
        # TODO the following assumes that we always only look
        peptide_score = pep_table[["opt_global_cv_MS:1000889_peptidoform_sequence", "best_search_engine_score[1]"]]
        self.peptide_search_score = peptide_score.groupby("opt_global_cv_MS:1000889_peptidoform_sequence").agg('min')["best_search_engine_score[1]"].to_dict()
        del peptide_score
        
        self.Total_Protein_Identified = len(prot.index)

        prot.dropna(how='all',subset=prot_abundance_cols, inplace=True)
        self.Total_Protein_Quantified = len(prot.index)

        self.pep_plot = Histogram('Number of peptides per protein', plot_category = 'frequency')

        # There probably are no shared peptides in the final quant results. We still do it to be safe.
        # I think the accessions here are actually also the leader accessions in case of groups
        counts_per_acc = pep_table['accession'].str.split(",").explode().value_counts()
        counts_per_acc.apply(self.pep_plot.addValue)
        #for c in counts_per_acc:
        #    self.pep_plot.addValue(c)
        
        categories = OrderedDict()
        categories['Frequency'] = {
            'name': 'Frequency',
            'description': 'Number of identified peptides per protein.'
        }
        self.pep_plot.to_dict(percentage = True, cats = categories)

        mL_spec_ident_final = {}

        for m, group in psm.groupby('stand_spectra_ref'):
            self.cal_num_table_data[m] = {'Sample Name': self.exp_design_table[m]['Sample']}
            self.cal_num_table_data[m]['Condition'] = self.exp_design_table[m]['MSstats_Condition']
            self.cal_num_table_data[m]['Fraction'] = self.exp_design_table[m]['Fraction']

            if config.kwargs['remove_decoy']:
                group = group[group['opt_global_cv_MS:1002217_decoy_peptide'] == 0]
            # Each loop resets the instance.
            self.oversampling_plot = Histogram('MS/MS counts per 3D-peak', plot_category = 'frequency', breaks = [1, 2, 3])

            for i in group.value_counts(['sequence', 'charge']).values:
                self.oversampling_plot.addValue(i)

            self.oversampling_plot.to_dict()
            self.oversampling[m] = self.oversampling_plot.dict['data']

            proteins = set(group['accession'])
            peptides = set(group['opt_global_cv_MS:1000889_peptidoform_sequence'])
            unique_peptides = set(group[group['unique'] == 1]['opt_global_cv_MS:1000889_peptidoform_sequence'])

            self.identified_spectrum[m] = list(map(lambda x: x.split(':')[1],
                                                group['spectra_ref']))
            self.mzml_peptide_map[m] = list(set(group['sequence'].tolist()))

            if None in proteins:
                proteins.remove(None)
            
            ## TODO this is not really the number of proteins but the number of protein groups
            self.cal_num_table_data[m]['protein_num'] = len(proteins)
            self.cal_num_table_data[m]['peptide_num'] = len(peptides)
            self.cal_num_table_data[m]['unique_peptide_num'] = len(unique_peptides)

            modified_pep = list(filter(lambda x: re.match(r'.*?\(.*\).*?', x) is not None, peptides))
            self.cal_num_table_data[m]['modified_peptide_num'] = len(modified_pep)

            mL_spec_ident_final[m] = len(set(self.identified_spectrum[m]))
        

        target_bin_data = {}
        decoy_bin_data = {}
        # TODO This is NOT a relative difference!
        psm['relative_diff'] = psm['exp_mass_to_charge'] - psm['calc_mass_to_charge']
        try:
            decoy_bin = psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] == 1]['relative_diff'].value_counts(sort=False, bins=1000)
            for index in decoy_bin.index:
                decoy_bin_data[float(index.mid)] = int(decoy_bin[index])
            self.delta_mass['decoy'] = decoy_bin_data
        except Exception as e:
            print("No decoy peptides found -> only showing target peptides")

        target_bin = psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] != 1]['relative_diff'].value_counts(sort=False, bins=1000)    
        for index in target_bin.index:
            target_bin_data[float(index.mid)] = int(target_bin[index])

        self.delta_mass['target'] = target_bin_data

        # extract delta mass
        self.mL_spec_ident_final = mL_spec_ident_final
        if config.kwargs['remove_decoy']:
            self.Total_ms2_Spectral_Identified = len(set(psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] != 1]
                                                        ['spectra_ref']))
            self.Total_Peptide_Count = len(set(psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] != 1]
                                                ['sequence']))
        else:
            self.Total_ms2_Spectral_Identified = len(set(psm['spectra_ref']))
            self.Total_Peptide_Count = len(set(psm['sequence']))

        log.warning("{}: Done aggregating mzTab file {}...".format(datetime.now().strftime("%H:%M:%S"),self.out_mzTab_path))

    def parse_diann_report(self):
        log.warning("Parsing {}...".format(self.diann_report_path))
        pattern = re.compile(r"\(.*?\)")
        report_data = pd.read_csv(self.diann_report_path, header=0, sep="\t")
        report_data["sequence"] = report_data.apply(lambda x: re.sub(pattern,"",x["Modified.Sequence"]), axis=1)
        self.Total_Protein_Quantified = len(set(report_data["Protein.Names"]))
        self.Total_Peptide_Count = len(set(report_data["sequence"]))
        protein_pep_map = report_data.groupby('Protein.Group').sequence.apply(list).to_dict()

        self.pep_plot = Histogram('number of peptides per proteins', plot_category = 'frequency')

        for _, peps in protein_pep_map.items():
            number = len(set(peps))
            self.pep_plot.addValue(number)

        self.peptide_search_score = dict()
        pattern = re.compile(r"\((.*?)\)")
        unimod_data = UnimodDatabase()
        for peptide, group in report_data.groupby("Modified.Sequence"):
            origianl_mods = re.findall(pattern, peptide)
            for mod in set(origianl_mods):
                name = unimod_data.get_by_accession(mod.upper()).get_name()
                peptide = peptide.replace(mod, name)
            if peptide.startswith("("):
                peptide = peptide + "."

            self.peptide_search_score[peptide] = np.min(group["Q.Value"])

        categorys = OrderedDict()
        categorys['Frequency'] = {
            'name': 'Frequency',
            'description': 'number of peptides per proteins'
        }
        self.pep_plot.to_dict(percentage = True, cats = categorys)

        for run_file, group in report_data.groupby("File.Name"):
            run_file = os.path.basename(run_file)
            self.cal_num_table_data[run_file] = {'Sample Name': self.exp_design_table[run_file]['Sample']}
            self.cal_num_table_data[run_file]['condition'] = self.exp_design_table[run_file]['MSstats_Condition']
            self.cal_num_table_data[run_file]['fraction'] = self.exp_design_table[run_file]['Fraction']
            self.cal_num_table_data[run_file]['protein_num'] = len(set(group["Protein.Ids"]))
            self.cal_num_table_data[run_file]['peptide_num'] = len(set(group["sequence"]))
            peptides = set(group["Modified.Sequence"])
            modified_pep = list(filter(lambda x: re.match(r".*?\(.*?\).*?", x) is not None, peptides))
            group_peptides = group.groupby('sequence')["Protein.Group"].apply(list).to_dict()
            unique_peptides = [pep for pep, prots in group_peptides.items() if len(set(prots)) == 1]
            self.cal_num_table_data[run_file]['unique_peptide_num'] = len(unique_peptides)
            self.cal_num_table_data[run_file]['modified_peptide_num'] = len(modified_pep)

    def parse_msstats_input(self):
        print("Parsing MSstats input file " + self.msstats_input_path)
        msstats_data = pd.read_csv(self.msstats_input_path)
        ## TODO we probably shouldn't even write out 0-intensity values to MSstats csv
        msstats_data = msstats_data[-(msstats_data["Intensity"] == 0)]

        # multiQC requires weird dicts
        msstats_data_dict_pep_full = dict()
        msstats_data_dict_pep_init  = dict()
        msstats_data_dict_prot_full = dict()
        msstats_data_dict_prot_init = dict()

        max_pep_intensity = 0.0

        repsPerCondition = self.sample_df.groupby("MSstats_Condition")["MSstats_BioReplicate"].agg(list).to_dict()
        conditions = list(self.sample_df["MSstats_Condition"].unique())
        conditions_str = [str(c) for c in conditions]
        conditions_dists = [str(c) + "_distribution" for c in conditions]
        cond_and_dist_cols = conditions_str + conditions_dists

        # TODO maybe aggregating in dicts is not the fastest. We also need to parse them again for proteins later.
        #  Maybe we should just introduce new pandas columns for every bioreplicate.
        def fillDict(g):
            d = dict.fromkeys(repsPerCondition[str(g.name)], 0)
            d.update(zip(g["BioReplicate"].astype(str), np.log10(g["Intensity"])))
            return json.dumps(d)

        def getIntyAcrossBioRepsAsStr(g):
            gdict = dict.fromkeys(conditions_str, 0.0)
            gdict.update(dict.fromkeys(conditions_dists, '{}'))
            gdict["ProteinName"] = g["ProteinName"].iloc[0]
            gdict["Average Intensity"] = np.log10(g["Intensity"].mean())
            ## TODO How to determine technical replicates? Should be same BioReplicate but different Fraction_Group (but fraction group is not annotated)
            condGrp = g.groupby(["Condition","BioReplicate"])["Intensity"].mean().reset_index().groupby("Condition").apply(fillDict)
            condGrp.index = [str(c) + "_distribution" for c in condGrp.index]
            gdict.update(condGrp.to_dict())
            mean = g.groupby(["Condition"])["Intensity"].mean()
            #print(mean)
            condGrpMean = np.log10(mean)
            condGrpMean.index = condGrpMean.index.map(str)
            gdict.update(condGrpMean.to_dict())
            return pd.Series(gdict)


        msstats_data_pep_agg = msstats_data.groupby(["PeptideSequence"]).apply(getIntyAcrossBioRepsAsStr)#.unstack()
        ## TODO Can we guarantee that the score was always PEP? I don't think so!
        msstats_data_pep_agg["BestSearchScore"] = 1 - msstats_data_pep_agg.index.map(self.peptide_search_score)
        msstats_data_dict_pep_full = msstats_data_pep_agg.to_dict('index')
        msstats_data_dict_pep_init = dict(itertools.islice(msstats_data_dict_pep_full.items(), 50))
        msstats_data_pep_agg.reset_index(inplace=True)


        cur.execute("CREATE TABLE PEPQUANT(PeptideSequence VARCHAR(100) PRIMARY KEY, ProteinName VARCHAR(100), BestSearchScore FLOAT(4,3), \"Average Intensity\" FLOAT(4,3))")
        con.commit()
        sql_col = "PeptideSequence,ProteinName,BestSearchScore, \"Average Intensity\""
        sql_t = "(" + ','.join(['?'] * (len(conditions) *2 + 4)) + ")"

        headers = OrderedDict()
        headers = { #'PeptideSequence': {'name': 'PeptideSequence'}, # this is the index
                    'ProteinName': {'name': 'ProteinName'},
                    'BestSearchScore': {'name': 'BestSearchScore', 'format': '{:,.5f}'},
                    'Average Intensity': {'name': 'Average Intensity', 'format': '{:,.3f}'}}

        
        for s in conditions:
            cur.execute("ALTER TABLE PEPQUANT ADD \"" + str(s) + "\" VARCHAR")
            con.commit()
            sql_col += ", \"" + str(s) + "\""
            headers[str(s)] = {'name': s, 'format': '{:,.5f}'}

        for s in list(map(lambda x: str(x) + "_distribution", conditions)):
            cur.execute("ALTER TABLE PEPQUANT ADD \"" + s + "\" VARCHAR(100)")
            con.commit()
            sql_col += ", \"" + s + "\""
            # we need a thousandsSep_format otherwise commas will be replaced
            headers[str(s)] = {'name': s, 'thousandsSep_format': ','}

        # PeptideSequence is index
        all_term = ["ProteinName", "BestSearchScore", "Average Intensity"] + list(map(str, conditions)) + list(map(lambda x: str(x) + "_distribution", conditions))
        cur.executemany("INSERT INTO PEPQUANT (" + sql_col + ") VALUES " + sql_t,
                        [(k, *itemgetter(*all_term)(v)) for k,v in msstats_data_dict_pep_full.items()])
        con.commit()


        pconfig = {
            'id': 'quantification_of_peptides',  # ID used for the table
            'table_title': 'quantification information of peptides',
            # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'raw_data_fn': 'multiqc_quantification_of_peptides_table',  # File basename to use for raw data file
            'sortRows': False,  # Whether to sort rows alphabetically
            'only_defined_headers': False,  # Only show columns that are defined in the headers config
            'col1_header': 'PeptideSequence',
            'thousandsSep_format': ",",
            'no_beeswarm': True,
            'shared_key': None
        }

        # only use the first 50 lines for the table
        table_html = sparklines.plot(msstats_data_dict_pep_init, headers, pconfig=pconfig, maxValue=max_pep_intensity)
        pattern = re.compile(r'<small id="quantification_of_peptides_numrows_text"')
        index = re.search(pattern, table_html).span()[0]
        t_html = table_html[:index] + '<input type="text" placeholder="search..." class="searchInput" ' \
                                    'onkeyup="searchQuantFunction()" id="quant_search">' \
                                    '<select name="quant_search_col" id="quant_search_col">'
        for key in ["ProteinName", "PeptideSequence"]:
            t_html += '<option>' + key + '</option>'
        table_html = t_html + '</select>' + '<button type="button" class="btn btn-default ' \
                                            'btn-sm" id="quant_reset" onclick="quantFirst()">Reset</button>' \
                    + table_html[index:]
        table_html = table_html + '''<div class="page_control"><span id="quantFirst">First Page</span><span 
        id="quantPre"> Previous Page</span><span id="quantNext">Next Page </span><span id="quantLast">Last 
        Page</span><span id="quantPageNum"></span>Page/Total <span id="quantTotalPage"></span>Pages <input 
        type="number" name="" id="pep_page" class="page" value="" oninput="this.value=this.value.replace(/\D/g);" 
        onkeydown="quant_page_jump()" min="1"/> </div> '''

        # Add a report section with the line plot
        self.add_section(
            name="Peptides Quantification Table",
            anchor="peptides_quant_result",
            description='This plot shows the quantification information of peptides'
                        'in quantms pipeline final result',
            helptext='''
                    The quantification information of peptides is obtained from the MSstats input file. 
                    The table shows the quantitative level and distribution of peptides in different study variables, run and peptiforms. The distribution show all the intensity values in a bar plot above and below the average intensity for all the fractions, runs and peptiforms.

                    * BestSearchScore: It is equal to 1 - min(Q.Value) for DIA datasets. Then it is equal to 1 - min(best_search_engine_score[1]), which is from best_search_engine_score[1] column in mzTab peptide table for DDA datasets.
                    * Average Intensity: Average intensity of each peptide sequence across all conditions with NA=0 or NA ignored.
                    * Peptide intensity in each condition (Eg. `CT=Mixture;CN=UPS1;QY=0.1fmol`): Summarize intensity of fractions, and then mean intensity in technical replicates/biological replicates separately. Click `distribution` to switch to bar plots.

                    ''',
            plot=table_html
        )

        # Helper functions for pandas
        def jsonToDict(s):
            if type(s) is str:
                return json.loads(s)
            else:
                return {}

        def reducer(accumulator, element):
            for key, value in jsonToDict(element).items():
                accumulator[key] = accumulator.get(key, 0) + value
            return accumulator

        def myDictSum(series):
            return json.dumps(reduce(reducer, series, {}))

        max_prot_intensity = 0
        agg_funs = dict.fromkeys(conditions_dists, myDictSum)
        agg_funs.update(dict.fromkeys(conditions_str, 'sum'))
        agg_funs["PeptideSequence"] = 'count'
        agg_funs["Average Intensity"] = 'sum'
        msstats_data_prot = msstats_data_pep_agg.groupby("ProteinName").agg(agg_funs)#.reset_index()
        del(msstats_data_pep_agg)
        msstats_data_prot.rename(columns={"PeptideSequence": "Peptides_Number"},inplace=True)
        msstats_data_dict_prot_full = msstats_data_prot.to_dict('index')

        msstats_data_dict_prot_init = dict(itertools.islice(msstats_data_dict_prot_full.items(), 50))

        headers = OrderedDict()
        # This will be the rowheader/index
        #headers['ProteinName'] = {
        #    'name': 'Protein Name',
        #    'description': 'Name/Identifier(s) of the protein (group)',
        #    'format': '{:,.0f}'
        #}
        headers['Peptides_Number'] = {
            'name': 'Number of Peptides',
            'description': 'Number of peptides per proteins',
            'format': '{:,.0f}'
        }
        headers['Average Intensity'] = {
            'name': 'Average Intensity',
            'description': 'Average intensity across all conditions',
            'format': '{:,.0f}'
        }

        # upload protein table to sqlite database
        cur.execute("CREATE TABLE PROTQUANT(ProteinName VARCHAR(100), Peptides_Number INT(100), \"Average Intensity\" VARCHAR)")
        con.commit()
        sql_col = "ProteinName,Peptides_Number,\"Average Intensity\""
        sql_t = "(" + ','.join(['?'] * (len(conditions)*2 + 3)) + ")"

        for s in conditions:
            cur.execute("ALTER TABLE PROTQUANT ADD \"" + str(s) + "\" VARCHAR")
            con.commit()
            sql_col += ", \"" + str(s) + "\""
            headers[str(s)] = {'name': s}

        for s in list(map(lambda x: str(x) + "_distribution", conditions)):
            cur.execute("ALTER TABLE PROTQUANT ADD \"" + s + "\" VARCHAR(100)")
            con.commit()
            sql_col += ", \"" + s + "\""
            headers[str(s)] = {'name': s}

        # ProteinName is index
        all_term = ["Peptides_Number", "Average Intensity"] + list(map(str, conditions)) + list(map(lambda x: str(x) + "_distribution", conditions))
        cur.executemany("INSERT INTO PROTQUANT (" + sql_col + ") VALUES " + sql_t,
        [(k, *itemgetter(*all_term)(v)) for k,v in msstats_data_dict_prot_full.items()])
        con.commit()

        pconfig = {
            'id': 'quantification_of_protein',  # ID used for the table
            'table_title': 'quantification information of protein',
            # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'raw_data_fn': 'multiqc_quantification_of_protein_table',  # File basename to use for raw data file
            'sortRows': False,  # Whether to sort rows alphabetically
            'only_defined_headers': False,  # Only show columns that are defined in the headers config
            'col1_header': 'ProteinName',
            'format': '{:,.5f}',
            'no_beeswarm': True
        }

        table_html = sparklines.plot(msstats_data_dict_prot_init, headers, pconfig=pconfig, maxValue=max_prot_intensity)
        pattern = re.compile(r'<small id="quantification_of_protein_numrows_text"')
        index = re.search(pattern, table_html).span()[0]
        t_html = table_html[:index] + '<input type="text" placeholder="search..." class="searchInput" ' \
                                    'onkeyup="searchProtFunction()" id="prot_search">' \
                                    '<select name="prot_search_col" id="prot_search_col">'
        for key in ["ProteinName"]:
            t_html += '<option>' + key + '</option>'
        table_html = t_html + '</select>' + '<button type="button" class="btn btn-default ' \
                                            'btn-sm" id="prot_reset" onclick="protFirst()">Reset</button>' \
                    + table_html[index:]
        table_html = table_html + '''<div class="page_control"><span id="protFirst">First Page</span><span 
        id="protPre"> Previous Page</span><span id="protNext">Next Page </span><span id="protLast">Last 
        Page</span><span id="protPageNum"></span>Page/Total <span id="protTotalPage"></span>Pages <input 
        type="number" name="" id="prot_page" class="page" value="" oninput="this.value=this.value.replace(/\D/g);" 
        onkeydown="prot_page_jump()" min="1"/> </div> '''

        self.add_section(
            name="Protein Quantification Table",
            anchor="protein_quant_result",
            description='This plot shows the quantification information of proteins'
                        ' in final result (mainly the mzTab file).',
            helptext='''
                    The quantification information of proteins is obtained from the msstats input file. 
                    The table shows the quantitative level and distribution of proteins in different study variables and run.

                    * Peptides_Number: The number of peptides for each protein.
                    * Average Intensity: Average intensity of each protein across all conditions with NA=0 or NA ignored.
                    * Protein intensity in each condition (Eg. `CT=Mixture;CN=UPS1;QY=0.1fmol`): Summarize intensity of peptides.
                    
                    Click `Show replicates` to switch to bar plots of quantities in each replicate.
                    ''',
            plot=table_html
        )

def read_openms_design(desfile): 
        with open(desfile, 'r') as f:
            data = f.readlines()
            s_row = False
            f_table = []
            s_table = []
            for row in data:
                if row == "\n":
                    continue                
                if "MSstats_Condition" in row:
                    s_row = True
                    s_header = row.replace('\n', '').split('\t')
                elif s_row:
                    s_table.append(row.replace('\n', '').split('\t'))
                elif "Spectra_Filepath" in row:
                    f_header = row.replace('\n', '').split('\t')
                else:
                    f_table.append(row.replace('\n', '').split('\t'))

            f_table = pd.DataFrame(f_table, columns=f_header)
            f_table["Run"] = f_table.apply(lambda x: os.path.splitext(os.path.basename(x["Spectra_Filepath"]))[0], axis=1)
            s_DataFrame = pd.DataFrame(s_table, columns=s_header)

        return s_DataFrame, f_table