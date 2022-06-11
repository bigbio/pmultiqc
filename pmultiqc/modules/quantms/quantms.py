#!/usr/bin/env python

""" MultiQC pmultiqc plugin module """

from __future__ import absolute_import
from collections import OrderedDict
from operator import itemgetter
import logging
from sdrf_pipelines.openms.openms import OpenMS
from multiqc import config
from multiqc.plots import table, bargraph, linegraph, heatmap
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.utils.mqc_colour import mqc_colour_scale
import pandas as pd
import re
from pyteomics import mztab
from pyopenms import IdXMLFile, MzMLFile, MSExperiment
import os
import sqlite3
import numpy as np
import math
import copy
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


class QuantMSModule(BaseMultiqcModule):
    def __init__(self):

        # Initialise the parent module Class object
        super(QuantMSModule, self).__init__(
            name='pmultiqc',
            target="pmultiqc",
            anchor='pmultiqc',
            href='https://github.com/bigbio/pmultiqc',
            info=" is an module to show the pipeline performance."
        )

        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_plugin', True):
            return None

        self.enable_exp = False
        self.enable_sdrf = False
        for f in self.find_log_files("quantms/exp_design"):
            self.exp_design = os.path.join(f["root"], f['fn'])
            self.enable_exp = True
        if self.enable_exp == False:
            for f in self.find_log_files("quantms/sdrf"):
                self.sdrf = os.path.join(f["root"], f['fn'])
                OpenMS().openms_convert(self.sdrf, config.kwargs['raw'],
                                    False, True, False, config.kwargs['condition'])
                self.enable_sdrf = True

        if self.enable_sdrf == False and self.enable_exp == False:
            raise AttributeError("exp_design and sdrf cannot be empty at the same time")

        self.PSM_table = dict()
        self.mzml_peptide_map = dict()
        self.identified_spectrum = dict()
        self.pep_quant_table = dict()
        self.mzml_table = OrderedDict()
        self.Total_ms2_Spectral = 0
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
        for f in self.find_log_files("quantms/mztab"):
            self.out_mzTab_path = os.path.join(f["root"], f['fn'])
            self.parse_out_mzTab()

        for report in self.find_log_files("quantms/diann_report"):
            self.diann_report_path = os.path.join(report["root"], report["fn"])
            self.enable_dia = True

        for msstats_input in self.find_log_files("quantms/msstats"):
            self.msstats_input_path = os.path.join(msstats_input["root"], msstats_input["fn"])

        self.mzML_paths = []
        for mzML_file in self.find_log_files("quantms/mzML"):
            self.mzML_paths.append(os.path.join(mzML_file["root"], mzML_file['fn']))

        mt = self.parse_mzml()
        self.idx_paths = []
        for idx_file in self.find_log_files("quantms/idXML"):
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
            description='This plot shows the pipeline performance overview',
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
        with open(self.exp_design, 'r') as f:
            data = f.readlines()
            s_row = False
            self.f_table = []
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
                    self.f_table.append(row.replace('\n', '').split('\t'))

            self.f_table = pd.DataFrame(self.f_table, columns=f_header)
            self.f_table["Run"] = self.f_table.apply(lambda x: os.path.splitext(os.path.basename(x["Spectra_Filepath"]))[0], axis=1)
            s_DataFrame = pd.DataFrame(s_table, columns=s_header)
            for file in np.unique(self.f_table["Spectra_Filepath"].tolist()):
                stand_file = os.path.basename(file)
                file_index = self.f_table[self.f_table["Spectra_Filepath"] == file]
                self.exp_design_table[stand_file] = {'Fraction_Group': file_index["Fraction_Group"].tolist()[0]}
                self.exp_design_table[stand_file]['Fraction'] = file_index["Fraction"].tolist()[0]
                self.exp_design_table[stand_file]['Label'] = '|'.join(file_index["Label"])
                sample = file_index["Sample"].tolist()
                self.exp_design_table[stand_file]['Sample'] = '|'.join(sample)
                self.exp_design_table[stand_file]['MSstats_Condition'] = ','.join(
                    [row['MSstats_Condition'] for _, row in s_DataFrame.iterrows()
                    if row['Sample'] in sample])
                self.exp_design_table[stand_file]['MSstats_BioReplicate'] = '|'.join(
                    [row['MSstats_BioReplicate'] for _, row in s_DataFrame.iterrows()
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
        maxnr = len(self.f_table.index)
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
            description='This plot shows the Proteomics Experimental Design',
            helptext='''
            This plot shows the Proteomics Experimental Design. You can see details about it in 
            https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/classOpenMS_1_1ExperimentalDesign.html
            ''',
            plot=table_html
        )

    def draw_summary_protein_ident_table(self):
        headers = OrderedDict()
        if self.enable_dia:
            summary_table = {
                self.Total_ms2_Spectral: {"Total Peptide Quantified": self.Total_Peptide_Count}}
            summary_table[self.Total_ms2_Spectral]['Total Protein Quantified'] = self.Total_Protein_Quantified    
        else:
            summary_table = {
                self.Total_ms2_Spectral: {"Total MS/MS Spectral Identified": self.Total_ms2_Spectral_Identified}}
            coverage = self.Total_ms2_Spectral_Identified / self.Total_ms2_Spectral * 100
            summary_table[self.Total_ms2_Spectral]['Identified MS/MS Spectral Coverage'] = coverage
            summary_table[self.Total_ms2_Spectral]['Total Peptide Identified'] = self.Total_Peptide_Count
            summary_table[self.Total_ms2_Spectral]['Total Protein Identified'] = self.Total_Protein_Identified
            summary_table[self.Total_ms2_Spectral]['Total Protein Quantified'] = self.Total_Protein_Quantified
            headers['Total MS/MS Spectral Identified'] = {
                'description': 'Total MS/MS Spectral Identified',
            }
            headers['Identified MS/MS Spectral Coverage'] = {
                'description': 'Identified MS/MS Spectral Coverage',
                'format': '{:,.2f}',
                "suffix": "%"
            }
        # Create table plot
        pconfig = {
            'id': 'identificaion summary table',  # ID used for the table
            'table_title': 'summary table',  # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'raw_data_fn': 'multiqc_summary_table_table',  # File basename to use for raw data file
            'sortRows': False,  # Whether to sort rows alphabetically
            'only_defined_headers': False,  # Only show columns that are defined in the headers config
            'col1_header': 'Total MS/MS Spectral',
            'format': '{:,.0f}',  # The header used for the first column
            "scale": "Set1"
        }
        table_html = table.plot(summary_table, headers, pconfig)

        # Add a report section with the line plot
        self.add_section(
            name="Summary Table",
            anchor="quantms_summary_table",
            description='This plot shows the quantms pipeline summary statistics',
            helptext='''
                    This plot shows the quantms pipeline summary statistics
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
        headers['Sample Name'] = {
            'description': 'Sample identifier',
            'color': "#ffffff"
        }
        headers['condition'] = {
                'description': 'Possible Study Variables',
                'color': "#ffffff"
        }
        headers['fraction'] = {
            'description': 'fraction identifier',
            'color': "#ffffff"
        }
        headers['peptide_num'] = {
            'description': 'identified the number of peptide in the pipeline',
            'color': "#ffffff"
        }
        headers['unique_peptide_num'] = {
            'description': 'identified the number of unique peptide in the pipeline',
            'color': "#ffffff"
        }
        headers['modified_peptide_num'] = {
            'description': 'identified the number of modified peptide in the pipeline',
            'color': "#ffffff"
        }
        headers['protein_num'] = {
            'description': 'identified the number of protein in the pipeline',
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
            'table_title': 'pipeline mzML tracking',  # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'raw_data_fn': 'multiqc_mzML_tracking_table',  # File basename to use for raw data file
            'sortRows': False,  # Whether to sort rows alphabetically
            'only_defined_headers': False,  # Only show columns that are defined in the headers config
            'col1_header': 'Spectra File',
            'format': '{:,.0f}'  # The header used for the first column
        }
        headers = OrderedDict()
        headers['MS1_Num'] = {
            'description': 'MS1 number',
            'color': "#ffffff"
        }
        headers['MS2_Num'] = {
            'description': 'MS2 number',
            'color': "#ffffff"
        }
        if 'MSGF' in self.mzml_table.values():
            headers['MSGF'] = {
                'description': 'Number of spectra identified by MSGF search engine',
                'color': "#ffffff"
            }
        if 'Comet' in self.mzml_table.values():
            headers['Comet'] = {
                'description': 'Number of spectra identified by Comet search engine',
                'color': "#ffffff"
            }
        headers['Final result of spectra'] = {
            'description': 'final number of spectra identified ',
            'color': "#ffffff"
        }
        table_html = table.plot(self.mzml_table, headers, pconfig)

        # Add a report section with the line plot
        self.add_section(
            name="Spectra Tracking",
            anchor="spectra_tracking",
            description='This plot shows the quantms pipeline mzML tracking',
            helptext='''
                    This table shows the changes in the number of spectra corresponding to each input file 
                    during the pipeline operation. And the number of peptides finally identified is obtained from 
                    the PSM table in the mzTab file. You can also remove the decoy with the `remove_decoy` parameter.:

                    * MS1_Num: The number of MS1 spectra extracted from mzMLs
                    * MS2_Num: The number of MS2 spectra extracted from mzMLs
                    * MSGF: The Number of spectra identified by MSGF search engine
                    * Comet: The Number of spectra identified by Comet search engine
                    * Final result of spectra: extracted from PSM table in mzTab file
                    * Final result of Peptides: extracted from PSM table in mzTab file
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
                "data_labels": [{"name": "Counts", "ylab": "Count"},
                                {"name": "Relative Frequency", "ylab": "Relative Frequency"}],
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
                "data_labels": [{"name": "Counts", "ylab": "Count"},
                                {"name": "Relative Frequency", "ylab": "Relative Frequency"}],
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
        log.warning("Calculate Heatmap Score")
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

        self.mzml_peak_distribution_plot = Histogram('Peak intensity', plot_category='range', breaks=[
                                                        0, 10, 100, 300, 500, 700, 900, 1000, 3000, 6000, 10000])

        self.mzml_charge_plot = Histogram('Precursor charge state', plot_category='frequency', breaks=[
                                            i for i in range(1, 8)])

        self.mzml_peaks_ms2_plot = Histogram('Number of Peaks per MS/MS spectrum', plot_category='range', breaks=[
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
            log.warning("Parsing {}...".format(m))
            MzMLFile().load(m, exp)            
            m = os.path.basename(m)
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
            self.Total_ms2_Spectral = self.Total_ms2_Spectral + ms2_number
            mzml_table[m] = {'MS1_Num': ms1_number}
            mzml_table[m]['MS2_Num'] = ms2_number
        
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
            if config.kwargs['remove_decoy']:
                identified_num = len(set([i.getMetaValue("spectrum_reference") for i in peptide_ids
                                                    if i.getHits()[0].getMetaValue("target_decoy") == 'target']))
            else:
                identified_num = len(peptide_ids)

            mzML_name = os.path.basename(protein_ids[0].getMetaValue("spectra_data")[0].decode("UTF-8"))
            search_engine = protein_ids[0].getSearchEngine()
            if search_engine== "MSGF+" or "msgf" in raw_id:
                mzml_table[mzML_name]['MSGF'] = identified_num
            elif search_engine == "Comet" or "comet" in raw_id:
                mzml_table[mzML_name]['Comet'] = identified_num
            else:
                mzml_table[mzML_name][search_engine] = identified_num

            mzml_table[mzML_name]['Final result of spectra'] = self.mL_spec_ident_final[mzML_name]
            mzml_table[mzML_name]['Final result of Peptides'] = len(self.mzml_peptide_map[mzML_name])

        # mzml file sorted based on experimental file
        for mzML_name in self.exp_design_table.keys():
            self.mzml_table[mzML_name] = mzml_table[mzML_name]

    def parse_out_mzTab(self):
        log.warning("Parsing mzTab file...")
        mztab_data = mztab.MzTab(self.out_mzTab_path)
        pep_table = mztab_data.peptide_table
        meta_data = dict(mztab_data.metadata)

        self.delta_mass['target'] = dict()
        self.delta_mass['decoy'] = dict()

        # PSM table data
        psm = mztab_data.spectrum_match_table
        prot = mztab_data.protein_table
        self.prot_search_score = dict()
        psm['stand_spectra_ref'] = psm.apply(
            lambda x: os.path.basename(meta_data[x.spectra_ref.split(':')[0] + '-location']), axis=1)

        pro_abundance = list(filter(lambda x: re.match(r'protein_abundance_assay.*?', x) is not None,
                                    prot.columns.tolist()))

        if config.kwargs['remove_decoy']:
            psm = psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] != 1]
            prot = prot[~prot['accession'].str.contains(config.kwargs['decoy_affix'])]
        prot = prot[prot['opt_global_result_type'] != 'protein_details']
        prot["protein_group"] = prot.apply(lambda x: x["ambiguity_members"].replace(",", ";"), axis=1)
        self.prot_search_score = dict(zip(prot["protein_group"], prot["best_search_engine_score[1]"]))
        self.Total_Protein_Identified = len(prot.index)
        prot.dropna(how='all',subset=pro_abundance, inplace=True)
        self.Total_Protein_Quantified = len(prot.index)

        self.pep_plot = Histogram('number of peptides per proteins', plot_category = 'frequency')

        for protein in prot['accession']:
            number = sum(pep_table.apply(lambda x: all(p in x['accession'] for p in protein.split(",")),axis=1))
            self.pep_plot.addValue(number)
        
        categorys = OrderedDict()
        categorys['Frequency'] = {
            'name': 'Frequency',
            'description': 'number of peptides per proteins'
        }
        self.pep_plot.to_dict(percentage = True, cats = categorys)

        mL_spec_ident_final = {}

        for m, group in psm.groupby('stand_spectra_ref'):
            self.cal_num_table_data[m] = {'Sample Name': self.exp_design_table[m]['Sample']}
            self.cal_num_table_data[m]['condition'] = self.exp_design_table[m]['MSstats_Condition']
            self.cal_num_table_data[m]['fraction'] = self.exp_design_table[m]['Fraction']

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
            self.cal_num_table_data[m]['protein_num'] = len(proteins)
            self.cal_num_table_data[m]['peptide_num'] = len(peptides)
            self.cal_num_table_data[m]['unique_peptide_num'] = len(unique_peptides)

            modified_pep = list(filter(lambda x: re.match(r'.*?\(.*\).*?', x) is not None, peptides))
            self.cal_num_table_data[m]['modified_peptide_num'] = len(modified_pep)

            mL_spec_ident_final[m] = len(set(self.identified_spectrum[m]))
        

        target_bin_data = {}
        decoy_bin_data = {}
        psm['relative_diff'] = psm['exp_mass_to_charge'] - psm['calc_mass_to_charge']
        try:
            decoy_bin = psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] == 1]['relative_diff'].value_counts(sort=False, bins=1000)
            for index in decoy_bin.index:
                decoy_bin_data[float(index.mid)] = int(decoy_bin[index])
            self.delta_mass['decoy'] = decoy_bin_data
        except Exception as e:
            print("No decoy peptides, only show target peptides")

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

        self.prot_search_score = dict()
        for prot, group in report_data.groupby("Protein.Names"):
            self.prot_search_score[prot] = 1 - np.min(group["Q.Value"])
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
        print("Parsing MSstats input file......")
        msstats_data = pd.read_csv(self.msstats_input_path)
        msstats_data["Run"] = msstats_data.apply(lambda x: os.path.splitext(os.path.basename(x["Reference"]))[0], axis=1)
        f_table = self.f_table.drop_duplicates(subset="Run", inplace=False)
        f_table = f_table[["Fraction_Group", "Run"]]
        msstats_data = msstats_data.merge(f_table, on="Run")
        msstats_data_dict = dict()
        msstats_data_pep  = dict()
        msstats_data_prot = dict()
        pep_intensity = dict()
        index = 1
        conditions = list(set(msstats_data["Condition"].tolist()))
        for pep, group in msstats_data.groupby("PeptideSequence"):
            if config.kwargs["remove_decoy"] and config.kwargs["decoy_affix"] in group["ProteinName"].values[0]:
                continue
            msstats_data_dict[index] = {'PeptideSequence': pep} 
            msstats_data_dict[index]['ProteinName'] = group["ProteinName"].values[0]

            # aggregate intensity of fraction and peptidoforms, then average technical/biological replicates intensity
            for condition, c_group in group.groupby("Condition"):
                condition_intensities = []
                if "Channel" in list(c_group):
                    for _, channel_group in c_group.groupby(["Channel","Fraction_Group"]):
                        condition_intensities.append(channel_group["Intensity"].sum())
                else:
                    for _, f_group in c_group.groupby("Fraction_Group"):
                        condition_intensities.append(f_group["Intensity"].sum())

                msstats_data_dict[index][str(condition)] = np.mean(condition_intensities)
            
                for prot in np.unique(group["ProteinName"]):
                    msstats_data_dict[index]["BestSearchScore"] = 1 - self.prot_search_score[prot]
                    if prot not in pep_intensity:
                        pep_intensity[prot] = {str(condition): [np.mean(condition_intensities)]}
                    elif str(condition) not in pep_intensity[prot]:
                        pep_intensity[prot][str(condition)] = [np.mean(condition_intensities)]
                    else:
                        pep_intensity[prot][str(condition)].append(np.mean(condition_intensities))
            
            rc = list(map(str,(set(group["Condition"]))))
            aic = itemgetter(*rc)(msstats_data_dict[index])
            try:
                msstats_data_dict[index]["Average Intensity"] = np.mean(list(aic))
            except Exception as e:
                msstats_data_dict[index]["Average Intensity"] = aic
            for c in conditions:
                if str(c) not in msstats_data_dict[index]:
                    msstats_data_dict[index][str(c)] = np.NaN
                    msstats_data_dict[index][str(c) + "_distribution"] = ""
                else:
                    gc = group[group["Condition"] == c]
                    intensity_distribution = (gc["Intensity"] - msstats_data_dict[index]["Average Intensity"]).tolist()
                    intensity_distribution_str = list(map(str, intensity_distribution))
                    msstats_data_dict[index][str(c) + "_distribution"] = ", ".join(intensity_distribution_str) + " ; column"
            
            if index < 51:
                msstats_data_pep[index] = msstats_data_dict[index]
            index += 1

        cur.execute("CREATE TABLE PEPQUANT(ID integer PRIMARY KEY AUTOINCREMENT, ProteinName VARCHAR(100), PeptideSequence VARCHAR(100), BestSearchScore FLOAT(4,3), \"Average Intensity\" FLOAT(4,3))")
        con.commit()
        sql_col = "ProteinName,PeptideSequence,BestSearchScore, \"Average Intensity\""
        sql_t = "(" + ','.join(['?'] * (len(conditions) *2 + 4)) + ")"

        headers = OrderedDict()
        headers = {'ProteinName': {'name': 'ProteinName'}, 'PeptideSequence': {'name': 'PeptideSequence'},
                    'BestSearchScore': {'name': 'BestSearchScore'}, 'Average Intensity': {'name': 'Average Intensity'}}
        
        for s in conditions:
            cur.execute("ALTER TABLE PEPQUANT ADD \"" + str(s) + "\" FLOAT")
            con.commit()
            sql_col += ", \"" + str(s) + "\""
            headers[str(s)] = {'name': s}

        for s in list(map(lambda x: str(x) + "_distribution", conditions)):
            cur.execute("ALTER TABLE PEPQUANT ADD \"" + s + "\" VARCHAR(100)")
            con.commit()
            sql_col += ", \"" + s + "\""
            headers[str(s)] = {'name': s}

        all_term = ["ProteinName", "PeptideSequence", "BestSearchScore", "Average Intensity"] + list(map(str, conditions)) + list(map(lambda x: str(x) + "_distribution", conditions))
        cur.executemany("INSERT INTO PEPQUANT (" + sql_col + ") VALUES " + sql_t,
                        [itemgetter(*all_term)(x[1]) for x in msstats_data_dict.items()])
        con.commit()


        pconfig = {
            'id': 'quantification_of_peptides',  # ID used for the table
            'table_title': 'quantification information of peptides',
            # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'raw_data_fn': 'multiqc_quantification_of_peptides_table',  # File basename to use for raw data file
            'sortRows': False,  # Whether to sort rows alphabetically
            'only_defined_headers': False,  # Only show columns that are defined in the headers config
            'col1_header': 'index',
            'format': '{:,.3f}',
            'no_beeswarm': True
        }

        table_html = sparklines.plot(msstats_data_pep, headers, pconfig=pconfig)
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
                            ''',
            plot=table_html
        )


        # plot protein table
        index = 1
        plot_prot = dict()
        for prot, group in msstats_data.groupby("ProteinName"):
            msstats_data_prot[prot] = {"Peptides_Number": len(np.unique(group["PeptideSequence"]))}
            for condition, c_group in group.groupby("Condition"):
                msstats_data_prot[prot][str(condition)] = np.sum(pep_intensity[prot][str(condition)])

            rc = list(map(str,(set(group["Condition"]))))
            aic = itemgetter(*rc)(msstats_data_prot[prot])
            try:
                msstats_data_prot[prot]["Average Intensity"] = np.mean(list(aic))
            except Exception as e:
                msstats_data_prot[prot]["Average Intensity"] = aic

            for c in conditions:
                if str(c) not in msstats_data_prot[prot]:
                    msstats_data_prot[prot][str(c)] = np.NaN
                    msstats_data_prot[prot][str(c) + "_distribution"] = ""
                else:
                    gc = group[group["Condition"] == c]
                    intensity_distribution_str = []
                    if "Channel" in list(msstats_data):
                        for _, ref_group in gc.groupby("Run"):
                            channel_intensities = []
                            for _, prot_channel_group in ref_group.groupby("Channel"):
                                channel_intensities.append(np.sum(prot_channel_group["Intensity"]))
                            intensity_distribution_str.append(str(np.mean(channel_intensities) - msstats_data_prot[prot]["Average Intensity"]))
                    else:
                        for _, ref_group in gc.groupby("Run"):
                            intensity_distribution =(np.sum(ref_group["Intensity"]) - msstats_data_prot[prot]["Average Intensity"]).tolist()
                            intensity_distribution_str.append(str(intensity_distribution))
                    msstats_data_prot[prot][str(c) + "_distribution"] = ", ".join(intensity_distribution_str) + " ; column"

            if index < 51:
                plot_prot[prot] = msstats_data_prot[prot]
            index += 1


        headers = OrderedDict()
        headers['Peptides_Number'] = {
            'name': 'Number Peptides',
            'description': 'number of peptides per proteins',
            'format': '{:,.0f}'
        }
        headers['Average Intensity'] = {
            'name': 'Average Intensity',
            'description': 'average intensity across all condition'
        }

        # upload protein table to sqlite database
        cur.execute("CREATE TABLE PROTQUANT(ProteinName VARCHAR(100), Peptides_Number INT(100), \"Average Intensity\" FLOAT(4,3))")
        con.commit()
        sql_col = "ProteinName,Peptides_Number,\"Average Intensity\""
        sql_t = "(" + ','.join(['?'] * (len(conditions)*2 + 3)) + ")"

        for s in conditions:
            cur.execute("ALTER TABLE PROTQUANT ADD \"" + str(s) + "\" FLOAT")
            con.commit()
            sql_col += ", \"" + str(s) + "\""
            headers[str(s)] = {'name': s}

        for s in list(map(lambda x: str(x) + "_distribution", conditions)):
            cur.execute("ALTER TABLE PROTQUANT ADD \"" + s + "\" VARCHAR(100)")
            con.commit()
            sql_col += ", \"" + s + "\""
            headers[str(s)] = {'name': s}

        all_term = ["Peptides_Number", "Average Intensity"] + list(map(str, conditions)) + list(map(lambda x: str(x) + "_distribution", conditions))
        cur.executemany("INSERT INTO PROTQUANT (" + sql_col + ") VALUES " + sql_t,
                        [tuple([x[0]]) + itemgetter(*all_term)(x[1]) for x in msstats_data_prot.items()])
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
            'format': '{:,.3f}',
            'no_beeswarm': True
        }

        table_html = sparklines.plot(plot_prot, headers, pconfig=pconfig)
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
                        'in quantms pipeline final result',
            helptext='''
                        The quantification information of proteins is obtained from the msstats input file. 
                        The table shows the quantitative level and distribution of proteins in different study variables and run.
                            ''',
            plot=table_html
        )

