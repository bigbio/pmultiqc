#!/usr/bin/env python

""" MultiQC example plugin module """

from __future__ import absolute_import
from collections import OrderedDict
import logging
from sdrf_pipelines.openms.openms import OpenMS
from multiqc import config
from multiqc.plots import table, bargraph, linegraph, heatmap
from multiqc.modules.base_module import BaseMultiqcModule
import pandas as pd
from collections import Counter
import re
from pyteomics import mztab, mzml, openms
import os
import sqlite3
import numpy as np
import math

# Initialise the main MultiQC logger
log = logging.getLogger(__name__)
if config.output_dir:
    if os.path.exists(config.output_dir):
        con = sqlite3.connect(os.path.join(config.output_dir, 'proteomicslfq.db'))
    else:
        os.makedirs(config.output_dir)
        con = sqlite3.connect(os.path.join(config.output_dir, 'proteomicslfq.db'))
else:
    con = sqlite3.connect('./proteomicslfq.db')

cur = con.cursor()
cur.execute("drop table if exists PSM")
con.commit()

cur.execute("drop table if exists QUANT")
con.commit()


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_plugin', True):
            return None

        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name='ProteomicsLFQ',
            target="proteomicslfq",
            anchor='proteomicslfq',
            href='https://github.com/MultiQC/example-plugin',
            info=" is an module to show the pipeline performance."
        )

        for file in self.find_log_files({'fn': "*.csv"}):
            self.out_csv_path = os.path.join(config.analysis_dir[0], file['fn'])

        for f in self.find_log_files({'fn': 'out.mzTab'}):
            self.out_mzTab_path = os.path.join(config.analysis_dir[0], f['fn'])

        if config.kwargs['exp_design']:
            self.exp_design = config.kwargs['exp_design']
        elif config.kwargs['sdrf']:
            OpenMS().openms_convert(config.kwargs['sdrf'], config.kwargs['raw'],
                                    False, True, False, config.kwargs['condition'])
            self.exp_design = './' + 'experimental_design.tsv'

        else:
            raise AttributeError("exp_design and sdrf cannot be empty at the same time")

        self.PSM_table = dict()
        self.mzml_peptide_map = dict()
        self.identified_spectrum = dict()
        self.pep_quant_table = dict()
        self.mzml_table = dict()
        self.Total_ms2_Spectral = 0
        self.Total_ms2_Spectral_Identified = 0
        self.Total_Peptide_Count = 0
        self.Total_Protein_Identified = 0
        self.num_pep_per_protein = OrderedDict()
        self.percent_pep_per_protein = OrderedDict()
        self.out_csv_data = dict()
        self.cal_num_table_data = dict()
        self.mL_spec_ident_final = dict()
        self.peak_intensity_distribution = dict()
        self.charge_state_distribution = dict()
        self.peak_per_ms2 = dict()
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

        # parse input data
        self.CalHeatMapScore()
        if config.kwargs['quant_method'] == 'lfq':
            self.parse_out_csv()
            self.parse_out_mzTab()
        else:
            self.parse_out_tmt()

        self.parse_mzml_idx()

        self.draw_heatmap()
        # draw the experimental design
        self.draw_exp_design()

        self.draw_summary_proten_ident_table()

        # draw_proteomicslfq_identi_num
        self.draw_proteomicslfq_identi_num()

        # draw number of peptides per protein
        self.draw_num_pep_per_protein()

        # draw peptides quantification information

        self.draw_pep_quant_info()

        self.draw_psm_table()
        self.draw_mzml_ms()
        self.draw_precursor_charge_distribution()
        self.draw_peaks_per_ms2()
        self.draw_peak_intensity_distribution()
        self.draw_delta_mass()

        self.css = {
            'assets/css/proteomicslfq.css':
                os.path.join(os.path.dirname(__file__), 'assets', 'css', 'proteomicslfq.css')
        }
        self.js = {
            'assets/js/proteomicslfq.js':
                os.path.join(os.path.dirname(__file__), 'assets', 'js', 'proteomicslfq.js'),
            'assets/js/axios.min.js':
                os.path.join(os.path.dirname(__file__), 'assets', 'js', 'axios.min.js'),
            'assets/js/sql-optimized.js':
                os.path.join(os.path.dirname(__file__), 'assets', 'js', 'sql-optimized.js')
        }

    def draw_heatmap(self):
        HeatMapScore = []
        xnames = ['Contaminants', 'Peptide Intensity', 'Charge', 'Missed Cleavages', 'Missed Cleavages Var',
                  'ID rate over RT', 'MS2 OverSampling', 'Pep Missing Values']
        ynames = []
        for k, v in self.heatmap_con_score.items():
            ynames.append(k)
            HeatMapScore.append([v, self.heatmap_pep_intensity[k], self.heatmap_charge_score[k],
                                 self.MissedCleavages_heatmap_score[k],
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
            anchor="proteomicslfq_heatmap",
            description='This plot shows the pipeline performance overview',
            helptext='''
                    This plot shows the pipeline performance overview. Some metrics are calculated.
                    
                    * Heatmap score[Contaminants]: as fraction of summed intensity with 0 = sample full of contaminants; 
                        1 = no contaminants
                    * Heatmap score[Pep Intensity (>23.0)]: Linear scale of the median intensity reaching the threshold, 
                        i.e. reaching 221 of 223 gives score 0.25.
                    * Heatmap score[Charge]: Deviation of the charge 2 proportion from a representative Raw file.
                    * Heatmap score [MC]: the fraction (0% - 100%) of fully cleaved peptides per Raw file
                    * Heatmap score [MC Var]: each Raw file is scored for its deviation from the ‘average’ digestion 
                        state of the current study.
                    * Heatmap score [ID rate over RT]: Scored using ‘Uniform’ scoring function. 
                    * Heatmap score [MS2 Oversampling]: The percentage of non-oversampled 3D-peaks.
                    * Heatmap score [Pep Missing]: Linear scale of the fraction of missing peptides.
                    
                    ''',
            plot=hm_html
        )

    def draw_exp_design(self):
        exp_design_table = dict()
        with open(self.exp_design, 'r') as f:
            data = f.readlines()
            empty_row = data.index('\n')
            f_table = [i.replace('\n', '').split('\t') for i in data[1:empty_row]]
            f_table = pd.DataFrame(f_table)
            s_table = [i.replace('\n', '').split('\t') for i in data[data.index('\n') + 1:]][1:]
            s_header = data[data.index('\n') + 1].replace('\n', '').split('\t')
            s_DataFrame = pd.DataFrame(s_table, columns=s_header)
            for file in np.unique(f_table[2].tolist()):
                exp_design_table[file] = {'Fraction_Group': f_table[f_table[2] == file][0].tolist()[0]}
                exp_design_table[file]['Fraction'] = f_table[f_table[2] == file][1].tolist()[0]
                exp_design_table[file]['Label'] = ','.join(f_table[f_table[2] == file][3])
                sample = f_table[f_table[2] == file][4].tolist()
                exp_design_table[file]['Sample'] = ','.join(sample)
                exp_design_table[file]['MSstats_Condition'] = ','.join(
                    [row['MSstats_Condition'] for _, row in s_DataFrame.iterrows()
                     if row['Sample'] in sample])
                exp_design_table[file]['MSstats_BioReplicate'] = ','.join(
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

        headers['Fraction_Group'] = {
            'description': 'Fraction_Group',
            'color': "#ffffff",
        }
        headers['Fraction'] = {
            'description': 'Fraction identifier',
            'color': "#ffffff",
        }
        headers['Label'] = {
            'description': 'Label',
            'color': "#ffffff",
        }
        headers['Sample'] = {
            'description': 'Sample Name',
            'color': "#ffffff",
        }
        headers['MSstats_Condition'] = {
            'description': 'MSstats Condition',
            'color': "#ffffff",
        }
        headers['MSstats_BioReplicate'] = {
            'description': 'MSstats BioReplicate',
            'color': "#ffffff",
        }
        table_html = table.plot(exp_design_table, headers, pconfig)

        # Add a report section with the line plot
        self.add_section(
            name="Experimental Design",
            anchor="proteomicslfq_exp_design",
            description='This plot shows the Proteomics Experimental Design',
            helptext='''
            This plot shows the Proteomics Experimental Design. You can see details about it in 
            https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/classOpenMS_1_1ExperimentalDesign.html
            ''',
            plot=table_html
        )

    def draw_summary_proten_ident_table(self):
        summary_table = {
            self.Total_ms2_Spectral: {"Total MS/MS Spectral Identified": self.Total_ms2_Spectral_Identified}}
        coverage = self.Total_ms2_Spectral_Identified / self.Total_ms2_Spectral * 100
        summary_table[self.Total_ms2_Spectral]['Identified MS/MS Spectral Coverage'] = coverage
        summary_table[self.Total_ms2_Spectral]['Total Peptide Count'] = self.Total_Peptide_Count
        summary_table[self.Total_ms2_Spectral]['Total Protein Identified'] = self.Total_Protein_Identified

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
        headers = OrderedDict()

        headers['Total MS/MS Spectral Identified'] = {
            'description': 'Total MS/MS Spectral Identified',
        }
        headers['Identified MS/MS Spectral Coverage'] = {
            'description': 'Identified MS/MS Spectral Coverage',
            'format': '{:,.2f}',
            "suffix": "%"
        }

        table_html = table.plot(summary_table, headers, pconfig)

        # Add a report section with the line plot
        self.add_section(
            name="Summary Table",
            anchor="proteomicslfq_summary_table",
            description='This plot shows the ProteomicsLFQ pipeline summary statistics',
            helptext='''
                    This plot shows the ProteomicsLFQ pipeline summary statistics
                   ''',
            plot=table_html
        )

    def draw_proteomicslfq_identi_num(self):
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
        if config.kwargs['quant_method'] == 'lfq':
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
            description='This plot shows the ProteomicsLFQ pipeline final result',
            helptext='''
            This plot shows the ProteomicsLFQ pipeline final result.
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
            'data_labels': [
                {'name': 'count', 'ylab': 'Frequency'},
                {'name': 'percentage', 'ylab': 'Percentage', 'format': '{:,.2f}', 'suffix': '%'}
            ]
        }
        headers = OrderedDict()
        headers['Frequency'] = {
            'name': 'Frequency',
            'description': 'number of peptides per proteins'
        }

        bar_html = bargraph.plot([self.num_pep_per_protein, self.percent_pep_per_protein], headers, pconfig)
        # Add a report section with the line plot
        self.add_section(
            name="Number of Peptides Per Proteins",
            anchor="num_of_pep_per_prot",
            description='This plot shows the number of peptides per proteins '
                        'in ProteomicsLFQ pipeline final result',
            helptext='''
                        This statistic is extracted from the out_msstats file. Proteins supported by more peptide 
                        identifications can constitute more confident results.
                    ''',
            plot=bar_html
        )

    def draw_pep_quant_info(self):
        pconfig = {
            'id': 'quantification_of_peptides',  # ID used for the table
            'table_title': 'quantification information of peptides',
            # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'raw_data_fn': 'multiqc_quantification_of_peptides_table',  # File basename to use for raw data file
            'sortRows': False,  # Whether to sort rows alphabetically
            'only_defined_headers': False,  # Only show columns that are defined in the headers config
            'col1_header': 'ID',
            'format': '{:,.6f}',  # The header used for the first column
            'scale': 'Set3',
            'no_beeswarm': True
        }

        headers = OrderedDict()
        headers['sequence'] = {
            'description': 'peptide_sequence',
            'color': "#ffffff",
            'format': '{:,.0f}'
        }
        table_html = table.plot(self.pep_quant_table, headers, pconfig)
        pattern = re.compile(r'<small id="quantification_of_peptides_numrows_text"')
        index = re.search(pattern, table_html).span()[0]
        t_html = table_html[:index] + '<input type="text" placeholder="search..." class="searchInput" ' \
                                      'onkeyup="searchQuantFunction()" id="quant_search">' \
                                      '<select name="quant_search_col" id="quant_search_col"><option value="">' \
                                      'ID</option>'
        for key, _ in self.pep_quant_table['1'].items():
            t_html += '<option>' + key.replace("[", "_").replace("]", "") + '</option>'
        table_html = t_html + '</select>' + '<button type="button" class="btn btn-default ' \
                                            'btn-sm" id="quant_reset" onclick="quantFirst()">Reset</button>' \
                     + table_html[index:]
        table_html = table_html + '''<div class="page_control"><span id="quantFirst">First Page</span><span 
        id="quantPre"> Previous Page</span><span id="quantNext">Next Page </span><span id="quantLast">Last 
        Page</span><span id="quantPageNum"></span>Page/Total <span id="quantTotalPage"></span>Pages <input 
        type="number" name="" id="pep_page" class="page" value="" oninput="this.value=this.value.replace(/\D/g);" 
        onkeydown="quant_page_jump()" min="1"/> </div> '''

        self.add_section(
            name="Quantification Result",
            anchor="quant_result",
            description='This plot shows the quantification information of peptides'
                        'in ProteomicsLFQ pipeline final result',
            helptext='''
                        The quantification information of peptides is obtained from the pep table in the mzTab file. 
                        The table shows the quantitative level of peptides in different study variables.
                            ''',
            plot=table_html
        )

    def draw_psm_table(self):
        pconfig = {
            'id': 'peptide_spectrum_match',  # ID used for the table
            'table_title': 'peptide spectrum match information',
            # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'raw_data_fn': 'multiqc_psm_table',  # File basename to use for raw data file
            'sortRows': False,  # Whether to sort rows alphabetically
            'only_defined_headers': False,  # Only show columns that are defined in the headers config
            'col1_header': 'PSM_ID',
            'format': '{:,.9f}',  # The header used for the first column
            'scale': 'Set5',
            'no_beeswarm': True
        }

        headers = OrderedDict()
        headers['sequence'] = {
            'description': 'peptide_sequence',
            'color': "#ffffff",
            'format': '{:,.0f}'
        }

        headers['unique'] = {
            'description': 'unique',
            'color': "#ffffff",
            'format': '{:,.0f}'
        }

        table_html = table.plot(self.PSM_table, headers, pconfig)

        pattern = re.compile(r'<small id="peptide_spectrum_match_numrows_text"')
        index = re.search(pattern, table_html).span()[0]
        t_html = table_html[:index] + '<input type="text" placeholder="search..." class="searchInput" ' \
                                      'onkeyup="searchPsmFunction()" id="psm_search">' \
                                      '<select name="psm_search_col" id="psm_search_col"><option value="">' \
                                      'PSM_ID</option>'

        for key, _ in self.PSM_table[1].items():
            t_html += '<option>' + key.replace("[", "_").replace("]", "") + '</option>'
        table_html = t_html + '</select>' + '<button type="button" class="btn btn-default ' \
                                            'btn-sm" id="psm_reset" onclick="psmFirst()">Reset</button>' \
                     + table_html[index:]
        table_html = table_html + '''<div class="page_control"><span id="psmFirst">First Page</span><span 
                id="psmPre"> Previous Page</span><span id="psmNext">Next Page </span><span id="psmLast">Last 
                Page</span><span id="psmPageNum"></span>Page/Total <span id="psmTotalPage"></span>Pages <input 
                type="number" name="" id="psm_page" class="page" value="" oninput="this.value=this.value.replace(/\D/g);" 
                onkeydown="psm_page_jump()" min="1"/> </div> '''

        # Add a report section with the line plot
        self.add_section(
            name="Peptide-Spectrum Matches",
            anchor="psm",
            description='This plot shows the PSM information'
                        'in ProteomicsLFQ pipeline final result',
            helptext='''
                        This table fully displays the peptide spectrum matching information in the mzTab file:

                        * sequence: peptide sequence
                        * unique
                        * search_engine_score
                            ''',
            plot=table_html
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
            description='This plot shows the ProteomicsLFQ pipeline mzML tracking',
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
        cats = OrderedDict()
        cats['0-10'] = {
            'name': '0-10',
            'description': 'Peak intensity is between 0 and 10'
        }
        cats['10-100'] = {
            'name': '10-100',
            'description': 'Peak intensity is between 10 and 100'
        }
        cats['100-300'] = {
            'name': '100-300',
            'description': 'Peak intensity is between 100 and 300'
        }
        cats['300-500'] = {
            'name': '300-500',
            'description': 'Peak intensity is between 300 and 500'
        }
        cats['500-700'] = {
            'name': '500-700',
            'description': 'Peak intensity is between 500 and 700'
        }
        cats['700-900'] = {
            'name': '700-900',
            'description': 'Peak intensity is between 700 and 900'
        }
        cats['900-1000'] = {
            'name': '900-1000',
            'description': 'Peak intensity is between 900 and 1000'
        }
        cats['1000-3000'] = {
            'name': '1000-3000',
            'description': 'Peak intensity is between 1000 and 3000'
        }
        cats['3000-6000'] = {
            'name': '3000-6000',
            'description': 'Peak intensity is between 3000 and 6000'
        }
        cats['6000-10000'] = {
            'name': '6000-10000',
            'description': 'Peak intensity is between 6000 and 10000'
        }
        cats['>10000'] = {
            'name': '>10000',
            'description': 'peak intensity > 10000'
        }

        bar_html = bargraph.plot(self.peak_intensity_distribution, cats, pconfig)
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
        cats = OrderedDict()
        cats['1'] = {
            'name': '1',
            'description': 'Precursor charge state is 1'
        }
        cats['2'] = {
            'name': '2',
            'description': 'Precursor charge state is 2'
        }
        cats['3'] = {
            'name': '3',
            'description': 'Precursor charge state is 3'
        }
        cats['4'] = {
            'name': '4',
            'description': 'Precursor charge state is 4'
        }
        cats['5'] = {
            'name': '5',
            'description': 'Precursor charge state is 5'
        }
        cats['6'] = {
            'name': '6',
            'description': 'Precursor charge state is 6'
        }
        cats['7'] = {
            'name': '7',
            'description': 'Precursor charge state is 7'
        }
        cats['>7'] = {
            'name': '>7',
            'description': 'Precursor charge state >7'
        }

        bar_html = bargraph.plot(self.charge_state_distribution, cats, pconfig)

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
        cats = OrderedDict()
        cats['0-100'] = {
            'name': '0-100',
            'description': 'Number of Peaks per MS/MS spectrum is between 0 and 100'
        }
        cats['100-200'] = {
            'name': '100-200',
            'description': 'Number of Peaks per MS/MS spectrum is between 100 and 200'
        }
        cats['200-300'] = {
            'name': '200-300',
            'description': 'Number of Peaks per MS/MS spectrum is between 200 and 300'
        }
        cats['300-400'] = {
            'name': '300-400',
            'description': 'Number of Peaks per MS/MS spectrum is between 300 and 400'
        }
        cats['400-500'] = {
            'name': '400-500',
            'description': 'Number of Peaks per MS/MS spectrum is between 400 and 500'
        }
        cats['500-600'] = {
            'name': '500-600',
            'description': 'Number of Peaks per MS/MS spectrum is between 500 and 600'
        }
        cats['600-700'] = {
            'name': '600-700',
            'description': 'Number of Peaks per MS/MS spectrum is between 600 and 700'
        }
        cats['700-800'] = {
            'name': '700-800',
            'description': 'Number of Peaks per MS/MS spectrum is between 700 and 800'
        }
        cats['800-900'] = {
            'name': '800-900',
            'description': 'Number of Peaks per MS/MS spectrum is between 800 and 900'
        }
        cats['900-1000'] = {
            'name': '900-1000',
            'description': 'Number of Peaks per MS/MS spectrum is between 900 and 1000'
        }
        cats['>1000'] = {
            'name': '>1000',
            'description': 'Number of Peaks per MS/MS spectrum > 1000'
        }

        bar_html = bargraph.plot(self.peak_per_ms2, cats, pconfig)

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

    def draw_delta_mass(self):

        self.delta_mass_percent['target'] = dict(zip(self.delta_mass['target'].keys(),
                                                     list(map(lambda v: v / len(self.delta_mass['target']),
                                                              self.delta_mass['target'].values()))))
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
        pep_table = mztab_data.peptide_table
        meta_data = dict(mztab_data.metadata)
        pep_table.loc[:, 'stand_spectra_ref'] = pep_table.apply(
            lambda x: os.path.basename(meta_data[x.spectra_ref.split(':')[0] + '-location']), axis=1)
        study_variables = list(filter(lambda x: re.match(r'peptide_abundance_study_variable.*?', x) is not None,
                                      pep_table.columns.tolist()))
        # total_intensity = np.sum(pep_table[study_variables])
        # cont_top = {} if set(pep_table[pep_table['opt_global_cv_MS:1002217_decoy_peptide'] == 1]['accession']) == {
        # None}: for i in set(pep_table[pep_table['opt_global_cv_MS:1002217_decoy_peptide'] == 1]['sequence']):
        # cont_top[i] = np.sum(np.sum(pep_table[pep_table['sequence'] == i][study_variables])) / total_intensity
        # else: for i in set(pep_table[pep_table['opt_global_cv_MS:1002217_decoy_peptide'] == 1]['accession']):
        # cont_top[i] = np.sum(np.sum(pep_table[pep_table['accession'] == i][study_variables]) / total_intensity)
        # cont_top_5 = sorted(cont_top.items(), key=lambda item: item[1], reverse=True) dictdata = {} for l in
        # cont_top_5: dictdata[l[0]] = l[1]
        for i in np.unique(pep_table['stand_spectra_ref']):
            self.heatmap_con_score[i] = 1 - np.sum(np.sum(
                pep_table[(pep_table['stand_spectra_ref'] == i) &
                          (pep_table['opt_global_cv_MS:1002217_decoy_peptide'] == 1)][study_variables])) / \
                                        np.sum(np.sum(pep_table[pep_table['stand_spectra_ref'] == i]
                                                      [study_variables]))
            if config.kwargs['remove_decoy']:
                T = sum(pep_table[(pep_table['stand_spectra_ref'] == i)
                                  & (pep_table['opt_global_cv_MS:1002217_decoy_peptide'] == 0)][study_variables].values. \
                        tolist(), [])
            else:
                T = sum(pep_table[pep_table['stand_spectra_ref'] == i][study_variables].values.tolist(), [])
            pep_median = np.median([j for j in T if not math.isnan(j) is True])
            self.heatmap_pep_intensity[i] = np.minimum(1.0, pep_median / (2 ** 23))  # Threshold

        #  HeatMapMissedCleavages
        global_peps = set(psm['opt_global_cv_MS:1000889_peptidoform_sequence'])
        global_peps_count = len(global_peps)
        if config.kwargs['remove_decoy'] and 'opt_global_cv_MS:1002217_decoy_peptide' in psm.columns:
            psm = psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] == 0]
        psm.loc[:, 'stand_spectra_ref'] = psm.apply(
            lambda x: os.path.basename(meta_data[x.spectra_ref.split(':')[0] + '-location']), axis=1)
        psm.loc[:, 'missed_cleavages'] = psm.apply(lambda x: self.cal_MissCleavages(x['sequence']), axis=1)
        msms_count_map = {}
        for name, group in psm.groupby(
                ['stand_spectra_ref', 'opt_global_cv_MS:1000889_peptidoform_sequence', 'charge']):
            if name[0] not in msms_count_map:
                msms_count_map[name[0]] = [len(group)]
            else:
                msms_count_map[name[0]].append(len(group))

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
            if np.max(msms_count_map[name]) >= 3:
                for i, value in enumerate(msms_count_map[name]):
                    if value >= 3:
                        msms_count_map[name][i] = "3+"
            tt = Counter(msms_count_map[name])
            self.HeatmapOverSamplingScore[name] = np.minimum(1.0, tt[1] * 1.0 / np.sum(list(tt.values())))

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

    def parse_out_csv(self):

        # result statistics table
        log.warning("Parsing out csv file...")
        data = pd.read_csv(self.out_csv_path, sep=',', header=0)
        data = data.astype(str)
        exp_data = pd.read_csv(self.exp_design, sep='\t', header=0, index_col=None, dtype=str)
        exp_data = exp_data.dropna(axis=0)
        Spec_File = exp_data['Spectra_Filepath'].tolist()

        if config.kwargs['remove_decoy']:
            data['DECOY'] = data.apply(lambda x: self.dis_decoy(x['ProteinName']), axis=1)

        for i in Spec_File:
            proteins = []
            Sample = list(exp_data[exp_data['Spectra_Filepath'] == i]['Sample'])[0]
            self.cal_num_table_data[i] = {'Sample Name': Sample}
            condition = '-'.join(list(set(data[data['Reference'] == i]['Condition'].tolist())))
            self.cal_num_table_data[i]['condition'] = condition
            fraction = exp_data[exp_data['Spectra_Filepath'] == i]['Fraction'].tolist()[0]
            self.cal_num_table_data[i]['fraction'] = fraction
            data[data['Reference'] == i]['ProteinName'].apply(lambda x: proteins.extend(x.split(';')))

            if config.kwargs['remove_decoy']:
                if config.kwargs['affix_type'] == 'prefix':
                    proteins = list(filter(lambda x: not x.startswith(config.kwargs['decoy_affix']), proteins))
                else:
                    proteins = list(filter(lambda x: not x.endswith(config.kwargs['decoy_affix']), proteins))
                peptides = set(data[(data['Reference'] == i) & (data['DECOY'] == 'TARGET')]['PeptideSequence'].tolist())
            else:
                peptides = set(data[data['Reference'] == i]['PeptideSequence'].tolist())

            self.cal_num_table_data[i]['protein_num'] = len(set(proteins))
            self.cal_num_table_data[i]['peptide_num'] = len(peptides)

            modified_pep = list(filter(lambda x: re.match(r'.*?\(.*\).*?', x) is not None, peptides))
            self.cal_num_table_data[i]['modified_peptide_num'] = len(modified_pep)

        #  table of peptide number per protein
        prot_pep_map = {}
        num_pep_per_protein = {}
        for _, row in data.iterrows():
            ProteinName = row['ProteinName']
            if len(ProteinName.split(';')) == 1:
                if config.kwargs['remove_decoy']:
                    if config.kwargs['affix_type'] == 'prefix':
                        if ProteinName.startswith(config.kwargs['decoy_affix']):
                            continue
                    else:
                        if ProteinName.endswith(config.kwargs['decoy_affix']):
                            continue
                if ProteinName not in prot_pep_map.keys():
                    prot_pep_map[ProteinName] = [row['PeptideSequence']]
                elif row['PeptideSequence'] not in prot_pep_map[ProteinName]:
                    prot_pep_map[ProteinName].append(row['PeptideSequence'])
            else:
                for p in ProteinName.split(';'):
                    if config.kwargs['remove_decoy']:
                        if config.kwargs['affix_type'] == 'prefix':
                            if ProteinName.startswith(config.kwargs['decoy_affix']):
                                continue
                        else:
                            if ProteinName.endswith(config.kwargs['decoy_affix']):
                                continue
                    if p not in prot_pep_map.keys():
                        prot_pep_map[p] = [row['PeptideSequence']]
                    elif row['PeptideSequence'] not in prot_pep_map[p]:
                        prot_pep_map[p].append(row['PeptideSequence'])

        for key, value in prot_pep_map.items():
            if str(len(value)) not in num_pep_per_protein.keys():
                num_pep_per_protein[str(len(value))] = {'Frequency': 1}
            else:
                num_pep_per_protein[str(len(value))]['Frequency'] += 1

        keys = sorted(num_pep_per_protein.items(), key=lambda x: int(x[0]))
        for k in keys:
            self.num_pep_per_protein[k[0]] = k[1]
            self.percent_pep_per_protein[k[0]] = {'Percentage': k[1]['Frequency'] / len(prot_pep_map)}
        self.Total_Protein_Identified = len(prot_pep_map)

    def parse_out_mzTab(self):
        log.warning("Parsing mzTab file...")
        mztab_data = mztab.MzTab(self.out_mzTab_path)
        pep_table = mztab_data.peptide_table
        meta_data = dict(mztab_data.metadata)
        self.delta_mass['target'] = dict()
        self.delta_mass['decoy'] = dict()
        # PSM table data
        psm = mztab_data.spectrum_match_table
        psm_table = dict()
        scores = list(filter(lambda x: re.match(r'search_engine_score.*?', x) is not None,
                             psm.columns.tolist()))
        t = (' float'.join(scores) + ' float').replace("[", "_").replace("]", "")
        cur.execute("CREATE TABLE PSM(PSM_ID int, sequence VARCHAR(100), PSM_UNIQUE int, " + t + ")")
        con.commit()

        mL_spec_ident_final = {}
        psm['stand_spectra_ref'] = psm.apply(
            lambda x: os.path.basename(meta_data[x.spectra_ref.split(':')[0] + '-location']), axis=1)
        for m in set(psm['stand_spectra_ref'].tolist()):
            if config.kwargs['remove_decoy']:
                self.identified_spectrum[m] = list(map(lambda x: x.split(':')[1],
                                                       psm[(psm['stand_spectra_ref'] == m) & (
                                                               psm['opt_global_cv_MS:1002217_decoy_peptide'] != 1)][
                                                           'spectra_ref'].tolist()))
                self.mzml_peptide_map[m] = list(set(psm[(psm['stand_spectra_ref'] == m) &
                                                        (psm['opt_global_cv_MS:1002217_decoy_peptide'] != 1)][
                                                        'sequence'].tolist()))
            else:
                self.identified_spectrum[m] = list(map(lambda x: x.split(':')[1],
                                                       psm[psm['stand_spectra_ref'] == m]['spectra_ref'].tolist()))
                self.mzml_peptide_map[m] = list(set(psm[psm['stand_spectra_ref'] == m]['sequence'].tolist()))
            mL_spec_ident_final[m] = len(self.identified_spectrum[m])
        psm['relative_diff'] = psm['exp_mass_to_charge'] - psm['calc_mass_to_charge']
        self.delta_mass['decoy'] = Counter(psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] == 1]
                                           ['relative_diff'].tolist())
        self.delta_mass['target'] = Counter(psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] != 1]
                                            ['relative_diff'].tolist())

        sql_t = "(" + ','.join(['?'] * (len(scores) + 3)) + ")"
        if config.kwargs['remove_decoy']:
            cur.executemany("INSERT INTO PSM (PSM_ID,sequence,PSM_UNIQUE," +
                            ','.join(scores).replace("[", "_").replace("]", "") + ") VALUES " + sql_t,
                            [tuple(x) for x in psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] != 1]
                            [['PSM_ID', 'sequence', 'unique'] + scores].values])
        else:
            cur.executemany("INSERT INTO PSM (PSM_ID,sequence,PSM_UNIQUE," +
                            ','.join(scores).replace("[", "_").replace("]", "") + ") VALUES " + sql_t,
                            [tuple(x) for x in psm[['PSM_ID', 'sequence', 'unique'] + scores].values])
        con.commit()
        for idx, row in psm.iterrows():
            if config.kwargs['remove_decoy'] and row['opt_global_cv_MS:1002217_decoy_peptide'] == 1:
                continue
            else:
                psm_table[row['PSM_ID']] = {'sequence': row['sequence']}
                # psm_table[row['PSM_ID']]['spectra_ref'] = row['spectra_ref']
                psm_table[row['PSM_ID']]['unique'] = row['unique']
                for score in scores:
                    psm_table[row['PSM_ID']][score] = row[score]
            if idx > 48:
                break

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
        self.PSM_table = psm_table

        # peptide quantification table data
        pep_quant = dict()
        study_variables = list(filter(lambda x: re.match(r'peptide_abundance_study_variable.*?', x) is not None,
                                      pep_table.columns.tolist()))
        cur.execute("CREATE TABLE QUANT(ID integer PRIMARY KEY AUTOINCREMENT, sequence VARCHAR(100))")
        con.commit()
        sql_col = "sequence"
        sql_t = "(" + ','.join(['?'] * (len(study_variables) + 1)) + ")"

        for s in study_variables:
            s = s.replace("[", "_").replace("]", "")
            cur.execute("ALTER TABLE QUANT ADD " + s + " FLOAT")
            con.commit()
            sql_col += "," + s
        cur.executemany("INSERT INTO QUANT (" + sql_col + ") VALUES " + sql_t,
                        [tuple(x) for x in pep_table[['sequence'] + study_variables].values])
        con.commit()
        for index, row in pep_table.iterrows():
            if config.kwargs['remove_decoy'] and row['opt_global_cv_MS:1002217_decoy_peptide'] == 1:
                continue
            else:
                pep_quant[str(index + 1)] = {'sequence': row['sequence']}  # keep id unique
                for s in study_variables:
                    pep_quant[str(index + 1)][s] = row[s]
            if index > 48:
                break

        self.pep_quant_table = pep_quant

    def parse_mzml_idx(self):
        self.peak_intensity_distribution['identified_spectra'] = dict(zip(
            ['>10000', '6000-10000', '3000-6000', '1000-3000', '900-1000', '700-900',
             '500-700', '300-500', '100-300', '10-100', '0-10'], [0] * 11
        ))
        self.peak_intensity_distribution['unidentified_spectra'] = dict(zip(
            ['>10000', '6000-10000', '3000-6000', '1000-3000', '900-1000', '700-900',
             '500-700', '300-500', '100-300', '10-100', '0-10'], [0] * 11
        ))
        self.charge_state_distribution['identified_spectra'] = dict(zip(
            ['>7', '7', '6', '5', '4', '3', '2', '1'], [0] * 8
        ))
        self.charge_state_distribution['unidentified_spectra'] = dict(zip(
            ['>7', '7', '6', '5', '4', '3', '2', '1'], [0] * 8
        ))
        self.peak_per_ms2['identified_spectra'] = dict(zip(
            ['>1000', '900-1000', '800-900', '700-800', '600-700', '500-600', '400-500',
             '300-400', '200-300', '100-200', '0-100'], [0] * 11
        ))

        self.peak_per_ms2['unidentified_spectra'] = dict(zip(
            ['>1000', '900-1000', '800-900', '700-800', '600-700', '500-600', '400-500',
             '300-400', '200-300', '100-200', '0-100'], [0] * 11
        ))

        mzml_table = {}
        mzmls_dir = config.kwargs['mzMLs']
        heatmap_charge = {}
        for m in os.listdir(mzmls_dir):
            # print(os.path.join(mzmls_dir, m))
            log.warning("Parsing {}...".format(os.path.join(mzmls_dir, m)))
            ms1_number = 0
            ms2_number = 0
            total_charge_num = 0
            charge_2 = 0
            for i in mzml.MzML(os.path.join(mzmls_dir, m)):
                if i['ms level'] == 1:
                    ms1_number += 1
                elif i['ms level'] == 2:
                    ms2_number += 1
                if i['id'] in self.identified_spectrum[m]:
                    if i['base peak intensity'] > 10000:
                        self.peak_intensity_distribution['identified_spectra']['>10000'] += 1
                    elif i['base peak intensity'] > 6000:
                        self.peak_intensity_distribution['identified_spectra']['6000-10000'] += 1
                    elif i['base peak intensity'] > 3000:
                        self.peak_intensity_distribution['identified_spectra']['3000-6000'] += 1
                    elif i['base peak intensity'] > 1000:
                        self.peak_intensity_distribution['identified_spectra']['1000-3000'] += 1
                    elif i['base peak intensity'] > 900:
                        self.peak_intensity_distribution['identified_spectra']['900-1000'] += 1
                    elif i['base peak intensity'] > 700:
                        self.peak_intensity_distribution['identified_spectra']['700-900'] += 1
                    elif i['base peak intensity'] > 500:
                        self.peak_intensity_distribution['identified_spectra']['500-700'] += 1
                    elif i['base peak intensity'] > 300:
                        self.peak_intensity_distribution['identified_spectra']['300-500'] += 1
                    elif i['base peak intensity'] > 100:
                        self.peak_intensity_distribution['identified_spectra']['100-300'] += 1
                    elif i['base peak intensity'] > 10:
                        self.peak_intensity_distribution['identified_spectra']['10-100'] += 1
                    else:
                        self.peak_intensity_distribution['identified_spectra']['0-10'] += 1

                    if 'precursorList' in i.keys():
                        try:
                            charge_state = i['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0][
                                'charge state']
                            if charge_state > 7:
                                self.charge_state_distribution['identified_spectra']['>7'] += 1
                            elif charge_state == 7:
                                self.charge_state_distribution['identified_spectra']['7'] += 1
                            elif charge_state == 6:
                                self.charge_state_distribution['identified_spectra']['6'] += 1
                            elif charge_state == 5:
                                self.charge_state_distribution['identified_spectra']['5'] += 1
                            elif charge_state == 4:
                                self.charge_state_distribution['identified_spectra']['4'] += 1
                            elif charge_state == 3:
                                self.charge_state_distribution['identified_spectra']['3'] += 1
                            elif charge_state == 2:
                                self.charge_state_distribution['identified_spectra']['2'] += 1
                            elif charge_state == 1:
                                self.charge_state_distribution['identified_spectra']['1'] += 1
                        except KeyError:
                            # log.warning("No charge state: {}".format(i['id']))
                            pass

                    if i['ms level'] == 2:
                        if i['defaultArrayLength'] >= 1000:
                            self.peak_per_ms2['identified_spectra']['>1000'] += 1
                        elif 900 <= i['defaultArrayLength'] < 1000:
                            self.peak_per_ms2['identified_spectra']['900-1000'] += 1
                        elif 800 <= i['defaultArrayLength'] < 900:
                            self.peak_per_ms2['identified_spectra']['800-900'] += 1
                        elif 700 <= i['defaultArrayLength'] < 800:
                            self.peak_per_ms2['identified_spectra']['700-800'] += 1
                        elif 600 <= i['defaultArrayLength'] < 700:
                            self.peak_per_ms2['identified_spectra']['600-700'] += 1
                        elif 500 <= i['defaultArrayLength'] < 600:
                            self.peak_per_ms2['identified_spectra']['500-600'] += 1
                        elif 400 <= i['defaultArrayLength'] < 500:
                            self.peak_per_ms2['identified_spectra']['400-500'] += 1
                        elif 300 <= i['defaultArrayLength'] < 400:
                            self.peak_per_ms2['identified_spectra']['300-400'] += 1
                        elif 200 <= i['defaultArrayLength'] < 300:
                            self.peak_per_ms2['identified_spectra']['200-300'] += 1
                        elif 100 <= i['defaultArrayLength'] < 200:
                            self.peak_per_ms2['identified_spectra']['100-200'] += 1
                        else:
                            self.peak_per_ms2['identified_spectra']['0-100'] += 1
                else:
                    if i['base peak intensity'] > 10000:
                        self.peak_intensity_distribution['unidentified_spectra']['>10000'] += 1
                    elif i['base peak intensity'] > 6000:
                        self.peak_intensity_distribution['unidentified_spectra']['6000-10000'] += 1
                    elif i['base peak intensity'] > 3000:
                        self.peak_intensity_distribution['unidentified_spectra']['3000-6000'] += 1
                    elif i['base peak intensity'] > 1000:
                        self.peak_intensity_distribution['unidentified_spectra']['1000-3000'] += 1
                    elif i['base peak intensity'] > 900:
                        self.peak_intensity_distribution['unidentified_spectra']['900-1000'] += 1
                    elif i['base peak intensity'] > 700:
                        self.peak_intensity_distribution['unidentified_spectra']['700-900'] += 1
                    elif i['base peak intensity'] > 500:
                        self.peak_intensity_distribution['unidentified_spectra']['500-700'] += 1
                    elif i['base peak intensity'] > 300:
                        self.peak_intensity_distribution['unidentified_spectra']['300-500'] += 1
                    elif i['base peak intensity'] > 100:
                        self.peak_intensity_distribution['unidentified_spectra']['100-300'] += 1
                    elif i['base peak intensity'] > 10:
                        self.peak_intensity_distribution['unidentified_spectra']['10-100'] += 1
                    else:
                        self.peak_intensity_distribution['unidentified_spectra']['0-10'] += 1

                    if 'precursorList' in i.keys():
                        try:
                            charge_state = i['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0][
                                'charge state']
                            if charge_state > 7:
                                self.charge_state_distribution['unidentified_spectra']['>7'] += 1
                            elif charge_state == 7:
                                self.charge_state_distribution['unidentified_spectra']['7'] += 1
                            elif charge_state == 6:
                                self.charge_state_distribution['unidentified_spectra']['6'] += 1
                            elif charge_state == 5:
                                self.charge_state_distribution['unidentified_spectra']['5'] += 1
                            elif charge_state == 4:
                                self.charge_state_distribution['unidentified_spectra']['4'] += 1
                            elif charge_state == 3:
                                self.charge_state_distribution['unidentified_spectra']['3'] += 1
                            elif charge_state == 2:
                                self.charge_state_distribution['unidentified_spectra']['2'] += 1
                            elif charge_state == 1:
                                self.charge_state_distribution['unidentified_spectra']['1'] += 1
                        except KeyError:
                            # log.warning("No charge state: {}".format(i['id']))
                            pass

                    if i['ms level'] == 2:
                        if i['defaultArrayLength'] >= 1000:
                            self.peak_per_ms2['unidentified_spectra']['>1000'] += 1
                        elif 900 <= i['defaultArrayLength'] < 1000:
                            self.peak_per_ms2['unidentified_spectra']['900-1000'] += 1
                        elif 800 <= i['defaultArrayLength'] < 900:
                            self.peak_per_ms2['unidentified_spectra']['800-900'] += 1
                        elif 700 <= i['defaultArrayLength'] < 800:
                            self.peak_per_ms2['unidentified_spectra']['700-800'] += 1
                        elif 600 <= i['defaultArrayLength'] < 700:
                            self.peak_per_ms2['unidentified_spectra']['600-700'] += 1
                        elif 500 <= i['defaultArrayLength'] < 600:
                            self.peak_per_ms2['unidentified_spectra']['500-600'] += 1
                        elif 400 <= i['defaultArrayLength'] < 500:
                            self.peak_per_ms2['unidentified_spectra']['400-500'] += 1
                        elif 300 <= i['defaultArrayLength'] < 400:
                            self.peak_per_ms2['unidentified_spectra']['300-400'] += 1
                        elif 200 <= i['defaultArrayLength'] < 300:
                            self.peak_per_ms2['unidentified_spectra']['200-300'] += 1
                        elif 100 <= i['defaultArrayLength'] < 200:
                            self.peak_per_ms2['unidentified_spectra']['100-200'] += 1
                        else:
                            self.peak_per_ms2['unidentified_spectra']['0-100'] += 1
                if 'precursorList' in i.keys():
                    total_charge_num += 1
                    try:
                        charge_state = i['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0][
                            'charge state']
                        if charge_state == 2:
                            charge_2 += 1
                    except Exception as e:
                        pass
            heatmap_charge[m] = charge_2 / total_charge_num
            self.Total_ms2_Spectral = self.Total_ms2_Spectral + ms2_number
            mzml_table[m] = {'MS1_Num': ms1_number}
            mzml_table[m]['MS2_Num'] = ms2_number
        median = np.median(list(heatmap_charge.values()))
        self.heatmap_charge_score = dict(zip(heatmap_charge.keys(),
                                             list(map(lambda v: 1 - np.abs(v - median),
                                                      heatmap_charge.values()))))

        raw_ids = config.kwargs['raw_ids']
        for raw_id in os.listdir(raw_ids):
            log.warning("Parsing {}...".format(raw_id))
            if 'msgf' in raw_id:
                mz = openms.idxml.IDXML(os.path.join(raw_ids, raw_id))

                # spectrum matched to target peptide
                if config.kwargs['remove_decoy']:
                    msgf_identified_num = len(set([i['spectrum_reference'] for i in mz
                                                   if i['PeptideHit'][0]['target_decoy'] == 'target']))
                else:
                    msgf_identified_num = len(set([i['spectrum_reference'] for i in mz]))

                mzML_name = raw_id.split('_msgf')[0] + '.mzML'
                mzml_table[mzML_name]['MSGF'] = msgf_identified_num

            if 'comet' in raw_id:
                mz = openms.idxml.IDXML(os.path.join(raw_ids, raw_id))
                if config.kwargs['remove_decoy']:
                    comet_identified_num = len(set([i['spectrum_reference'] for i in mz
                                                    if i['PeptideHit'][0]['target_decoy'] == 'target']))
                else:
                    comet_identified_num = len(set([i['spectrum_reference'] for i in mz]))
                mzML_name = raw_id.split('_comet')[0] + '.mzML'
                mzml_table[mzML_name]['Comet'] = comet_identified_num

            mzml_table[mzML_name]['Final result of spectra'] = self.mL_spec_ident_final[mzML_name]
            mzml_table[mzML_name]['Final result of Peptides'] = len(self.mzml_peptide_map[mzML_name])

        self.mzml_table = mzml_table

    def parse_out_tmt(self):
        log.warning("Parsing mzTab file...")
        exp_data = pd.read_csv(self.exp_design, sep='\t', header=0, index_col=None, dtype=str)
        exp_data = exp_data.dropna(axis=0)
        Spec_File = exp_data['Spectra_Filepath'].tolist()
        mztab_data = mztab.MzTab(self.out_mzTab_path)
        pep_table = mztab_data.peptide_table
        meta_data = dict(mztab_data.metadata)
        self.delta_mass['target'] = dict()
        self.delta_mass['decoy'] = dict()

        # PSM table data
        psm = mztab_data.spectrum_match_table
        psm_table = dict()
        psm['stand_spectra_ref'] = psm.apply(
            lambda x: os.path.basename(meta_data[x.spectra_ref.split(':')[0] + '-location']), axis=1)
        pep_table['stand_spectra_ref'] = pep_table.apply(
            lambda x: os.path.basename(meta_data[x.spectra_ref.split(':')[0] + '-location']), axis=1)

        if config.kwargs['remove_decoy']:
            Total_Protein_Identified = set(psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] == 0]
                                           ['accession'])
            self.Total_Protein_Identified = len(Total_Protein_Identified)
        else:
            Total_Protein_Identified = set(psm['accession'])
            self.Total_Protein_Identified = len(Total_Protein_Identified)

        num_pep_per_protein = dict()
        percent_pep_per_protein = dict()
        for protein in Total_Protein_Identified:
            number = str(len(set(psm[psm['accession'] == protein]['sequence'])))
            if number in num_pep_per_protein:
                num_pep_per_protein[number]['Frequency'] += 1
                percent_pep_per_protein[number]['Percentage'] = num_pep_per_protein[number]['Frequency'] / \
                                                                self.Total_Protein_Identified
            else:
                num_pep_per_protein[number] = {'Frequency': 1}
                percent_pep_per_protein[number] = {'Percentage': 1 / self.Total_Protein_Identified}

        keys = sorted(num_pep_per_protein.items(), key=lambda x: int(x[0]))  # sort
        for key in keys:
            self.num_pep_per_protein[key[0]] = key[1]
            self.percent_pep_per_protein[key[0]] = percent_pep_per_protein[key[0]]

        # PSM information
        scores = list(filter(lambda x: re.match(r'search_engine_score.*?', x) is not None,
                             psm.columns.tolist()))
        t = (' float'.join(scores) + ' float').replace("[", "_").replace("]", "")
        cur.execute("CREATE TABLE PSM(PSM_ID int, sequence VARCHAR(100), PSM_UNIQUE int, " + t + ")")
        con.commit()

        mL_spec_ident_final = {}
        for m in set(psm['stand_spectra_ref'].tolist()):
            Sample = list(exp_data[exp_data['Spectra_Filepath'] == m]['Sample'])
            self.cal_num_table_data[m] = {'Sample Name': '|'.join(Sample)}
            fraction = exp_data[exp_data['Spectra_Filepath'] == m]['Fraction'].tolist()[0]
            self.cal_num_table_data[m]['fraction'] = fraction

            if config.kwargs['remove_decoy']:
                proteins = set(psm[(psm['stand_spectra_ref'] == m) &
                                   (psm['opt_global_cv_MS:1002217_decoy_peptide'] == 0)]['accession'])
                peptides = set(psm[(psm['stand_spectra_ref'] == m) &
                                   (psm['opt_global_cv_MS:1002217_decoy_peptide'] == 0)]
                               ['opt_global_cv_MS:1000889_peptidoform_sequence'])

                self.identified_spectrum[m] = list(map(lambda x: x.split(':')[1],
                                                       psm[(psm['stand_spectra_ref'] == m) & (
                                                               psm['opt_global_cv_MS:1002217_decoy_peptide'] != 1)][
                                                           'spectra_ref'].tolist()))
                self.mzml_peptide_map[m] = list(set(psm[(psm['stand_spectra_ref'] == m) &
                                                        (psm['opt_global_cv_MS:1002217_decoy_peptide'] != 1)][
                                                        'sequence'].tolist()))
            else:
                proteins = set(psm[psm['stand_spectra_ref'] == m]['accession'])
                peptides = set(psm[psm['stand_spectra_ref'] == m]
                               ['opt_global_cv_MS:1000889_peptidoform_sequence'])

                self.identified_spectrum[m] = list(map(lambda x: x.split(':')[1],
                                                       psm[psm['stand_spectra_ref'] == m]['spectra_ref'].tolist()))
                self.mzml_peptide_map[m] = list(set(psm[psm['stand_spectra_ref'] == m]['sequence'].tolist()))

            if None in proteins:
                proteins.remove(None)
            self.cal_num_table_data[m]['protein_num'] = len(set(proteins))
            self.cal_num_table_data[m]['peptide_num'] = len(peptides)

            modified_pep = list(filter(lambda x: re.match(r'.*?\(.*\).*?', x) is not None, peptides))
            self.cal_num_table_data[m]['modified_peptide_num'] = len(modified_pep)

            mL_spec_ident_final[m] = len(self.identified_spectrum[m])

        psm['relative_diff'] = psm['exp_mass_to_charge'] - psm['calc_mass_to_charge']
        self.delta_mass['decoy'] = Counter(psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] == 1]
                                           ['relative_diff'].tolist())
        self.delta_mass['target'] = Counter(psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] != 1]
                                            ['relative_diff'].tolist())

        sql_t = "(" + ','.join(['?'] * (len(scores) + 3)) + ")"
        if config.kwargs['remove_decoy']:
            cur.executemany("INSERT INTO PSM (PSM_ID,sequence,PSM_UNIQUE," +
                            ','.join(scores).replace("[", "_").replace("]", "") + ") VALUES " + sql_t,
                            [tuple(x) for x in psm[psm['opt_global_cv_MS:1002217_decoy_peptide'] != 1]
                            [['PSM_ID', 'sequence', 'unique'] + scores].values])
        else:
            cur.executemany("INSERT INTO PSM (PSM_ID,sequence,PSM_UNIQUE," +
                            ','.join(scores).replace("[", "_").replace("]", "") + ") VALUES " + sql_t,
                            [tuple(x) for x in psm[['PSM_ID', 'sequence', 'unique'] + scores].values])
        con.commit()
        for idx, row in psm.iterrows():
            if config.kwargs['remove_decoy'] and row['opt_global_cv_MS:1002217_decoy_peptide'] == 1:
                continue
            else:
                psm_table[row['PSM_ID']] = {'sequence': row['sequence']}
                # psm_table[row['PSM_ID']]['spectra_ref'] = row['spectra_ref']
                psm_table[row['PSM_ID']]['unique'] = row['unique']
                for score in scores:
                    psm_table[row['PSM_ID']][score] = row[score]
            if idx > 48:
                break

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
        self.PSM_table = psm_table

        # peptide quantification table data
        pep_quant = dict()
        study_variables = list(filter(lambda x: re.match(r'peptide_abundance_study_variable.*?', x) is not None,
                                      pep_table.columns.tolist()))
        cur.execute("CREATE TABLE QUANT(ID integer PRIMARY KEY AUTOINCREMENT, sequence VARCHAR(100))")
        con.commit()
        sql_col = "sequence"
        sql_t = "(" + ','.join(['?'] * (len(study_variables) + 1)) + ")"

        for s in study_variables:
            s = s.replace("[", "_").replace("]", "")
            cur.execute("ALTER TABLE QUANT ADD " + s + " FLOAT")
            con.commit()
            sql_col += "," + s
        cur.executemany("INSERT INTO QUANT (" + sql_col + ") VALUES " + sql_t,
                        [tuple(x) for x in pep_table[['sequence'] + study_variables].values])
        con.commit()
        for index, row in pep_table.iterrows():
            if config.kwargs['remove_decoy'] and row['opt_global_cv_MS:1002217_decoy_peptide'] == 1:
                continue
            else:
                pep_quant[str(index + 1)] = {'sequence': row['sequence']}  # keep id unique
                for s in study_variables:
                    pep_quant[str(index + 1)][s] = row[s]
            if index > 48:
                break

        self.pep_quant_table = pep_quant
