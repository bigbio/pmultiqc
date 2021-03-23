#!/usr/bin/env python

""" MultiQC example plugin module """

from __future__ import absolute_import
from collections import OrderedDict
import logging
import json
from multiqc import config
from multiqc.plots import table, bargraph
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.utils import report
import pandas as pd
import re
from pyteomics import mztab, mzml, openms
import os

# Initialise the main MultiQC logger
log = logging.getLogger('__name__')


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
            info=" is an example module to show how the MultQC pluginm system works."
        )

        for file in self.find_log_files({'fn': "out_msstats.csv"}):
            self.out_csv_path = os.path.join(config.analysis_dir[0], file['fn'])

        for f in self.find_log_files({'fn': 'out.mzTab'}):
            self.out_mzTab_path = os.path.join(config.analysis_dir[0], f['fn'])

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
        # parse input data
        self.parse_out_csv()
        self.parse_out_mzTab()
        self.parse_mzml_idx()

        # draw the experimental design
        self.draw_exp_design()
        #
        self.draw_summary_proten_ident_table()
        #
        # # draw_proteomicslfq_identi_num
        self.draw_proteomicslfq_identi_num()

        # draw number of peptides per protein
        self.draw_num_pep_per_protein()

        # draw peptides quantification information

        self.draw_pep_quant_info()
        #
        self.draw_psm_table()
        self.draw_mzml_ms()
        self.draw_peak_intensity_distribution()
        # exist some problems
        # self.css = {
        #     './assets/css/proteomicslfq.css':
        #         os.path.join(os.path.dirname(__file__), 'assets', 'css', 'proteomicslfq.css')
        # }

    def draw_exp_design(self):
        exp_design_table = dict()
        with open(config.kwargs['exp_design'], 'r') as f:
            data = f.readlines()
            empty_row = data.index('\n')
            f_table = data[1:empty_row]
            s_table = [i.replace('\n', '').split('\t') for i in data[data.index('\n') + 1:]][1:]
            s_header = data[data.index('\n') + 1].replace('\n', '').split('\t')
            s_DataFrame = pd.DataFrame(s_table, columns=s_header)
            for index in range(len(f_table)):
                # TODO be improved
                exp_design_table[f_table[index].split('\t')[2]] = {'Fraction_Group': f_table[index].split('\t')[0]}
                exp_design_table[f_table[index].split('\t')[2]]['Fraction'] = f_table[index].split('\t')[1]
                exp_design_table[f_table[index].split('\t')[2]]['Label'] = f_table[index].split('\t')[3]
                sample = f_table[index].split('\t')[4].replace('\n', '')
                exp_design_table[f_table[index].split('\t')[2]]['Sample'] = sample
                exp_design_table[f_table[index].split('\t')[2]]['MSstats_Condition'] = \
                    s_DataFrame[s_DataFrame['Sample'] == sample]['MSstats_Condition'].tolist()[0]
                exp_design_table[f_table[index].split('\t')[2]]['MSstats_BioReplicate'] = \
                    s_DataFrame[s_DataFrame['Sample'] == sample]['MSstats_BioReplicate'].tolist()[0]

        # Create table plot
        pconfig = {
            'id': 'Experimental_Design',  # ID used for the table
            'table_title': 'Experimental Design',  # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'raw_data_fn': 'multiqc_Experimental_Design_table',  # File basename to use for raw data file
            'sortRows': False,  # Whether to sort rows alphabetically
            'only_defined_headers': False,  # Only show columns that are defined in the headers config
            'col1_header': 'Spectra File',
            'format': '{:,.0f}'  # The header used for the first column
        }
        headers = OrderedDict()

        headers['Fraction_Group'] = {
            'description': 'Fraction_Group',
            'color': "#ffffff"
        }
        headers['Fraction'] = {
            'description': 'Fraction identifier',
            'color': "#ffffff"
        }
        headers['Label'] = {
            'description': 'Label',
            'color': "#ffffff"
        }
        headers['Sample'] = {
            'description': 'Sample Name',
            'color': "#ffffff"
        }
        headers['MSstats_Condition'] = {
            'description': 'MSstats Condition',
            'color': "#ffffff"
        }
        headers['MSstats_BioReplicate'] = {
            'description': 'MSstats BioReplicate',
            'color': "#ffffff"
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
            and identified the number of modified peptide in the pipeline, eg.
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

        # write to file
        self.write_data_file(self.num_pep_per_protein, 'num_pep_per_protein')

        bar_html = bargraph.plot([self.num_pep_per_protein, self.percent_pep_per_protein], headers, pconfig)
        # Add a report section with the line plot
        self.add_section(
            name="Number Of Peptides Per Proteins",
            anchor="num_of_pep_per_prot",
            description='This plot shows the number of peptides per proteins '
                        'in ProteomicsLFQ pipeline final result',
            helptext='''
                    This longer description explains what exactly the numbers mean
                    and supports markdown formatting. This means that we can do _this_:

                    * Something important
                    * Something else important
                    * Best of all - some `code`

                    Doesn't matter if this is copied from documentation - makes it
                    easier for people to find quickly.
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
            'col1_header': 'Index',
            'format': '{:,.6f}',  # The header used for the first column
            'scale': 'Set3'
        }

        headers = OrderedDict()
        headers['accession'] = {
            'description': 'accession number of sequence',
            'color': "#ffffff",
            'format': '{:,.0f}',
            'hidden': True
        }

        headers['peptide_sequence'] = {
            'description': 'peptide_sequence',
            'color': "#ffffff",
            'format': '{:,.0f}'
        }

        headers['spectra_ref'] = {
            'description': 'spectra_reference',
            'color': "#ffffff",
            'format': '{:,.0f}'
        }

        # write to file
        self.write_data_file(self.pep_quant_table, 'pep_quant')

        if len(self.pep_quant_table) > 500:
            L = 0
            new_pep_quant = dict()
            for key, value in self.pep_quant_table.items():
                new_pep_quant[key] = value
                L += 1
                if L == 500:
                    break
        else:
            new_pep_quant = self.pep_quant_table
        table_html = table.plot(new_pep_quant, headers, pconfig)

        self.add_section(
            name="Quantification Result",
            anchor="quant_result",
            description='This plot shows the quantification information of peptides'
                        'in ProteomicsLFQ pipeline final result',
            helptext='''
                            This longer description explains what exactly the numbers mean
                            and supports markdown formatting. This means that we can do _this_:

                            * Something important
                            * Something else important
                            * Best of all - some `code`

                            Doesn't matter if this is copied from documentation - makes it
                            easier for people to find quickly.
                            ''',
            plot=table_html
        )

    def draw_psm_table(self):
        pconfig = {
            'id': 'peptide spectrum match',  # ID used for the table
            'table_title': 'peptide spectrum match information',
            # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'raw_data_fn': 'multiqc_psm_table',  # File basename to use for raw data file
            'sortRows': False,  # Whether to sort rows alphabetically
            'only_defined_headers': False,  # Only show columns that are defined in the headers config
            'col1_header': 'PSM_ID',
            'format': '{:,.9f}',  # The header used for the first column
            'scale': 'Set5'
        }

        headers = OrderedDict()
        headers['sequence'] = {
            'description': 'peptide_sequence',
            'color': "#ffffff",
            'format': '{:,.0f}'
        }

        headers['spectra_ref'] = {
            'description': 'spectra_reference',
            'color': "#ffffff",
            'format': '{:,.0f}'
        }
        headers['unique'] = {
            'description': 'unique',
            'color': "#ffffff",
            'format': '{:,.0f}'
        }

        # write to file
        self.write_data_file(self.PSM_table, 'psm_table')

        if len(self.PSM_table) > 500:
            L = 0
            new_psm_table = {}
            for key, value in self.PSM_table.items():
                new_psm_table[key] = value
                L += 1
                if L == 500:
                    break
        else:
            new_psm_table = self.PSM_table
        table_html = table.plot(new_psm_table, headers, pconfig)
        # Add a report section with the line plot
        self.add_section(
            name="Peptide-Spectrum Matches",
            anchor="psm",
            description='This plot shows the PSM information'
                        'in ProteomicsLFQ pipeline final result',
            helptext='''
                            This longer description explains what exactly the numbers mean
                            and supports markdown formatting. This means that we can do _this_:

                            * Something important
                            * Something else important
                            * Best of all - some `code`

                            Doesn't matter if this is copied from documentation - makes it
                            easier for people to find quickly.
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
        headers['MSGF'] = {
            'description': 'Number of spectra identified by MSGF search engine',
            'color': "#ffffff"
        }
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
                    This longer description explains what exactly the numbers mean
                    and supports markdown formatting. This means that we can do _this_:

                    * Something important
                    * Something else important
                    * Best of all - some `code`

                    Doesn't matter if this is copied from documentation - makes it
                    easier for people to find quickly.
                    ''',
            plot=table_html
        )

    def draw_peak_intensity_distribution(self):
        # Create table plot
        pconfig = {
            'id': 'peak_intensity_distribution',  # ID used for the table
            'cpswitch': False,
            'title': 'Peak Intensity Distribution',
            'xlab': 'Peak Intensity'
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

        # write to file
        self.write_data_file(self.peak_intensity_distribution, 'peak_intensity_distribution')
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

    def parse_out_csv(self):

        # result statistics table
        data = pd.read_csv(self.out_csv_path, sep=',', header=0)
        data = data.astype(str)
        exp_data = pd.read_csv(config.kwargs['exp_design'], sep='\t', header=0, index_col=None, dtype=str)
        exp_data = exp_data.dropna(axis=0)
        Spec_File = exp_data['Spectra_Filepath'].tolist()
        for i in Spec_File:
            proteins = []
            tag = True
            index = 0
            Sample = list(exp_data[exp_data['Spectra_Filepath'] == i]['Sample'])[0]
            self.cal_num_table_data[i] = {'Sample Name': Sample}
            condition = '-'.join(list(set(data[data['Reference'] == i]['Condition'].tolist())))
            self.cal_num_table_data[i]['condition'] = condition
            fraction = exp_data[exp_data['Spectra_Filepath'] == i]['Fraction'].tolist()[0]
            self.cal_num_table_data[i]['fraction'] = fraction
            ProteinNames = data[data['Reference'] == i]['ProteinName'].tolist()
            if config.kwargs['remove_decoy']:
                for p in ProteinNames:
                    if len(p.split(';')) == 1:
                        if p.startswith('DECOY_'):
                            continue
                        elif p not in proteins:
                            proteins.append(p)
                    else:
                        for j in p.split(';'):
                            if j.startswith('DECOY_'):
                                continue
                            elif j not in proteins:
                                proteins.append(j)
            else:
                for p in ProteinNames:
                    if len(p.split(';')) == 1:
                        if p not in proteins:
                            proteins.append(p)
                    else:
                        for j in p.split(';'):
                            if j not in proteins:
                                proteins.append(j)
            self.cal_num_table_data[i]['protein_num'] = len(proteins)

            peptides = list(set(data[data['Reference'] == i]['PeptideSequence'].tolist()))
            if config.kwargs['remove_decoy']:
                while index < len(peptides):
                    ProteinNames = data[data[data['Reference'] == i]['PeptideSequence'] == peptides[i]]['ProteinName'].tolist()
                    for p in ProteinNames:
                        if len(p.split(';')) == 1:
                            if p.startswith('DECOY_'):
                                continue
                            else:
                                tag = False
                                break
                        else:
                            for j in p.split(';'):
                                if j.startswith('DECOY_'):
                                    continue
                                else:
                                    tag = False
                                    break
                            if not tag:
                                break
                    if tag:
                        peptides.remove(peptides[i])
                    else:
                        index += 1

            else:
                self.cal_num_table_data[i]['peptide_num'] = len(peptides)

            modified_pep = list(filter(lambda x: re.match(r'.*?\(.*\).*?', x) is not None, peptides))
            self.cal_num_table_data[i]['modified_peptide_num'] = len(modified_pep)

        self.write_data_file(self.cal_num_table_data, 'result_statistics_table')

        #  table of peptide number per protein
        prot_pep_map = {}
        num_pep_per_protein = {}
        for _, row in data.iterrows():
            ProteinName = row['ProteinName']
            if len(ProteinName.split(';')) == 1:
                if config.kwargs['remove_decoy'] and ProteinName.startswith('DECOY_'):
                    continue
                if ProteinName not in prot_pep_map.keys():
                    prot_pep_map[ProteinName] = [row['PeptideSequence']]
                elif row['PeptideSequence'] not in prot_pep_map[ProteinName]:
                    prot_pep_map[ProteinName].append(row['PeptideSequence'])
            else:
                for p in ProteinName.split(';'):
                    if config.kwargs['remove_decoy'] and p.startswith('DECOY_'):
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
        mztab_data = mztab.MzTab(self.out_mzTab_path)
        pep_table = mztab_data.peptide_table
        meta_data = dict(mztab_data.metadata)

        # PSM table data
        psm = mztab_data.spectrum_match_table
        psm_table = dict()
        scores = list(filter(lambda x: re.match(r'search_engine_score.*?', x) is not None,
                             psm.columns.tolist()))
        spectrum = []
        seq = []
        mL_spec_ident_final = {}
        for _, row in psm.iterrows():
            if config.kwargs['remove_decoy'] and row['opt_global_cv_MS:1002217_decoy_peptide'] == 1:
                continue
            else:
                psm_table[row['PSM_ID']] = {'sequence': row['sequence']}
                psm_table[row['PSM_ID']]['spectra_ref'] = row['spectra_ref']
                psm_table[row['PSM_ID']]['unique'] = row['unique']
                for score in scores:
                    psm_table[row['PSM_ID']][score] = row[score]

                mzML_Name = os.path.basename(meta_data[row['spectra_ref'].split(':')[0] + '-location'])
                if row['spectra_ref'] not in spectrum:
                    spectrum.append(row['spectra_ref'])
                    if mzML_Name not in self.identified_spectrum.keys():
                        self.identified_spectrum[mzML_Name] = [row['spectra_ref'].split(':')[1]]
                    elif row['spectra_ref'].split(':')[1] not in self.identified_spectrum[mzML_Name]:
                        self.identified_spectrum[mzML_Name].append(row['spectra_ref'].split(':')[1])

                    if mzML_Name not in mL_spec_ident_final:
                        mL_spec_ident_final[mzML_Name] = 1
                    else:
                        mL_spec_ident_final[mzML_Name] += 1
                if row['sequence'] not in seq:
                    seq.append(row['sequence'])
                    if mzML_Name not in self.mzml_peptide_map:
                        self.mzml_peptide_map[mzML_Name] = 1
                    else:
                        self.mzml_peptide_map[mzML_Name] += 1

        self.mL_spec_ident_final = mL_spec_ident_final
        self.Total_ms2_Spectral_Identified = len(spectrum)
        self.Total_Peptide_Count = len(seq)
        self.PSM_table = psm_table

        # peptide quantification table data
        pep_quant = dict()
        quant_method = meta_data['quantification_method']
        study_variables = list(filter(lambda x: re.match(r'peptide_abundance_study_variable.*?', x) is not None,
                                      pep_table.columns.tolist()))
        for index, row in pep_table.iterrows():
            if config.kwargs['remove_decoy'] and row['opt_global_cv_MS:1002217_decoy_peptide'] == 1:
                continue
            else:
                pep_quant[str(index)] = {'peptide_sequence': row['sequence']}  # keep id unique
                pep_quant[str(index)]['accession'] = row['accession']
                pep_quant[str(index)]['spectra_ref'] = row['spectra_ref']
                # pep_quant[str(index)]['quantification_method'] = quant_method
                for s in study_variables:
                    pep_quant[str(index)][s] = row[s]

        self.pep_quant_table = pep_quant

    def parse_mzml_idx(self):
        ms1_number = 0
        ms2_number = 0
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

        mzml_table = {}
        mzmls_dir = config.kwargs['mzMLs']

        for m in os.listdir(mzmls_dir):
            # print(os.path.join(mzmls_dir, m))
            log.info("Parsing {}...".format(os.path.join(mzmls_dir, m)))
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
                        charge_state = i['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state']
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
                            charge_state = i['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state']
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
                            log.warning("No charge state: {}".format(i['id']))
                            # print("No charge state: " + i['id'])
                            
            self.Total_ms2_Spectral = self.Total_ms2_Spectral + ms2_number
            mzml_table[m] = {'MS1_Num': ms1_number}
            mzml_table[m]['MS2_Num'] = ms2_number

        raw_ids = config.kwargs['raw_ids']
        for raw_id in os.listdir(raw_ids):
            log.info("Parsing {}...".format(raw_ids))
            if 'msgf' in raw_id:
                mz = openms.idxml.IDXML(os.path.join(raw_ids, raw_id))

                # spectrum matched to target peptide
                msgf_identified_num = len(set([i['spectrum_reference'] for i in mz
                                               if i['PeptideHit'][0]['target_decoy'] == 'target']))
                mzML_name = raw_id.split('_msgf')[0] + '.mzML'
                mzml_table[mzML_name]['MSGF'] = msgf_identified_num

            if 'comet' in raw_id:
                mz = openms.idxml.IDXML(os.path.join(raw_ids, raw_id))
                comet_identified_num = len(set([i['spectrum_reference'] for i in mz
                                                if i['PeptideHit'][0]['target_decoy'] == 'target']))
                mzML_name = raw_id.split('_comet')[0] + '.mzML'
                mzml_table[mzML_name]['Comet'] = comet_identified_num

            mzml_table[mzML_name]['Final result of spectrum'] = self.mL_spec_ident_final[mzML_name]
            mzml_table[mzML_name]['Final result of Peptides'] = self.mzml_peptide_map[mzML_name]
        self.mzml_table = mzml_table
        # write to file
        self.write_data_file(mzml_table, 'mzMLs_statistics_table')
