#!/usr/bin/env python

""" MultiQC example plugin module """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table
from multiqc.modules.base_module import BaseMultiqcModule
import pandas as pd
import re
from pyteomics import mztab, mzml, openms
import os

# Initialise the main MultiQC logger
log = logging.getLogger('__name__')


class MultiqcModule(BaseMultiqcModule):
    css = None

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

        self.num_pep_per_protein = dict()
        self.out_csv_data = dict()
        self.cal_num_table_data = dict()

        # draw_proteomicslfq_identi_num
        self.draw_proteomicslfq_identi_num()

        # draw number of peptides per protein
        self.draw_num_pep_per_protein()

        # draw peptides quantification information

        self.draw_pep_quant_info()

        self.draw_psm_table()
        self.draw_mzml_ms()

        # exist some problems
        # self.css = {
        #     './assets/css/proteomicslfq.css':
        #         os.path.join(os.path.dirname(__file__), 'assets', 'css', 'proteomicslfq.css')
        # }

    def draw_proteomicslfq_identi_num(self):
        self.parse_input_data_v2(self.out_csv_path, config.kwargs['exp_design'])
        self.write_data_file(self.cal_num_table_data, 'result_statistics_table')
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
            description='This plot shows the ProteomicsLFQ pipeline final result',
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

    # draw number of peptides per proteins
    def draw_num_pep_per_protein(self):
        data = pd.read_csv(self.out_csv_path, sep=',', header=0)
        data = data.astype(str)
        PeptideSequences = data['PeptideSequence'].tolist()
        pep_count = {}
        for i in set(PeptideSequences):
            pep_count[i] = PeptideSequences.count(i)

        for pep, count in pep_count.items():
            ProteinNames = data[data['PeptideSequence'] == pep]['ProteinName'].tolist()
            for protein in ProteinNames:
                if len(protein.split(';')) == 1:
                    value = {'number of peptides per proteins': count}
                    self.num_pep_per_protein[protein] = value
                else:
                    for p in protein.split(';'):
                        value = {'number of peptides per proteins': count}
                        self.num_pep_per_protein[p] = value
        # Create table plot
        pconfig = {
            'id': 'number_of_peptides_per_proteins',  # ID used for the table
            'table_title': 'number_of_peptides_per_proteins',  # Title of the table. Used in the column config modal
            'save_file': False,  # Whether to save the table data to a file
            'raw_data_fn': 'multiqc_number_of_peptides_per_proteins_table',  # File basename to use for raw data file
            'sortRows': True,  # Whether to sort rows alphabetically
            'only_defined_headers': False,  # Only show columns that are defined in the headers config
            'col1_header': 'ProteinName',
            'format': '{:,.0f}',  # The header used for the first column
            'scale': 'BuGn'
        }
        headers = OrderedDict()
        headers['number of peptides per proteins'] = {
            'description': 'number of peptides per proteins',
            'color': "#ffffff"
        }
        # write to file
        self.write_data_file(self.num_pep_per_protein, 'num_pep_per_protein')

        if len(self.num_pep_per_protein) > 500:
            L = 0
            new_num_pep_per_protein = dict()
            for key, value in self.num_pep_per_protein.items():
                new_num_pep_per_protein[key] = value
                L += 1
                if L == 500:
                    break
        else:
            new_num_pep_per_protein = self.num_pep_per_protein

        table_html = table.plot(new_num_pep_per_protein, headers, pconfig)
        # Add a report section with the line plot
        self.add_section(
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
            plot=table_html
        )

    def draw_pep_quant_info(self):
        pep_table = mztab.MzTab(self.out_mzTab_path).peptide_table
        meta_data = mztab.MzTab(self.out_mzTab_path).metadata
        pep_quant = dict()
        quant_method = dict(meta_data)['quantification_method']
        study_variables = list(filter(lambda x: re.match(r'peptide_abundance_study_variable.*?', x) is not None,
                                      pep_table.columns.tolist()))
        for index, row in pep_table.iterrows():
            pep_quant[str(index)] = {'peptide_sequence': row['sequence']}  # keep id unique
            pep_quant[str(index)]['accession'] = row['accession']
            pep_quant[str(index)]['spectra_ref'] = row['spectra_ref']
            # pep_quant[str(index)]['quantification_method'] = quant_method
            for s in study_variables:
                pep_quant[str(index)][s] = row[s]
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
        self.write_data_file(pep_quant, 'pep_quant')

        if len(pep_quant) > 500:
            L = 0
            new_pep_quant = dict()
            for key, value in pep_quant.items():
                new_pep_quant[key] = value
                L += 1
                if L == 500:
                    break
        else:
            new_pep_quant = pep_quant
        table_html = table.plot(new_pep_quant, headers, pconfig)
        # Add a report section with the line plot
        self.add_section(
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
        psm = mztab.MzTab(self.out_mzTab_path).spectrum_match_table
        psm_table = dict()
        scores = list(filter(lambda x: re.match(r'search_engine_score.*?', x) is not None,
                             psm.columns.tolist()))
        for _, row in psm.iterrows():
            psm_table[row['PSM_ID']] = {'sequence': row['sequence']}
            psm_table[row['PSM_ID']]['spectra_ref'] = row['spectra_ref']
            psm_table[row['PSM_ID']]['unique'] = row['unique']
            for score in scores:
                psm_table[row['PSM_ID']][score] = row[score]

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
        self.write_data_file(psm_table, 'psm_table')

        if len(psm_table) > 500:
            L = 0
            new_psm_table = {}
            for key, value in psm_table.items():
                new_psm_table[key] = value
                L += 1
                if L == 500:
                    break
        else:
            new_psm_table = psm_table
        table_html = table.plot(new_psm_table, headers, pconfig)
        # Add a report section with the line plot
        self.add_section(
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

    def parse_input_data_v2(self, out_csv, exp_design):
        data = pd.read_csv(out_csv, sep=',', header=0)
        data = data.astype(str)
        exp_data = pd.read_csv(exp_design, sep='\t', header=0, index_col=None, dtype=str)
        exp_data = exp_data.dropna(axis=0)
        Spec_File = exp_data['Spectra_Filepath'].tolist()
        for i in Spec_File:
            Sample = list(exp_data[exp_data['Spectra_Filepath'] == i]['Sample'])[0]
            self.cal_num_table_data[i] = {'Sample Name': Sample}
            condition = '-'.join(list(set(data[data['Reference'] == i]['Condition'].tolist())))
            self.cal_num_table_data[i]['condition'] = condition
            fraction = exp_data[exp_data['Spectra_Filepath'] == i]['Fraction'].tolist()[0]
            self.cal_num_table_data[i]['fraction'] = fraction
            peptides = list(set(data[data['Reference'] == i]['PeptideSequence'].tolist()))
            self.cal_num_table_data[i]['peptide_num'] = len(peptides)
            modified_pep = list(filter(lambda x: re.match(r'.*?\(.*\).*?', x) is not None, peptides))
            self.cal_num_table_data[i]['modified_peptide_num'] = len(modified_pep)
            protein_num = len(list(set(data[data['Reference'] == i]['ProteinName'].tolist())))
            self.cal_num_table_data[i]['protein_num'] = protein_num

    def draw_mzml_ms(self):
        ms1_number = 0
        ms2_number = 0
        mzml_table = {}
        mzmls_dir = config.kwargs['mzMLs']

        for m in os.listdir(mzmls_dir):
            print(os.path.join(mzmls_dir, m))
            for i in mzml.MzML(os.path.join(mzmls_dir, m)):
                if i['ms level'] == 1:
                    ms1_number += 1
                elif i['ms level'] == 2:
                    ms2_number += 1
            mzml_table[m] = {'MS1_Num': ms1_number}
            mzml_table[m]['MS2_Num'] = ms2_number

        raw_ids = config.kwargs['raw_ids']
        for raw_id in os.listdir(raw_ids):
            if 'msgf' in raw_id:
                mz = openms.idxml.IDXML(os.path.join(raw_ids, raw_id))
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

        # write to file
        self.write_data_file(mzml_table, 'mzMLs_statistics_table')

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
        table_html = table.plot(mzml_table, headers, pconfig)

        # Add a report section with the line plot
        self.add_section(
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
