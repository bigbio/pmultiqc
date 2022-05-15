from multiqc.plots import bargraph
from collections import OrderedDict
from multiqc.modules.base_module import BaseMultiqcModule

class Histogram(BaseMultiqcModule):
    def __init__(self):

        # Initialise the parent module Class object
        super(Histogram, self).__init__()

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
                        'in quantms pipeline final result',
            helptext='''
                        This statistic is extracted from the out_msstats file. Proteins supported by more peptide 
                        identifications can constitute more confident results.
                    ''',
            plot=bar_html
        )

    def draw_peak_intensity_distribution(self):
        # Create table plot
        pconfig = {
            'id': 'peak_intensity_distribution',  # ID used for the table
            'cpswitch': False,
            'title': 'Peak Intensity Distribution',
            # 'xlab': 'Peak Intensity'
        }
        peak_intensity_list = ['0-10', '10-100'] + [str(i) + '-' + str(i + 200) 
                                for i in range(100, 900, 200)] + ['900-1000', '1000-3000', '3000-6000', '6000-10000', '>10000']
        cats = OrderedDict()
        for i in peak_intensity_list:
            cats[i] = {
                    'name': i,
                    'description': 'Peak intensity ' + 
                                    ('> 10000' if i == '>10000' else 'is between ' +  i.split('-')[0] + ' and ' + i.split('-')[1])
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

        charge_list = [str(i) for i in range(1, 8)]
        charge_list.append('>7')
        cats = OrderedDict()
        for i in charge_list:
            cats[i] = {
                    'name': i,
                    'description': 'Precursor charge state ' + ('>7' if i == '>7' else 'is ' + i) 
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
        peaks_list = [str(i) + '-' + str(i + 100) for i in range(0, 1000, 100)]
        peaks_list.append('>1000')
        cats = OrderedDict()
        for i in peaks_list:
            cats[i] = {
                    'name': i,
                    'description': 'Number of Peaks per MS/MS spectrum ' + 
                                    ('> 1000' if i == '>1000' else 'is between ' + i.split('-')[0] + ' and ' + i.split('-')[1])
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

    def draw_oversampling(self):

        # Create bar plot
        pconfig = {
            'id': 'Oversampling_Distribution',
            'cpswitch': False,
            'title': 'MS2 counts per 3D-peak',
            'scale': "set3"
        }
        oversampling_list = [str(i) for i in range(1, 4)]
        oversampling_list[2] = oversampling_list[2] + '+'
        cats = OrderedDict()
        for i in oversampling_list:
            cats[i] = {
            'name': i,
             'description': f'A peak whose peptide ion (same sequence and same charge state) was identified by {i} '
                        'distinct MS2 spectra'
            }

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