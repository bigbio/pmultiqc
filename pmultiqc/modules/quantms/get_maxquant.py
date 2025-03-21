import pandas as pd
import re
import numpy as np
import os

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .histogram import Histogram

def read(file_path, file_type=None, filter_type=None):
    MQ_data = pd.read_csv(file_path, sep='\t', low_memory=False)

    if filter_type == 'Reverse':
        # Reverse
            # When marked with ‘+’, this particular protein group contains no protein, 
            # made up of at least 50% of the peptides of the leading protein, with a 
            # peptide derived from the reversed part of the decoy database. 
            # These should be removed for further data analysis. 
            # The 50% rule is in place to prevent spurious protein hits to erroneously flag the protein group as reverse.
        if 'Reverse' in MQ_data.columns:
            MQ_data = MQ_data[MQ_data['Reverse'] != '+'].reset_index(drop=True)

        # TODO Contaminant: Potential contaminant

    # proteinGroups.txt
    if file_type == 'protgroup':

        if 'Mol. weight [kDa]' in MQ_data.columns:

            intensity_cols = [col for col in MQ_data.columns if re.search(r'intensity', col, re.IGNORECASE)]
            new_intensity_cols = [re.sub(r"intensity", 'AbInd', col, flags=re.IGNORECASE) for col in intensity_cols]
            
            for col, new_col in zip(intensity_cols, new_intensity_cols):
                MQ_data[new_col] = MQ_data[col] / MQ_data['Mol. weight [kDa]']
        else:
            raise ValueError("The column 'Mol. weight [kDa]' could not be found!")

    # summary.txt
    if file_type == 'summary':

        row_index = MQ_data.notna().sum(axis=1) > 2
        row_index.iloc[-1] = False

        MQ_data = MQ_data[row_index]

    # evidence.txt
    # if file_type == 'evidence':
        # print('evidence.txt')
        # columns: 
        #   *Type: The type of the feature. 
        #           ‘MSMS’ - for an MS/MS spectrum without an MS1 isotope pattern assigned.
        #           ‘ISO-MSMS’ - MS1 isotope cluster identified by MS/MS.
        #           ‘MULTI-MSMS’ - MS1 labeling cluster identified by MS/MS.
        #           ‘MULTI-SECPEP’ - MS1 labeling cluster identified by MS/MS as second peptide.
        #           ‘MULTI-MATCH’ - MS1 labeling cluster identified by matching between runs.
        #           In case of label-free data there is no difference between ‘MULTI’ and ‘ISO’.
        #   *Modified sequence: Sequence representation including the post-translational modifications 
        #                       (abbreviation of the modification in brackets before the modified AA).
        #                       The sequence is always surrounded by underscore characters (’_’).

    return MQ_data


# 1. proteinGroups.txt
def get_protegroups(file_path):

    MQ_data = read(file_path, 'protgroup', filter_type='Reverse')

    # 'intensity'
    intensity_hlm_exp_cols = [col for col in MQ_data.columns if re.match(r"^Intensity [HLM] \S+$", col)]
    intensity_exp_cols = [col for col in MQ_data.columns if re.match(r"^Intensity \S+$", col)]
    intensity_cols = [col for col in MQ_data.columns if re.match(r"Intensity$", col)]

    intensity_columns = []

    if len(intensity_hlm_exp_cols):
        intensity_columns = intensity_hlm_exp_cols
    elif len(intensity_exp_cols):
        intensity_columns = intensity_exp_cols
    elif len(intensity_cols):
        intensity_columns = intensity_cols
    
    # 1: PG: ~Contaminants
    contaminant_percent_dict = pg_contaminants(MQ_data, intensity_columns)
    
    # 2. PG: ~Intensity distribution
    intensity_distr_dict = pg_intensity_distr(MQ_data, intensity_columns)

    ## 'LFQ intensity'
    lfq_intensity_hlm_exp_cols = [col for col in MQ_data.columns if re.match(r"^LFQ intensity [HLM] \S+$", col)]
    lfq_intensity_exp_cols = [col for col in MQ_data.columns if re.match(r"^LFQ intensity \S+$", col)]

    lfq_intensity_columns = []

    if len(lfq_intensity_hlm_exp_cols):
        lfq_intensity_columns = lfq_intensity_hlm_exp_cols
    elif len(lfq_intensity_exp_cols):
        lfq_intensity_columns = lfq_intensity_exp_cols

    # LFQ
    if lfq_intensity_columns:
        lfq_intensity_distr = pg_intensity_distr(MQ_data, lfq_intensity_columns)
    else:
        lfq_intensity_distr = None

    # PCA
    if len(intensity_columns) > 1:
        raw_intensity_pca = pg_PCA(MQ_data, intensity_columns)
    else:
        raw_intensity_pca = None

    if len(lfq_intensity_columns) > 1:
        lfq_intensity_pca = pg_PCA(MQ_data, lfq_intensity_columns)
    else:
        lfq_intensity_pca = None

    return {
        'pg_contaminant': contaminant_percent_dict, 
        'pg_intensity_distri': intensity_distr_dict, 
        # 'intensity_cols': intensity_columns,
        # 'lfq_intensity_cols': lfq_intensity_columns,
        # 'pg_data': MQ_data, 
        'pg_lfq_intensity_distri': lfq_intensity_distr, 
        'raw_intensity_pca': raw_intensity_pca, 
        'lfq_intensity_pca': lfq_intensity_pca
    }

# 1-1. PG:~Contaminants
def pg_contaminants(MQ_data, intensity_cols):
    
    df1 = MQ_data[intensity_cols].sum().to_frame().reset_index()
    df1.columns = ['group', 'total_intensity']
    
    df2 = MQ_data[MQ_data['Potential contaminant'] == '+'][intensity_cols].sum().to_frame().reset_index()
    df2.columns = ['group', 'contaminant_total_intensity']

    result_df = pd.merge(df1, df2, on='group', how='inner')
    result_df['contaminant_percent'] = result_df['contaminant_total_intensity'] / result_df['total_intensity'] * 100.00

    result_dict = dict()
    for k, v in dict(zip(result_df['group'], result_df['contaminant_percent'])).items():
        result_dict[k] = {'Potential Contaminants': v}

    return result_dict


# 1-2. PG: ~Intensity distribution
def pg_intensity_distr(MQ_data, intensity_cols):

    if any(column not in MQ_data.columns for column in ['Potential contaminant']):
        return None
    
    raw_df = MQ_data[intensity_cols]
    contaminant_df = MQ_data[MQ_data['Potential contaminant'] == '+'][intensity_cols]

    # Histogram plot data
    def histogram_fun(intensity_df):
        histogram_dict = {}
        for col in intensity_df.columns:
            plot_data = Histogram('Intensity distribution', plot_category='range', 
                                    breaks=[0, 10, 100, 300, 500, 700, 900, 1000, 3000, 6000, 10000])

            for intensity in intensity_df[col].to_list():
                plot_data.addValue(intensity)
            
            plot_data.to_dict()
            histogram_dict[col] = plot_data.dict['data']
        return histogram_dict
    
    histogram_data = [histogram_fun(raw_df), histogram_fun(contaminant_df)]

    # Take the logarithm and remove zero values
    def box_fun(raw_df):
        # log_df = raw_df.map(lambda x: 1 if (pd.isna(x) or x == 0) else x)
        # log_df = log_df.map(lambda x: np.log2(x)).reset_index(drop = True)
        log_df = raw_df.applymap(lambda x: 1 if (pd.isna(x) or x == 0) else x)
        # log_df = log_df.applymap(lambda x: np.log2(x)).reset_index(drop = True)
        log_df = np.log2(log_df).reset_index(drop=True)

        log_df_dict = log_df.to_dict(orient='list')
        log_df_dict = {key: [value for value in values if value != 0] for key, values in log_df_dict.items()}

        return log_df_dict

    # Box plot data
    boxplot_data = [box_fun(raw_df), box_fun(contaminant_df)]

    integrate_data = {'histogram': histogram_data, 
                      'box': boxplot_data}

    return integrate_data

# 1-3.proteinGroups.txt: PCA
def pg_PCA(pg_data, cols_name):

    if any(column not in pg_data.columns for column in ['Potential contaminant']):
        return None

    pg_data = pg_data[pg_data['Potential contaminant'] != '+'].copy()
    pca_df = pg_data[cols_name].copy().T

    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(pca_df)

    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_data)

    pca_result_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
    pca_result_df['Group'] = pca_df.index

    pca_dict = {}
    for raw_name, group in pca_result_df.groupby('Group'):
        pca_dict[raw_name] = {'x': group.iloc[0, 0], 'y': group.iloc[0, 1]}

    return pca_dict


# 2. summary.txt
def get_summary(file_path):

    MQ_data = read(file_path, 'summary')

    if any(column not in MQ_data.columns for column in ['Raw file', 'MS/MS Identified [%]']):
        return None
    
    if all(MQ_data['MS/MS Identified [%]'] == 0):
        return None

    msms_identified = dict()
    for k, v in dict(zip(MQ_data['Raw file'], MQ_data['MS/MS Identified [%]'])).items():
        msms_identified[k] = {'MS/MS Identified': v}

    return msms_identified


# 3. evidence.txt
def get_evidence(file_path):

    MQ_data = read(file_path, filter_type='Reverse')

    if not all(column in MQ_data.columns for column in ['Type']):
        raise ValueError('Missing required columns (#Type) in "evidence.txt"!')

    MQ_data['is_transferred'] = MQ_data['Type'] == 'MULTI-MATCH'

    evidence_df = MQ_data[MQ_data['Type'] != 'MULTI-MATCH'].copy()
    evidence_df_tf = MQ_data[MQ_data['Type'] == 'MULTI-MATCH'].copy()

    # Top Contaminants per Raw file
    top_cont_dict = evidence_top_contaminants(evidence_df, Top_N=5)

    # peptide intensity distribution
    peptide_intensity_dict = evidence_peptide_intensity(evidence_df)

    # Distribution of precursor charges
    charge_counts_dict = evidence_charge_distribution(evidence_df)

    # Modifications per Raw file
    modified_dict = evidence_modified(evidence_df)

    # rt_count_dict
    rt_count_dict = evidence_rt_count(evidence_df)

    # Peak width over RT
    peak_rt_dict = evidence_peak_width_rt(evidence_df)

    # Oversampling
    oversampling_dict = evidence_oversampling(evidence_df)

    # Uncalibrated Mass Error
    uncalibrated_mass_error = evidence_uncalibrated_mass_error(evidence_df)

    # Uncalibrated Mass Error
    calibrated_mass_error = evidence_calibrated_mass_error(evidence_df)

    # Peptide ID count
    peptide_id_count = evidence_peptide_count(evidence_df, evidence_df_tf)

    # ProteinGroups count
    protein_group_count = evidence_protein_count(evidence_df, evidence_df_tf)

    return {
        'top_contaminants': top_cont_dict, 
        'peptide_intensity': peptide_intensity_dict, 
        'charge_counts': charge_counts_dict, 
        'modified_percentage': modified_dict, 
        'rt_counts': rt_count_dict, 
        # 'all_evidence': MQ_data, 
        'evidence_df': evidence_df, 
        # 'evidence_df_tf': evidence_df_tf, 
        'peak_rt': peak_rt_dict, 
        'oversampling': oversampling_dict, 
        'uncalibrated_mass_error': uncalibrated_mass_error, 
        'calibrated_mass_error': calibrated_mass_error, 
        'peptide_id_count': peptide_id_count, 
        'protein_group_count': protein_group_count
        }

# 3-1. evidence.txt: Top Contaminants per Raw file
def evidence_top_contaminants(evidence_df, Top_N):

    if any(column not in evidence_df.columns for column in ['Potential contaminant', 
                                                            'Proteins', 
                                                            'Intensity', 
                                                            'Raw file']):
        return None
    
    evidence_data = evidence_df.copy()

    if 'Protein Names' in evidence_data.columns:
        evidence_data['protein_name'] = evidence_data['Protein Names'].combine_first(evidence_data['Proteins'])
    else:
        evidence_data['protein_name'] = evidence_data['Proteins']

    sum_intensity = evidence_data['Intensity'].sum()

    contaminant_df = evidence_data[evidence_data['Potential contaminant'] == '+']

    contaminant_count = contaminant_df.groupby('protein_name')['Intensity'].sum() / sum_intensity * 100
    contaminant_count = contaminant_count.sort_values(ascending=False)
    top_contaminant = list(contaminant_count.head(Top_N).index)

    contaminant_df.loc[~contaminant_df['protein_name'].isin(top_contaminant), 'protein_name'] = 'Other'

    intensity_per_file = evidence_data.groupby('Raw file', as_index=False)['Intensity'].sum()
    intensity_per_file.rename(columns={'Intensity': 'total_intensity'}, inplace=True)

    intensity_per_file_protein = contaminant_df.groupby(['Raw file', 'protein_name'], as_index=False)['Intensity'].sum()
    intensity_per_file_protein.rename(columns={'Intensity': 'contaminant_intensity'}, inplace=True)

    intensity_per_file_protein = pd.merge(intensity_per_file_protein, intensity_per_file, on='Raw file')

    intensity_per_file_protein['intensity_percent'] = intensity_per_file_protein['contaminant_intensity'] / intensity_per_file_protein['total_intensity'] * 100

    # def shorten_protein_name(name):
    #     p_split = name.split(';')
    #     p_split_s = [s[: 13] + '..' if len(s) > 15 else s for s in p_split[: min(2, len(p_split))]]
    #     result = ';'.join(p_split_s)
    #     if len(p_split) > 2:
    #         result += ';..'
    #     return result

    # intensity_per_file_protein['processed_protein_name'] = intensity_per_file_protein['protein_name'].apply(shorten_protein_name)

    top_contaminant_dict = {}

    plot_dict = {}
    for raw_file, group in intensity_per_file_protein.groupby('Raw file'):
        plot_dict[raw_file] = group[['protein_name', 'intensity_percent']].set_index('protein_name')['intensity_percent'].to_dict()

    top_contaminant_dict['plot_data'] = plot_dict
    top_contaminant_dict['cats'] = list(intensity_per_file_protein['protein_name'].unique())

    return top_contaminant_dict

# 3-2. evidence.txt: peptide intensity distribution
def evidence_peptide_intensity(evidence_df):

    if any(column not in evidence_df.columns for column in ['Potential contaminant', 
                                                            'Intensity', 
                                                            'Raw file']):
        return None
    
    evidence_data = evidence_df.copy()
    contaminant_df = evidence_data[evidence_data['Potential contaminant'] == '+'].copy()

    # Histogram plot data
    def histogram_fun(intensity_df):
        histogram_dict = {}
        for raw_file, group in intensity_df.groupby('Raw file'):
            plot_data = Histogram('Peptide intensity distribution', plot_category='range', 
                                    breaks=[0, 10, 100, 300, 500, 700, 900, 1000, 3000, 6000, 10000])
            for intensity_i in group['Intensity']:
                plot_data.addValue(intensity_i)
            plot_data.to_dict()
            histogram_dict[raw_file] = plot_data.dict['data']
        return histogram_dict

    histogram_data = [histogram_fun(evidence_data), histogram_fun(contaminant_df)]

    # Take the logarithm and remove zero values
    def box_fun(intensity_df):

        intensity_df['intensity_processed'] = intensity_df['Intensity'].apply(lambda x: 1 if (pd.isna(x) or x == 0) else x)
        # intensity_df['log_intensity_processed'] = intensity_df['intensity_processed'].apply(lambda x: np.log2(x))
        intensity_df['log_intensity_processed'] = np.log2(intensity_df['intensity_processed'])

        box_dict = {}
        for raw_file, group in intensity_df.groupby('Raw file'):
            log_intensity = group['log_intensity_processed']
            box_dict[raw_file] = list(log_intensity[log_intensity != 0])
        return box_dict
    
    # Box plot data
    boxplot_data = [box_fun(evidence_data), box_fun(contaminant_df)]

    integrate_data = {'histogram': histogram_data, 
                      'box': boxplot_data}

    return integrate_data

# 3-3.evidence.txt: charge distribution
def evidence_charge_distribution(evidence_data):

    if any(column not in evidence_data.columns for column in ['Potential contaminant', 
                                                              'Charge', 
                                                              'Raw file']):
        return None
    
    evidence_data = evidence_data[evidence_data['Potential contaminant'] != '+'].copy()

    charge_counts = evidence_data.groupby('Raw file')['Charge'].value_counts().reset_index(name='count')

    plot_dict = {}
    for raw_file, group in charge_counts.groupby('Raw file'):
        charge_counts_sorted = group.sort_values(by='Charge')
        charge_counts_sorted['Charge'] = charge_counts_sorted['Charge'].astype(str)
        plot_dict[raw_file] = dict(zip(charge_counts_sorted['Charge'], charge_counts_sorted['count']))

    charge_dict = {}
    charge_dict['plot_data'] = plot_dict
    charge_dict['cats'] = list(map(str, sorted(charge_counts['Charge'].unique())))

    return charge_dict

# 3-4.evidence.txt: Modifications per Raw file
def evidence_modified(evidence_data):

    if any(column not in evidence_data.columns for column in ['Potential contaminant', 
                                                              'Modifications', 
                                                              'Raw file']):
        return None
    
    evidence_data = evidence_data[evidence_data['Potential contaminant'] != '+'].copy()

    def group_percentage(group):
        counts = group['Modifications'].str.split(',').explode().value_counts()
        percentage_df = (counts / len(group['Modifications']) * 100).reset_index()
        percentage_df.columns = ['Modifications', 'Percentage']

        # Modified (Total)
        percentage_df.loc[percentage_df['Modifications'] == 'Unmodified', 'Percentage'] = 100 - percentage_df.loc[percentage_df['Modifications'] == 'Unmodified', 'Percentage']
        percentage_df.loc[percentage_df['Modifications'] == 'Unmodified', 'Modifications'] = 'Modified (Total)'

        return percentage_df

    plot_dict = {}
    modified_cats = []

    for raw_file, group in evidence_data.groupby('Raw file'):
        group_processed = group_percentage(group)
        plot_dict[raw_file] = dict(zip(group_processed['Modifications'], group_processed['Percentage']))
        modified_cats.extend(group_processed['Modifications'])

    modified_dict = {}
    modified_dict['plot_data'] = plot_dict
    modified_dict['cats'] = list(sorted(modified_cats, key=lambda x: (x == 'Modified (Total)', x)))

    return modified_dict

# 3-5.evidence.txt: IDs over RT
def evidence_rt_count(evidence_data):

    if any(column not in evidence_data.columns for column in ['Potential contaminant', 
                                                              'Retention time', 
                                                              'Raw file']):
        return None

    evidence_data = evidence_data[evidence_data['Potential contaminant'] != '+'].copy()

    rt_range = [evidence_data['Retention time'].min(), evidence_data['Retention time'].max()]

    def hist_compute(rt_list, rt_range):
        counts, bin_edges = np.histogram(rt_list, bins=np.arange(rt_range[0] - 3, rt_range[1] + 3, 1))
        bin_mid = (bin_edges[: -1] + bin_edges[1: ]) / 2
        rt_counts = pd.DataFrame({'retention_time': bin_mid, 'counts': counts})
        
        return dict(zip(rt_counts['retention_time'], rt_counts['counts']))

    rt_count_dict = {}
    for raw_file, group in evidence_data.groupby('Raw file'):
        rt_count_dict[raw_file] = hist_compute(group['Retention time'], rt_range)

    return rt_count_dict

# 3-6.evidence.txt: Peak width over RT
def evidence_peak_width_rt(evidence_data):

    if any(column not in evidence_data.columns for column in ['Potential contaminant', 
                                                              'Retention length', 
                                                              'Retention time', 
                                                              'Raw file']):
        return None

    evidence_data = evidence_data[evidence_data['Potential contaminant'] != '+'][['Retention length', 'Retention time', 'Raw file']].copy()

    rt_range = [evidence_data['Retention time'].min(), evidence_data['Retention time'].max()]
    breaks = np.arange(rt_range[0] - 3, rt_range[1] + 3, 1)

    def RT_RL_compute(group_df):
        group_df['bin'] = np.digitize(group_df['Retention time'], breaks, right=False)
        bin_group = group_df.groupby('bin')
        rt_rl_df = bin_group['Retention length'].median().reset_index()
        rt_rl_df.rename(columns={'Retention length': 'median_RL'}, inplace=True)
        rt_rl_df['bin_RT'] = breaks[rt_rl_df['bin'] - 1]

        return rt_rl_df

    peak_width_rt_dict = {}
    for raw_file, group in evidence_data.groupby('Raw file'):
        peak_width_rt_dict[raw_file] = dict(zip(RT_RL_compute(group)['bin_RT'], RT_RL_compute(group)['median_RL']))

    return peak_width_rt_dict

# 3-7.evidence.txt: Oversampling (MS/MS counts per 3D-peak)
def evidence_oversampling(evidence_data):

    if any(column not in evidence_data.columns for column in ['Potential contaminant', 
                                                              'MS/MS Count', 
                                                              'Raw file']):
        return None

    evidence_data = evidence_data[evidence_data['Potential contaminant'] != '+'][['MS/MS Count', 'Raw file']].copy()

    evidence_data['MS/MS Count'] = evidence_data['MS/MS Count'].apply(lambda x: '>=3' if x >= 3 else x)
    oversampling_df = evidence_data.groupby('Raw file')['MS/MS Count'].value_counts().reset_index()
    oversampling_df['MS/MS Count'] = oversampling_df['MS/MS Count'].astype(str)

    plot_dict = {}
    for raw_file, group in oversampling_df.groupby('Raw file'):
        group['fraction'] = group['count'] / group['count'].sum() * 100
        plot_dict[raw_file] = dict(zip(group['MS/MS Count'], group['fraction']))

    oversampling = {}
    oversampling['plot_data'] = plot_dict
    # oversampling['plot_data'] = {raw_file: dict(zip(group['MS/MS Count'], group['count'])) for raw_file, group in oversampling_df.groupby('Raw file')}
    oversampling['cats'] = list(oversampling_df['MS/MS Count'].unique())

    return oversampling

# 3-8.evidence.txt: Uncalibrated mass error
def evidence_uncalibrated_mass_error(evidence_data):

    if any(column not in evidence_data.columns for column in ['Potential contaminant', 
                                                              'Uncalibrated Mass Error [ppm]', 
                                                              'Raw file']):
        return None
    
    evidence_data = evidence_data[evidence_data['Potential contaminant'] != '+'].copy()

    uncalibrated_mass_error = {}
    for raw_file, group in evidence_data.groupby('Raw file'):
        mass_error = list(group['Uncalibrated Mass Error [ppm]'].map(lambda x: 0 if pd.isna(x) else x))

        uncalibrated_mass_error[raw_file] = [value for value in mass_error if value != 0]

    return uncalibrated_mass_error

# 3-8.evidence.txt: Calibrated mass error
def evidence_calibrated_mass_error(evidence_data):
    
    if any(column not in evidence_data.columns for column in ['Potential contaminant', 
                                                              'Mass Error [ppm]', 
                                                              'Raw file']):
        return None
    
    evidence_data = evidence_data[evidence_data['Potential contaminant'] != '+'].copy()

    calibrated_mass_error = {}
    for raw_file, group in evidence_data.groupby('Raw file'):
        mass_error = list(group['Mass Error [ppm]'].map(lambda x: 0 if pd.isna(x) else x))

        calibrated_mass_error[raw_file] = [value for value in mass_error if value != 0]

    return calibrated_mass_error

# 3-9.evidence.txt: Peptide ID count
def evidence_peptide_count(evidence_df, evidence_df_tf):

    if any(column not in evidence_df.columns for column in ['Potential contaminant', 
                                                            'Modified sequence', 
                                                            'is_transferred', 
                                                            'Raw file']):
        return None
    
    if any(column not in evidence_df_tf.columns for column in ['Potential contaminant', 
                                                               'Modified sequence', 
                                                               'is_transferred', 
                                                               'Raw file']):
        return None
    
    evidence_data = evidence_df.copy()
    evidence_data_tf = evidence_df_tf.copy()

    evidence_data = evidence_data[evidence_data['Potential contaminant'] != '+'].copy()
    evidence_data_tf = evidence_data_tf[evidence_data_tf['Potential contaminant'] != '+'].copy()

    required_cols = ['Raw file', 'is_transferred', 'Modified sequence']
    evid_df = pd.concat([evidence_data[required_cols], evidence_data_tf[required_cols]], axis=0, ignore_index=True)

    def get_peptide_counts(evd_df):
        
        peptide_counts = pd.DataFrame()
        for raw_file, group in evd_df.groupby('Raw file'):
            pep_set_genuineUnique = group[~group['is_transferred']]['Modified sequence'].unique()
            pep_set_allMBRunique = group[group['is_transferred']]['Modified sequence'].unique()
            pep_count_GenAndMBR = len(set(pep_set_genuineUnique).intersection(pep_set_allMBRunique))
            pep_count_newMBR = len(pep_set_allMBRunique) - pep_count_GenAndMBR
            pep_count_onlyGenuine = len(pep_set_genuineUnique) - pep_count_GenAndMBR

            if any(evd_df['is_transferred']):
                file_peptide_counts = pd.DataFrame({
                    'Raw file': [raw_file, raw_file, raw_file], 
                    'counts': [pep_count_onlyGenuine, pep_count_GenAndMBR, pep_count_newMBR],
                    'category': ['Genuine (Exclusive)', 'Genuine + Transferred', 'Transferred (Exclusive)'],
                    'MBRgain': [None, None, pep_count_newMBR / (pep_count_onlyGenuine + pep_count_GenAndMBR) * 100]
                })
                categorys = ['Genuine (Exclusive)', 'Genuine + Transferred', 'Transferred (Exclusive)']
            else:
                file_peptide_counts = pd.DataFrame({
                    'Raw file': [raw_file], 
                    'counts': [pep_count_onlyGenuine],
                    'category': ['Genuine'],
                    'MBRgain': [None]
                })
                categorys = ['Genuine']

            peptide_counts = pd.concat([peptide_counts, file_peptide_counts], axis=0, ignore_index=True)
        return peptide_counts, categorys

    peptide_counts_df, cats = get_peptide_counts(evid_df)

    plot_data = {}
    for raw_file, group in peptide_counts_df.groupby('Raw file'):
        plot_data[raw_file] = dict(zip(group['category'], group['counts']))

    peptide_id_count = {}
    peptide_id_count['plot_data'] = plot_data
    peptide_id_count['cats'] = cats
    peptide_id_count['title_value'] = 'MBR gain: +{}%'.format(round(peptide_counts_df['MBRgain'].mean(), 2)) if any(evid_df['is_transferred']) else ''
    
    return peptide_id_count

# 3-10.evidence.txt: ProteinGroups count
def evidence_protein_count(evidence_df, evidence_df_tf):

    if any(column not in evidence_df.columns for column in ['Potential contaminant', 
                                                            'Protein group IDs', 
                                                            'is_transferred', 
                                                            'Raw file']):
        return None
    
    if any(column not in evidence_df_tf.columns for column in ['Potential contaminant', 
                                                               'Protein group IDs', 
                                                               'is_transferred', 
                                                               'Raw file']):
        return None

    evidence_data = evidence_df.copy()
    evidence_data_tf = evidence_df_tf.copy()

    evidence_data = evidence_data[evidence_data['Potential contaminant'] != '+'].copy()
    evidence_data_tf = evidence_data_tf[evidence_data_tf['Potential contaminant'] != '+'].copy()

    required_cols = ['Raw file', 'is_transferred', 'Protein group IDs']
    evid_df = pd.concat([evidence_data[required_cols], evidence_data_tf[required_cols]], axis=0, ignore_index=True)

    def get_protein_group_counts(evd_df):
        
        protein_group_counts = pd.DataFrame()
        for raw_file, group in evd_df.groupby('Raw file'):

            group['protein_group_mtd'] = group['Protein group IDs'] + '_' + group['is_transferred'].astype(str)
            duplicated_protein_group = group[~group['protein_group_mtd'].duplicated()]

            protein_groups = duplicated_protein_group['Protein group IDs'].apply(lambda x: x.split(';'))
            protein_group_genuine_unique = protein_groups[~duplicated_protein_group['is_transferred']].explode().unique()
            protein_group_MBR_unique = protein_groups[duplicated_protein_group['is_transferred']].explode().unique()
            protein_group_Gen_and_MBR = len(set(protein_group_genuine_unique).intersection(protein_group_MBR_unique))
            protein_group_MBR_only = len(protein_group_MBR_unique) - protein_group_Gen_and_MBR
            protein_group_Genuine_only = len(protein_group_genuine_unique) - protein_group_Gen_and_MBR

            if any(evd_df['is_transferred']):
                file_protein_group_counts = pd.DataFrame({
                    'Raw file': [raw_file, raw_file, raw_file], 
                    'counts': [protein_group_Genuine_only, protein_group_Gen_and_MBR, protein_group_MBR_only],
                    'category': ['Genuine (Exclusive)', 'Genuine + Transferred', 'Transferred (Exclusive)'],
                    'MBRgain': [None, None, protein_group_MBR_only / (protein_group_Genuine_only + protein_group_Gen_and_MBR) * 100]
                })
                categorys = ['Genuine (Exclusive)', 'Genuine + Transferred', 'Transferred (Exclusive)']
            else:
                file_protein_group_counts = pd.DataFrame({
                    'Raw file': [raw_file], 
                    'counts': [protein_group_Genuine_only],
                    'category': ['Genuine'],
                    'MBRgain': [None]
                })
                categorys = ['Genuine']

            protein_group_counts = pd.concat([protein_group_counts, file_protein_group_counts], axis=0, ignore_index=True)
        return protein_group_counts, categorys

    protein_group_counts_df, cats = get_protein_group_counts(evid_df)

    plot_data = {}
    for raw_file, group in protein_group_counts_df.groupby('Raw file'):
        plot_data[raw_file] = dict(zip(group['category'], group['counts']))

    protein_group_count = {}
    protein_group_count['plot_data'] = plot_data
    protein_group_count['cats'] = cats
    protein_group_count['title_value'] = 'MBR gain: +{}%'.format(round(protein_group_counts_df['MBRgain'].mean(), 2)) if any(evid_df['is_transferred']) else ''
    
    return protein_group_count


# 4.msms.txt
def get_msms(file_path, evidence_df=None):

    MQ_data = read(file_path)

    if evidence_df is None:
        return {'missed_cleavages': None}

    # Missed cleavages per Raw file
    missed_cleavages = msms_missed_cleavages(MQ_data, evidence_df)

    return {
        # 'MQ_data': MQ_data, 
        'missed_cleavages': missed_cleavages
    }

# 4-1.msms.txt: Missed cleavages per Raw file
def msms_missed_cleavages(msms_df, evidence_df):

    if any(column not in msms_df.columns for column in ['Evidence ID', 
                                                        'Missed cleavages', 
                                                        'Raw file']):
        return None
    
    if any(column not in evidence_df.columns for column in ['Potential contaminant', 'id']):
        return None
    
    # excludes contaminants & excludes Type == 'MULTI-MATCH'
    evidence_not_contaminant_ids = evidence_df[evidence_df['Potential contaminant'] != '+']['id'].copy()
    not_contaminant_index = msms_df['Evidence ID'][msms_df['Evidence ID'].isin(evidence_not_contaminant_ids)].index

    msms_not_contaminant = msms_df.loc[not_contaminant_index, ['Raw file', 'Missed cleavages', 'Evidence ID']]
    msms_not_contaminant['Missed cleavages'] = msms_not_contaminant['Missed cleavages'].astype('str')

    plot_dict = {}
    for raw_file, group in msms_not_contaminant.groupby('Raw file'):

        missed_cleavages_df = pd.DataFrame(group['Missed cleavages'].value_counts().reset_index())
        missed_cleavages_df['percentage'] = missed_cleavages_df['count'] / missed_cleavages_df['count'].sum() * 100
        
        plot_dict[raw_file] = dict(zip(missed_cleavages_df['Missed cleavages'], missed_cleavages_df['percentage']))

    missed_cleavages_dict = {}
    missed_cleavages_dict['plot_data'] = plot_dict
    missed_cleavages_dict['cats'] = sorted(list(msms_not_contaminant['Missed cleavages'].unique()))

    return missed_cleavages_dict


# 5.msScans.txt
def get_msmsScans(file_path):

    MQ_data = read(file_path)

    # if not all(column in MQ_data.columns for column in ['Scan event number', 'Retention time', 'Raw file', 'Scan event number']):
    #     raise ValueError('Missing required columns in "msmsScans.txt"!')
    
    # TODO check 'Scan event number'

    MQ_data['round_RT'] = MQ_data['Retention time'].apply(lambda x: round(x / 2) * 2)

    # Ion Injection Time over RT
    ion_injec_time_RT = msmsScans_ion_injec_time_RT(MQ_data)

    # TopN over RT
    top_overRT = msmsScans_top_over_RT(MQ_data)

    # TopN
    top_n = msmsScans_top_N(MQ_data)

    return {
        # 'MQ_data': MQ_data, 
        'ion_injec_time_RT': ion_injec_time_RT, 
        'top_n': top_n,
        'top_overRT': top_overRT
        }

# 5-1.msmsScans.txt: Ion Injection Time over RT
def msmsScans_ion_injec_time_RT(msmsScans_df):

    if any(column not in msmsScans_df.columns for column in ['Raw file', 'Ion injection time']):
        return None
    
    if msmsScans_df['Ion injection time'].isna().all():
        return None
    
    median_ion_injec_time_df = msmsScans_df.groupby(['Raw file', 'round_RT'])['Ion injection time'].median().reset_index()
    median_ion_injec_time_df = median_ion_injec_time_df.rename(columns={'Ion injection time': 'median_ion_injection_time'})

    mean_ion_injec_time_df = msmsScans_df.groupby(['Raw file'])['Ion injection time'].mean().reset_index()
    mean_ion_injec_time_df = mean_ion_injec_time_df.rename(columns={'Ion injection time': 'mean_ion_injection_time'})
    mean_ion_injec_time_df['int_mean_ion_injection_time'] = mean_ion_injec_time_df['mean_ion_injection_time'].apply(lambda x: int(x) if not pd.isna(x) else 0)
    mean_ion_injec_time_df['int_mean_ion_injection_time'] = mean_ion_injec_time_df['int_mean_ion_injection_time'].astype(str)
    mean_ion_injec_time_df['raw_file_mean_ion_time'] = mean_ion_injec_time_df['Raw file'] + \
         ' (~' + mean_ion_injec_time_df['int_mean_ion_injection_time'].astype(str) + 'ms)'

    result_df = pd.merge(median_ion_injec_time_df, mean_ion_injec_time_df[['Raw file', 'raw_file_mean_ion_time']], on='Raw file')

    ion_injec_time_dict = {}
    for raw_file, group in result_df.groupby('raw_file_mean_ion_time'):

        ion_injec_time_dict[raw_file] = dict(zip(group['round_RT'], group['median_ion_injection_time']))

    return ion_injec_time_dict

# 5-2.msmsScans.txt: TopN over RT
def msmsScans_top_over_RT(msmsScans_df):

    if any(column not in msmsScans_df.columns for column in ['Raw file', 
                                                             'Retention time', 
                                                             'Scan event number']):
        return None
    
    # Scan event number: 
    #   This number indicates which MS/MS scan this one is in the consecutive order of the MS/MS scans that are acquired after an MS scan.
    # Retention time: 
    #   Time point along the elution profile at which the MS/MS data was recorded.
    msmsScans_data = msmsScans_df[['Raw file', 'Retention time', 'round_RT', 'Scan event number']].copy()
    msmsScans_data = msmsScans_data.sort_values(by=['Raw file', 'Retention time'])

    def find_local_maxima(list_data):
        local_max_indices = np.zeros(len(list_data), dtype=bool)
        for i in range(0, len(list_data) - 1):
            if np.isnan(list_data[i]):
                continue
            if list_data[i] > list_data[i+1] and list_data[i] >= 0:
                local_max_indices[i] = True
        return local_max_indices

    local_max_SE_number_df = pd.DataFrame()
    for raw_file, group in msmsScans_data.groupby('Raw file'):
        local_max_SE_number_df = pd.concat([local_max_SE_number_df, 
                                            group[find_local_maxima(list(group['Scan event number']))]], 
                                            axis=0, ignore_index=True)

    local_max_SE_number_df = local_max_SE_number_df.rename(columns={'Scan event number': 'local_max_scan_event_number'})
    median_SE_number_df = local_max_SE_number_df.groupby(['Raw file', 'round_RT'])['local_max_scan_event_number'].median().reset_index()

    scan_event_number_dict = {}
    for raw_file, group in median_SE_number_df.groupby('Raw file'):
        scan_event_number_dict[raw_file] = dict(zip(group['round_RT'], group['local_max_scan_event_number']))

    return scan_event_number_dict

# 5-3.msmsScans.txt: TopN
def msmsScans_top_N(msmsScans_df):

    if any(column not in msmsScans_df.columns for column in ['Raw file', 
                                                             'Scan event number']):
        return None
    
    msmsScans_data = msmsScans_df[['Scan event number', 'Raw file']].copy()

    while True:
        scan_event_index = 1 + np.where(np.diff(msmsScans_data['Scan event number']) > 1)[0]
        if len(scan_event_index) == 0: 
            break
        msmsScans_data.loc[scan_event_index, 'Scan event number'] -= 1

    file_SE_count_df = msmsScans_data[['Scan event number', 'Raw file']].value_counts().reset_index(name='count')
    max_SE_number = max(file_SE_count_df['Scan event number'])

    def process_group(group_df, max_scan_event_number, raw_file):
        event_count = group_df['count'].values

        if not all(event_count[i] >= event_count[i + 1] for i in range(len(event_count) - 1)):
            raise ValueError("Scan event distribution is not monotonically increasing!")
        if max(group_df['Scan event number']) != len(group_df):
            raise ValueError("Scan event distribution has unexpected holes...!")
        
        event_pre = np.append(event_count[1: ], 0)
        event_diff = event_count - event_pre

        SE_number = group_df['Scan event number'].values
        if max(SE_number) < max_scan_event_number:
            event_diff = list(event_diff) + [0] * (max_scan_event_number - max(SE_number))
            SE_number = list(SE_number) + list(range(max(SE_number) + 1, max_scan_event_number + 1))

        result_df = pd.DataFrame({'Raw file': raw_file, 'Scan event number': SE_number, 'count': event_diff})
        return result_df

    file_SE_count_ratio = pd.DataFrame()
    for raw_file, group in file_SE_count_df.groupby('Raw file'):
        file_SE_count_ratio = pd.concat([file_SE_count_ratio, process_group(group, max_SE_number, raw_file)], axis=0, ignore_index=True)

    file_SE_count_ratio['Scan event number'] = file_SE_count_ratio['Scan event number'].astype(str)

    plot_dict = {}
    for raw_file, group in file_SE_count_ratio.groupby('Raw file'):
        plot_dict[raw_file] = dict(zip(group['Scan event number'], group['count']))

    SE_count_dict = {}
    SE_count_dict['plot_data'] = plot_dict
    SE_count_dict['cats'] = [str(x) for x in reversed(file_SE_count_ratio['Scan event number'].unique())]

    return SE_count_dict


# 6.parameters.txt
def get_parameters(file_path):

    MQ_data = read(file_path)

    # if not all(column in MQ_data.columns for column in ['Parameter', 'Value']):
    #     raise ValueError('Missing required columns in "parameters.txt"!')

    parameters_tb_dict = parameters_table(MQ_data)

    return {
        # 'MQ_data': MQ_data, 
        'parameters_tb_dict': parameters_tb_dict
    }

# 6-1.parameters.txt: Parameters
def parameters_table(parameters_df):

    if any(column not in parameters_df.columns for column in ['Parameter', 
                                                              'Value']):
        return None
    
    parameters_data = parameters_df[~parameters_df['Parameter'].str.startswith('AIF')]

    fasta_files = parameters_data[parameters_data['Parameter'] == 'Fasta file']['Value'].values[0]
    fasta_files = fasta_files.split(';')
    def parse_location(location):
        if '\\' in location:
            location = location.replace('\\', '/')
        return os.path.basename(location)
    fasta_file_list = [parse_location(fasta_file) for fasta_file in fasta_files]
    fasta_file_list = ';'.join(fasta_file_list)

    table_data = parameters_data.drop_duplicates(subset='Parameter', keep='first').reset_index(drop=True)
    table_data.loc[table_data['Parameter'] == 'Fasta file', 'Value'] = fasta_file_list

    parameters_dict = {}
    for index, row in table_data.iterrows():
        parameters_dict[index + 1] = row.to_dict()

    return parameters_dict