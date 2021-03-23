import pandas as pd

from pathlib import Path
from scipy import stats
import numpy as np

import math
import subprocess

from scipy.stats import norm
import statsmodels
import scipy
import streamlit as st
from datetime import date


def load_files(samples, results_dir):
    df_list = []
    cntrl_list = []
    for sample in samples:
        df_list.append(pd.read_csv(Path(results_dir)/sample/f"{sample}_counts.csv", index_col = 0))
        cntrl_list.append(pd.read_csv(Path(results_dir)/sample/f"{sample}_counts.control.csv", index_col = 0))
    return pd.concat([pd.concat(df_list), pd.concat(cntrl_list)])


def subset_experiment(df,  dnaid,  exp, col1='experiment', col2='dnaid'):
    '''
    example query string : '(exp=="TV5490A") & (dnaid == "dnaid2023")'
    '''
    query_string = f'({col1} == "{exp}") & ({col2} == "{dnaid}")'
    return df.copy().query(query_string)


# Analyze a DataSet
## Looking at acontrols

def calculate_correlation(exp_df, control_file, for_each='sampleID', how='log', cutoff=0.9):
    """
    Subset counts for control barcodes
    Calculate correlation on log counts (log), log counts, but keep 0 (log_w_0), or raw data (raw)

    """
    controls = pd.read_table(control_file, names=['barcode', 'phenotype', 'conc'])
    control_cnts = exp_df[exp_df.barcode.isin(controls.barcode)].copy()
    if how == 'raw':
        col1 = 'conc'
        col2 = 'cnt'
    else:
        control_cnts['logConc'] = np.log10(control_cnts['conc'])
        col1 = 'logConc'
        if how == 'log':
            control_cnts['logCnts'] = np.log10(control_cnts['cnt'])
        elif how == 'log_w_0':
            control_cnts['logCnts'] = np.log10(control_cnts['cnt'].replace({0: 1}))
        col2 = 'logCnts'
    corr_df = control_cnts.groupby(['phenotype', for_each])[[col1, col2]].corr()
    corr_df = corr_df.reset_index()
    corr_df = corr_df[corr_df['level_2'] == col1].drop(['level_2', col1], axis=1)
    corr_df.columns = ['phenotype', 'sampleID', 'R']

    good_samples = corr_df[(corr_df.R > cutoff) & (corr_df.phenotype == 'wt')].sampleID.values
    return corr_df, good_samples


## Filter Data

def filter_inoculum(exp_df, filter_below=0):
    filt_df = (exp_df.copy()
               .drop(['ShortName', 'locus_tag'], axis=1)
               .drop_duplicates()
               .pivot(index='barcode', columns='sampleID', values='cnt'))
    filt_df = filt_df.fillna(0)
    columns_to_filter = [f for f in filt_df.columns if 'inoculum' in f]
    filt_df = filt_df[(filt_df[columns_to_filter] >= filter_below).all(1)]
    #filt_df = filt_df[(filt_df['inoculum_d0'] >=filter_below)& (filt_df['unenriched_inoculum_d0'] >= filter_below)]
    return filt_df


def filter_samples(exp_df, good_samples):
    return exp_df.copy()[exp_df.sampleID.isin(good_samples)]


def generate_DE_dataset(exp_df, good_samples, filter_below = 1):
    sample_data = exp_df[['sampleID', 'mouse', 'day', 'organ', 'dnaid']].set_index('sampleID').drop_duplicates()
    sample_data = sample_data.loc[sample_data.index.intersection(good_samples)]
    expr_data = filter_samples(exp_df, good_samples)
    expr_data = filter_inoculum(expr_data, filter_below=filter_below)
    expr_data = expr_data[list(sample_data.index)].reset_index()
    return sample_data, expr_data


def run_command(args):
    """Run command, transfer stdout/stderr"""
    result = subprocess.run(args)
    try:
        result.check_returncode()
    except subprocess.CalledProcessError as e:
        raise e


def get_fintess_results(fitness_dir, dnaid, experiment, sdf, edf):
    sdf_path = Path(fitness_dir) / f"{dnaid}_{experiment}_sdf.csv"
    edf_path = Path(fitness_dir) / f"{dnaid}_{experiment}_edf.csv"
    sdf.to_csv(sdf_path)
    edf.set_index('barcode').to_csv(edf_path)
    r = run_command(['Rscript', 'DEseq.R', sdf_path, edf_path])
    fitness = pd.concat(
        [pd.read_table(f, sep=' ').assign(day=f.stem.split("_")[3]) for f in Path(fitness_dir).iterdir() if
         f"{dnaid}_{experiment}_fitness" in f.stem])
    n_samples = sdf.groupby('day').mouse.nunique().to_dict()
    fitness['n_samples'] = fitness.day.map(n_samples)
    fitness = fitness.reset_index().rename({'index': 'barcode'}, axis=1)
    return fitness


def sigma(lfcSE):
    return np.sqrt(lfcSE.pow(2).sum()) / len(lfcSE)


def calculate_2dist_zscore(u1, s1, u2, s2):
    return (u1 - u2) / np.sqrt((s1 ** 2) + (s2 ** 2))


def calculte_comparisons(fitness, df, control_file):
    """

    fitness: DESeq2 output, log2FoldChange value for each barcode comparing each time point with inoculum
    df: df for 1 experiment and 1 dnaid
    controls: control meta df?
    """
    days = sorted(list(fitness['day'].unique()))

    # days.remove('d0')
    controls = pd.read_table(control_file, names=['barcode', 'phenotype', 'conc'])
    # Get all entries that were mapped to a gene
    gene_bc = df[df.locus_tag.notnull()].barcode.values
    gene_df = fitness[fitness.barcode.isin(gene_bc)]  # subsetting only on barcodes present in fitness table
    # Add gene annotation to the fitness table
    gene_df = gene_df.merge(df[['barcode', 'ShortName']], how='left', on='barcode').drop_duplicates()

    # Calculate mean log2FoldChange and sigma for each gene (for all )
    gene_mean = gene_df.groupby(['ShortName', 'day']).agg(
        {'log2FoldChange': ['mean', 'median'], 'lfcSE': [sigma]}).reset_index()

    gene_mean.columns = ['gene', 'day', 'gene_FC', 'gene_FC_median', 'sigma']

    # Get all the WITS barcodes
    controls_bc = controls[controls.phenotype == 'wt'].barcode.values
    cntrl_df = fitness[fitness.barcode.isin(controls_bc)]
    # Calculate mean log2FoldChange and sigma for the control barcodes (for all barcodes)
    cntrl_mean = cntrl_df.groupby(['day']).agg({'log2FoldChange': ['mean', 'median'], 'lfcSE': [sigma]})
    cntrl_mean.columns = ['cntrl_FC', 'cntrl_FC_median', 'cntrl_sigma']
    cntrl_mean = cntrl_mean.reset_index()

    print(cntrl_mean)
    # Calculate zscore and competitive index (CI) for each gene
    gene_mean = gene_mean.merge(cntrl_mean, how='left', on='day')
    gene_mean['zscore'] = gene_mean.apply(
        lambda x: calculate_2dist_zscore(x['gene_FC'], x['sigma'], x['cntrl_FC'], x['cntrl_sigma']), axis=1)
    gene_mean['ci'] = gene_mean.apply(lambda x: 2 ** x['gene_FC'] / 2 ** x['cntrl_FC'], axis=1)

    gene_mean = gene_mean[['gene', 'day', 'zscore', 'ci']]

    # Get all barcodes that were not mapped to a gene
    others_bc = df[(df.locus_tag.isna())].barcode.values
    other_df = fitness[fitness.barcode.isin(others_bc)]
    other_df = other_df.merge(cntrl_mean, how='left', on='day')

    # Calculate zscore and CI for each barcode
    other_df['zscore'] = other_df.apply(
        lambda x: calculate_2dist_zscore(x['log2FoldChange'], x['lfcSE'], x['cntrl_FC'], x['cntrl_sigma']), axis=1)
    other_df['ci'] = other_df.apply(lambda x: 2 ** x['log2FoldChange'] / 2 ** x['cntrl_FC'], axis=1)
    other_df = other_df[['barcode', 'day', 'zscore', 'ci']].rename({'barcode': 'gene'}, axis=1)

    # Concatenate the gene and barcode results
    results = pd.concat([gene_mean, other_df])

    # Calculate p-values for the genes/barcodes
    results['pval'] = results.zscore.apply(lambda x: scipy.stats.norm.sf(abs(x)) * 2)
    # Spread results over days

    results = results.drop_duplicates()
    results = (results.pivot(index=['gene'], columns=['day'], values=['zscore', 'ci', 'pval'])
               .reset_index())

    # Rename columns
    results.columns = ['gene'] + [f'{day}_zscore' for day in days] + [f'{day}_ci' for day in days] + [f'{day}_pval' for
                                                                                                      day in days]
    # Adjust p-values for multiple testing
    for day in days:
        results[f'{day}_padj'] = \
            statsmodels.stats.multitest.multipletests(results[f'{day}_pval'], alpha=0.05, method='fdr_bh')[1]
    return results.set_index('gene')


# Final processing

def to_list(x):
    bc_list = list(x)
    if len(bc_list) == 1:
        return bc_list[0]
    return ", ".join(list(x))


def get_control_fitness(fitness, control_file):
    controls = pd.read_table(control_file, names=['barcode', 'phenotype', 'conc'])
    cnrls = controls.merge(fitness, how='left', on='barcode').assign(control='yes')
    cnrls = cnrls[['barcode', 'phenotype', 'conc', 'log2FoldChange', 'n_samples', 'day', 'control']]
    cnrls = cnrls[cnrls.day.notnull()]
    days = list(cnrls.day.unique())
    cnrls = cnrls.pivot(index=['barcode', 'phenotype', 'conc', 'control'], columns=['day'],
                        values=['log2FoldChange', 'n_samples']).reset_index()
    cnrls['barcode'] = cnrls['phenotype'] + "-" + cnrls['conc'].astype(str)
    day_columns = [f'{day}_fitness' for day in days] + [f'{day}_num_samples' for day in days]
    cnrls.columns = ['barcode', 'phenotype', 'conc', 'control'] + day_columns
    to_keep = ['barcode'] + day_columns
    cnrls = cnrls[to_keep]
    mean_fit = cnrls.groupby(['barcode']).agg(
        {d: ['mean', 'std'] for d in [f'{day}_fitness' for day in days]}).reset_index()
    mean_fit.columns = [f'{i[0]}_{i[1]}' if i[1] else f'{i[0]}' for i in mean_fit.columns]
    mean_fit = mean_fit.merge(cnrls[['barcode'] + [f'{day}_num_samples' for day in days]], how='left', on='barcode')
    return mean_fit.set_index('barcode')


def final_fitness_table(fitness, exp_df, control_file, results):
    barcode_info = exp_df[['barcode', 'locus_tag', 'ShortName', 'library']].drop_duplicates()
    fitness = fitness.merge(barcode_info, how='left', on='barcode')

    fitness_gene = fitness[(fitness.locus_tag.notnull()) & (fitness.locus_tag != '-')]
    fit_summary = fitness_gene.groupby(['library', 'ShortName', 'locus_tag', 'day', 'n_samples', ]).agg(
        {'barcode': ['count', to_list], 'log2FoldChange': ['mean', 'std']}).reset_index()
    fit_summary.columns = ['library', 'gene', 'locus', 'day', 'num_samples', 'num_barcodes', 'barcode', 'mean_fitness',
                           'std_fitness']

    fit_summary = (fit_summary.pivot(index=['gene', 'locus', 'num_barcodes', 'barcode', 'library'], columns=['day'],
                                     values=['mean_fitness', 'std_fitness', 'num_samples'])
                   .reset_index())

    fitness_no_gene = fitness[(fitness.locus_tag.isna()) | (fitness.locus_tag == '-')]
    fitness_no_gene = fitness_no_gene.groupby(['barcode', 'library', 'day', 'n_samples', ]).agg(
        {'barcode': ['count', to_list], 'log2FoldChange': ['mean', 'std']}).reset_index()
    fitness_no_gene.columns = ['gene', 'library', 'day', 'num_samples', 'num_barcodes', 'barcode', 'mean_fitness',
                               'std_fitness']

    fitness_no_gene = (fitness_no_gene.pivot(index=['gene', 'num_barcodes', 'barcode', 'library'], columns=['day'],
                                             values=['mean_fitness', 'std_fitness', 'num_samples'])
                       .reset_index())

    days = sorted(list(fitness.day.unique()))
    # days.remove('d0')
    day_columns = [f'{day}_fitness_mean' for day in days] + [f'{day}_fitness_std' for day in days] + [
        f'{day}_num_samples' for day in days]
    fit_summary.columns = ['gene', 'locus', 'num_barcodes', 'barcode', 'library'] + day_columns
    fitness_no_gene.columns = ['gene', 'num_barcodes', 'barcode', 'library'] + day_columns

    fit = pd.concat([fit_summary, fitness_no_gene])
    fit = fit.merge(exp_df[['barcode', 'sseqid', 'sstart']].drop_duplicates(), how='left', on='barcode')

    cnrls = get_control_fitness(fitness, control_file)
    cnrls['gene'] = cnrls.index
    fit_tab = pd.concat([fit, cnrls]).set_index('gene')
    final = fit_tab.merge(results, how='outer', left_index=True, right_index=True)
    col_order = col_order = ['num_samples', 'fitness_mean', 'fitness_std', 'ci', 'zscore', 'pval', 'padj']
    final_cols = ['locus', 'num_barcodes', 'library', 'barcode', 'sstart', 'sseqid'] + [f'{d}_{s}' for d in days for s
                                                                                        in col_order]
    return final[final_cols]

@st.cache
def analyze_experiment(exp_df, dnaid, experiment, good_samples, control_file, to_filter, outdir):


    sdf, edf = generate_DE_dataset(exp_df, good_samples, filter_below=to_filter)
    # Run DESeq2
    #st.write("Running DESeq2")
    fitness = get_fintess_results(outdir, dnaid, experiment, sdf, edf)
    # Calculate z-scores
    #st.write('Calculating z-scores')
    results = calculte_comparisons(fitness, exp_df, control_file)
    # Summarize results
    #st.write('Summarizing')
    final = final_fitness_table(fitness, exp_df, control_file, results)
    return final