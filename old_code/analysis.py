import streamlit as st

import numpy as np
import pandas as pd
from pathlib import Path

#import plotnine as p9
from scipy import stats
from scipy.stats import norm
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels as sm
import plotnine as p9
#import skmisc
#from matplotlib.patches import Ellipse
#from sklearn.decomposition import PCA



import subprocess
st.title('TNSeq Data Analysis')


# # These need to be input by the user
results_dir = "../data/processed/results"
controls = pd.read_table('../data/metadata/controls.txt', header=None,
                        names = ['DN', 'barcode', 'phenotype', 'conc'])
meta_dir = "../data/metadata"
# # Loading the data
#


@st.cache
def read_count_files(sample, exp, results_dir=results_dir):
    df = pd.read_table(Path(results_dir)/f'{sample}/{sample}_counts_{exp}.txt', sep=" ").assign(exp=exp)
    df = (df.reset_index().rename({'index':'barcode'}, axis=1)
          .melt(id_vars=['barcode','Position', 'Element', 'Strand', 'Feature', 'ShortName', 'exp']))
    df['proportion'] = df['value']/ df.groupby('variable')['value'].transform('sum')
    expansion = df['variable'].str.split('-', expand=True)
    df['mouse'], df['day'],df['organ'] = expansion[0], expansion[1], expansion[2]
    df = df.rename({'variable':'sample', 'value':'cnts'}, axis=1)
    return df

@st.cache
def load_sample(sample, meta_dir = meta_dir, results_dir = results_dir):
    meta = (pd.read_table(Path(meta_dir)/f'{sample}_metadata.txt', header=None,
                        names = ['DN', 'lib', 'exp', 'DN2', 'sample', 'day', 'organ'])
            .drop(['DN', 'DN2'], axis=1))
    exps = meta.exp.unique()
    dfs = [read_count_files(sample, exp) for exp in exps]
    fdf = pd.concat(dfs).assign(dnaid=sample)
    fdf['sample'] = fdf['sample'].str.replace('-', '_')
    #fdf = fdf.merge(meta[['lib', 'exp']], how='left', left_on='exp', right_on='exp')
    return fdf

@st.cache
def load_samples(samples, meta_dir=meta_dir, results_dir=results_dir):
    dfs = [load_sample(sample, meta_dir, results_dir) for sample in samples]
    return pd.concat(dfs)




def subset_experiment(df, query_string):
    '''
    example query string : '(exp=="TV5490A") & (dnaid == "dnaid2023")'
    '''
    return df.query(query_string)


def calculate_2dist_zscore(u1, s1, u2, s2):
    return (u1 - u2) / np.sqrt((s1 ** 2) + (s2 ** 2))


@st.cache
def get_fitness_results(fdf, query,  experiment):
    exp = subset_experiment(fdf, query)
    days = list(exp['day'].unique())
    days.remove('d0')
    all_results = pd.concat([pd.read_table(f"../data/tmp_{experiment}_{day}.txt", sep=' ').assign(day=day) for day in days])
    return all_results

@st.cache
def calculate_comparisons(fitness, df, controls, cntrl_type='wt'):
    genes = set(df.Feature.values)
    genes.remove('-')
    other_barcodes = set(df[df.Feature == '-'].barcode.values)
    control_barcodes = set(controls[controls.phenotype == cntrl_type].barcode.values)
    days = fitness.day.unique()
    all_comps = []
    for day in days:
        control_fits = fitness[fitness.day == day].loc[fitness.index.intersection(control_barcodes)]
        control_mu = control_fits.log2FoldChange.mean()
        control_sigma = np.sqrt(control_fits.lfcSE.pow(2).sum()) / control_fits.shape[0]
        gene_comps = {}
        for gene in genes:
            gene_barcodes = set(df[df.Feature == gene].barcode.values)
            gene_fits = fitness[fitness.day == day].loc[fitness.index.intersection(gene_barcodes)]
            if gene_fits.shape[0] > 0:
                gene_mu = gene_fits.log2FoldChange.mean()
                gene_sigma = np.sqrt(gene_fits.lfcSE.pow(2).sum()) / gene_fits.shape[0]
                zscore = calculate_2dist_zscore(gene_mu, gene_sigma, control_mu, control_sigma)
                ci = 2 ** gene_mu / 2 ** control_mu
                num_bc = gene_fits.shape[0]
                meanExp = gene_fits.baseMean.mean()
                std = np.std(gene_fits.baseMean)

                gene_comps[gene] = [zscore, num_bc, meanExp, std, ci]
        for barcode in other_barcodes:
            other_fit = fitness[fitness.day == day].loc[fitness.index.intersection([barcode])]
            if not other_fit.empty:
                zscore = calculate_2dist_zscore(other_fit.log2FoldChange.values[0],
                                                other_fit.lfcSE.values[0], control_mu, control_sigma)
                ci = 2 ** other_fit.log2FoldChange.values[0] / 2 ** control_mu
                gene_comps[barcode] = [zscore, 1, other_fit.baseMean.values[0], 0, ci]

        comp_df = pd.DataFrame(gene_comps,
                               index=[day + '_zscore', day + '_num_bc', day + "_meanExpr", day + "_std", day + "_ci"]).T

        all_comps.append(comp_df)
    return pd.concat(all_comps, axis=1)


def comp_stats(comp):
    pvalues = 2 * norm.cdf(-np.abs(comp), 0, 1)
    p_adjust = sm.stats.multitest.multipletests(pvalues, method='fdr_bh')[1]
    s = pd.DataFrame(comp)
    s[comp.name + '_pval'] = pvalues
    s[comp.name + '_padj'] = p_adjust
    s = s.rename({comp.name: comp.name + "_zscore"})
    return s






@st.cache
def subset_experiment(df, query_string):
    '''
    example query string : '(exp=="TV5490A") & (dnaid == "dnaid2023")'
    '''
    return df.query(query_string)

@st.cache
def good_mice(df, controls, cutoff):
    control_cnts = controls.merge(df, left_on='barcode', right_on = 'barcode').drop(['DN','Position', 'Element',
                                                                                 'Strand', 'Feature', 'ShortName'], axis=1)

    wt_cnts = control_cnts.copy()[control_cnts.phenotype == 'wt']
    wt_cnts['logConc'] = np.log10(wt_cnts['conc'])
    wt_cnts['logCnts'] = np.log10(wt_cnts['cnts'])
    log_corr = calculate_correlation(wt_cnts, 'sample', 'logConc', 'logCnts').set_index('sample')
    log_corr.columns = ['R']
    return log_corr[log_corr.R > cutoff].index

@st.cache
def generate_DE_dataset(df, samples_to_keep, to_filter=0):  # Assumes already subset to 1 experiment, 1 dnaid
    sample_data = df[['sample', 'mouse', 'day', 'organ', 'dnaid']].set_index('sample').drop_duplicates()
    sample_data = sample_data.loc[sample_data.index.intersection(samples_to_keep)]
    expr_data = df.pivot(index='barcode', columns='sample', values='cnts')
    expr_data = expr_data[
        (expr_data['inoculum_d0_inoculum'] >= to_filter) & (expr_data['unenriched_inoculum_d0_inoculum'] >= to_filter)]
    expr_data = expr_data[list(sample_data.index)].reset_index()
    return sample_data, expr_data




def analyze_experiment(fdf, query, to_keep, controls, to_filter=0, cntrl_type = 'wt'):
    exp1 = subset_experiment(fdf, query)

    sdf, edf = generate_DE_dataset(exp1, to_keep, to_filter)
    sdf.to_csv(f"../data/{experiment}_sdf.test")
    edf.set_index('barcode').to_csv(f"../data/{experiment}_edf.test")
    run_command(['Rscript', 'DEseq.R', f"../data/{experiment}_sdf.test", f"../data/{experiment}_edf.test"])
    fitness = get_fitness_results(fdf, query, experiment)
    barcode_info = exp1[['barcode', 'Feature', 'ShortName']].drop_duplicates().set_index('barcode')
    fitness_annot = fitness.merge(barcode_info, how='left', left_index=True, right_index = True)
    comp_to_wt = calculate_comparisons(fitness, exp1, controls, cntrl_type)
    final_list = [comp_stats(comp_to_wt[c]) for c in comp_to_wt.columns if 'zscore' in c]
    return fitness_annot, comp_to_wt, pd.concat(final_list, axis=1)



def counts_overtime(fdf, exp, gene):
    df = fdf[(fdf.ShortName == gene) & (fdf.exp == exp)]
    nbc = df.barcode.nunique()
    if nbc == 0:
        return f"{gene} not found"
    print(nbc / 4)
    inoculum = df[(df.day == 'd0') & (df.mouse == 'inoculum')]
    df = df[df.day != 'd0']
    if nbc / 4 < 1:

        xdim = 4 * nbc
        print(xdim)
        ydim = 5
    else:
        xdim = 16
        ydim = 5 * nbc / 4
    p9.options.figure_size = (xdim, ydim)
    g = (p9.ggplot(df, p9.aes(x='day', y='cnts', color='mouse', group='mouse'))
         + p9.geom_point(size=6)
         + p9.geom_line()
         + p9.theme_classic()
         + p9.ylab("Counts")
         + p9.xlab("Day")
         + p9.scale_y_log10()
         + p9.geom_hline(inoculum, p9.aes(yintercept='cnts', color='dnaid'), linetype="dashed", size=1)
         + p9.facet_wrap("~barcode")

         )

    return g



# Load all the count data for specified dnaids.

fdf = load_samples(['dnaid2023', 'dnaid2024']) # Do not mutate this
gene_info = fdf.copy()[['Feature', 'ShortName']].drop_duplicates().set_index('Feature')
st.write("This is where you would input your data...")
st.write("Loading Count Data...")
st.write(f"Loaded data for {', '.join(fdf.dnaid.unique())} ")
st.write(f"The following experiments are available: {', '.join(fdf.exp.unique())}")

st.header("Step 1: Choose Experiment")

col1, col2 = st.beta_columns(2)
with col1:
    experiment = st.selectbox("Experiment: ", tuple(['Choose an experiment'] + list(fdf.exp.unique()))) # Example 'TV5490A'
with col2:
    dnaid = st.selectbox("dnaid: ", tuple(['Choose a dnaid'] + list(fdf.dnaid.unique())))


if (experiment == 'Choose an experiment') | (dnaid == 'Choose a dnaid'):
    st.write('No experiment selected')
    st.stop()
else:
    st.write(f'You selected {experiment}, {dnaid}')

#####################################################
query = f'(exp=="{experiment}") & (dnaid == "{dnaid}")'
exp_df = subset_experiment(fdf, query)
#####################################################

st.header("Step 2: Let's Look at the controls")

# WT Standard Curve:

# todo make this a function
control_cnts = controls.merge(exp_df, left_on='barcode', right_on='barcode').drop(['DN', 'Position', 'Element',
                                                                                    'Strand', 'Feature', 'ShortName'],
                                                                                   axis=1)

wt_cnts = control_cnts[control_cnts.phenotype == 'wt']
wt_cnts['logConc'] = np.log10(wt_cnts['conc'])
wt_cnts['logCnts'] = np.log10(wt_cnts['cnts'])


def calculate_correlation(df, groupby, v1, v2):
    corr_df = df.groupby(groupby)[[v1, v2]].corr()
    corr_df = corr_df.reset_index()
    corr_df = corr_df[corr_df['level_1'] == v1].drop(['level_1', v1], axis=1).set_index('sample')
    corr_df.columns = ['R']
    return corr_df

wt_corr_df = calculate_correlation(wt_cnts, 'sample', 'logConc', 'logCnts')




# Overall Correlation:
# todo split by days
overall_corr_df = exp_df.pivot(index='barcode', columns='sample', values='proportion').corr()



control_option = st.selectbox("What would you like to see? ", ('Options', 'Standard Curve', 'Overall Correlation',  'Done with Controls'))
if control_option == 'Options':
    st.write('No Control Selected')
    st.stop()
elif control_option == 'Standard Curve':
    ctrl_type = st.selectbox("Which controls?", tuple(controls.phenotype.unique()))
    ctrl_cnts = control_cnts[control_cnts.phenotype == ctrl_type]
    ctrl_cnts['logConc'] = np.log10(wt_cnts['conc'])
    ctrl_cnts['logCnts'] = np.log10(wt_cnts['cnts'])

    if st.button('Show R'):
        st.write(calculate_correlation(wt_cnts, 'sample', 'logConc', 'logCnts'))
    if st.button('Show graphs'):
    # Plotting the correlations
        p9.options.figure_size = (10, 40)
        g = (p9.ggplot(ctrl_cnts, p9.aes(x='conc', y='cnts'))
             + p9.geom_point()
             + p9.geom_smooth(method="lm")
             + p9.theme_classic()
             + p9.ylab("Count")
             + p9.xlab("Conc")
             + p9.scale_y_log10()
             + p9.scale_x_log10()
             + p9.facet_grid('mouse~day'))
        st.pyplot(p9.ggplot.draw(g))

elif control_option == 'Overall Correlation':
    if st.button('Show R'):
        mean_corr = overall_corr_df.mean()
        mean_corr.columns = ['R']
        st.write(mean_corr)
    if st.button('Show graphs'):
        sns.set_style("white")
        sns.set_context("notebook", font_scale=2.0)
        mmap = sns.color_palette("Blues", as_cmap=True)
        g = sns.clustermap(overall_corr_df, linewidths=0.5, linecolor='black', figsize=(15, 15), vmin=0.6, vmax=1, cmap=mmap,
                           cbar_kws={'label': 'Correlation Coefficient'})
        st.pyplot(g)

else:
    st.write("Moving on...")

st.header("Step 3: Analyze Data")

st.write("Choose Mice passing controls")
pass_controls = st.selectbox("Mice passing controls: ", ('Choose Mice', 'All', 'Standard Curve cutoff', 'Overall Correlation cutoff'))

if pass_controls == 'Choose Mice':
    st.write("No Mice selected")
    st.stop()

elif pass_controls == 'All':
    to_keep = list(exp_df['sample'].unique())


elif pass_controls == 'Standard Curve cutoff':
    cutoff  = st.number_input('Cutoff to use: ')
    to_keep = list(wt_corr_df[wt_corr_df['R'] > cutoff].index)


elif pass_controls == 'Overall Correlation cutoff':
    cutoff  = st.number_input('Cutoff to use: ')
    to_keep = list(overall_corr_df[overall_corr_df.mean() > cutoff].index)

else:
    to_keep = []
    st.write("OhOh")
    st.stop()


to_filter = st.number_input("Filter counts below: ")

analysis = st.selectbox("Analyze", ('Options', 'Analyze', 'Load previously analyzed'))


if analysis == 'Analyze':
    fitness_annot, comp_to_wt, zstats = analyze_experiment(fdf, query, to_keep, controls, to_filter=to_filter)

elif analysis == 'Load previously analyzed':
    fitness = get_fitness_results(fdf, query, experiment)
    barcode_info = exp_df[['barcode', 'Feature', 'ShortName']].drop_duplicates().set_index('barcode')
    fitness_annot = fitness.merge(barcode_info, how='left', left_index=True, right_index=True)
    comp_to_wt = calculate_comparisons(fitness, exp_df, controls, 'wt')
    final_list = [comp_stats(comp_to_wt[c]) for c in comp_to_wt.columns if 'zscore' in c]
    zstats = pd.concat(final_list, axis=1)
else:
    st.write('Ready?')
    st.stop()



st.header("Step 4: Visualize Results")

results = (comp_to_wt.merge(zstats[[f'{day}_zscore_padj' for day in exp_df.day.unique() if day != 'd0']], left_index=True, right_index=True)
           .merge(gene_info, how='left', left_index=True, right_index=True))


for day in exp_df.day.unique():
    if day == 'd0':
        continue
    cols_to_show = [f'{day}_zscore', f'{day}_ci']
    st.write(day)
    sig_results = results[[f'{day}_zscore', f'{day}_ci', f'{day}_zscore_padj', 'ShortName']]
    sig_results.columns = ['z-score', 'CI', 'padj', 'ShortName']

    st.dataframe(sig_results[sig_results.padj < 0.05])



goi =  st.text_input('Gene of Interest')
if not goi:
    st.stop()


g = counts_overtime(fdf, experiment, goi)
st.pyplot(p9.ggplot.draw(g))

# Volcano Plot



#
# @st.cache
# def vis_fitness(fitness, controls, gene, day):
#     fit = fitness[fitness.day == day]
#     cf = (controls[controls.phenotype == 'wt'].set_index('barcode').merge(fit, how='left', left_index=True,
#                                                                           right_index=True)
#           .drop(['DN'], axis=1))
#     cf = cf.drop_duplicates()
#     gf = fit[fit.ShortName == gene]
#
#     sns.set_style("white")
#     sns.set_context("notebook", font_scale=1.5)
#     g = plt.figure(figsize=(4, 3))
#     cf.log2FoldChange.hist(bins=25, label='Control Barcodes')
#     gf.log2FoldChange.hist(bins=15, label=f'{gene} Barcodes')
#     plt.legend()
#     plt.xlabel('log2FoldChange compared to inoculum')
#     return cf, gf, g


for day in exp_df.day.unique():
    if day == 'd0':
        continue
    st.subheader(f'Day: {day}')

    #cf, gf, g2 = vis_fitness(fitness_annot, controls, goi, day)
    #st.write(gf)
    #st.pyplot(g2)

    #st.write(f'Fitness: {gf.log2FoldChange.mean()}')
    z_score   = results[results["ShortName"] == goi][[f"{day}_zscore"]].values[0][0]

    padj = results[results["ShortName"] == goi][[f"{day}_zscore_padj"]].values[0][0]
    CI = results[results["ShortName"] == goi][[f"{day}_ci"]].values[0][0]
    st.write( f"Z-score: {z_score}, {padj}")
    st.write(f"CI: {CI}")



#
#
#




### One Experiment ####





# Output files will be for now in '../data/tmp_{day}.txt change that to smth more meaningful


#statsmodels.stats.multitest.multipletests(pvals, alpha=0.05, method='hs', is_sorted=False, returnsorted=False)


#

# # Visualization
#
# original_fitness, original_comps, original_results = analyze_experiment(fdf, 'TV5490A', controls, 'wt', dnaid='dnaid2023')
#
#
#
#
