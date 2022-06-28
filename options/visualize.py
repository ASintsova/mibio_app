#import plotnine as p9
import pandas as pd
import numpy as np


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


def graph_wits(exp_df, controls, corr_df, phenotype, what='all', file_path=''):
    cntrl_counts = (controls.merge(exp_df[['barcode', 'cnt', 'sampleID', 'day', 'mouse']], how='left', on='barcode')
                        .drop(['DN'], axis=1))
    cntrl_counts = cntrl_counts.merge(corr_df, on=['phenotype', 'sampleID'])
    cntrl_counts = cntrl_counts[cntrl_counts.phenotype == phenotype]
    cntrl_counts['R_lab'] = cntrl_counts.R.apply(lambda x: f'R = {round(x, 2)}')
    if what == 'all':
        to_graph = cntrl_counts.copy()
    elif what == 'inoculum':
        inoc_ids = [f for f in cntrl_counts.sampleID.values if 'inoculum' in f]
        to_graph = cntrl_counts[cntrl_counts.sampleID.isin(inoc_ids)].copy()
    else:
        to_graph = cntrl_counts[cntrl_counts.mouse.isin(what)].copy()

    x = to_graph.day.nunique()
    y = to_graph.mouse.nunique()

    p9.options.figure_size = (x*3, y*3)
    g = (p9.ggplot(to_graph, p9.aes(x='conc', y='cnt'))
      + p9.geom_point()
      + p9.geom_smooth(method="lm")
      + p9.theme_classic()
      + p9.theme(text=p9.element_text(size=14),
                 axis_text_x=p9.element_text(rotation=90, hjust=1))
      + p9.ylab("Count")
      + p9.xlab("Expected Abundance")
      + p9.ggtitle(phenotype)
      + p9.scale_y_log10()
      + p9.scale_x_log10()
      + p9.geom_text(p9.aes(label='R_lab', x=0.0001, y=.1))
      + p9.facet_grid('mouse~day'))
    if file_path:
        p9.save_as_pdf_pages([g + p9.theme(figure_size=(x*3, y*3))], file_path)
    return g


def view_inoculum_counts(fd):
    g = (p9.ggplot(fd, p9.aes(x='cnt', fill='mouse'))
          + p9.geom_histogram(bins=100)
          + p9.theme_classic()
          + p9.theme(text = p9.element_text(size=14),
                     axis_text_x=p9.element_text(rotation=90, hjust=1))
          + p9.xlab("Count")
          + p9.facet_grid('mouse~day'))
    return g


def filter_inoculum(exp_df, filter_below=0):
    filt_df = (exp_df.copy()
               .drop(['ShortName', 'locus_tag'], axis=1)
               .drop_duplicates()
               .pivot(index='barcode', columns='sampleID', values='cnt'))
    columns_to_filter = [f for f in filt_df.columns if 'inoculum' in f]
    filt_df = filt_df[(filt_df[columns_to_filter] >= filter_below).all(1)]
    return filt_df


#def filter_samples(exp_df, good_samples):
#    return exp_df.copy()[exp_df.sampleID.isin(good_samples)]


def filter_all_exps(df, to_filter=1000):
    filt_dfs = []
    for i,g in df.groupby(['dnaid', 'experiment']):
        fildf = filter_inoculum(g, to_filter).reset_index().rename({'index':'barcode'}, axis=1)
        fildf = fildf.melt(id_vars=['barcode'], value_name='cnt').assign(dnaid=i[0], experiment=i[1])
        filt_dfs.append(fildf)
    fdf = pd.concat(filt_dfs)
    tdf = fdf.merge(df[['barcode', 'ShortName', 'locus_tag', 'dnaid','experiment', 'sampleID', 'mouse', 'day', 'organ']],
                  on=['dnaid', 'experiment', 'sampleID', 'barcode'], how='left')
    return tdf


def view_barcodes(df, gene, to_filter=0):
    count_df = filter_all_exps(df, to_filter)

    gene_df = count_df[count_df.ShortName == gene]
    nbc = gene_df.barcode.nunique()
    if nbc == 0:
        return f"{gene} not found"
    inoculum = gene_df[(gene_df.day == 'd0') & (gene_df.mouse == 'inoculum')]
    to_graph = gene_df[gene_df.day != 'd0']
    if nbc / 4 < 1:
        xdim = 4 * nbc
        print(xdim)
        ydim = 5
    else:
        xdim = 16
        ydim = 5 * nbc / 1.5

    p9.options.figure_size = (xdim, ydim)

    g = (p9.ggplot(to_graph, p9.aes(x='day', y='cnt', color='mouse', shape='experiment', group='mouse'))
         + p9.geom_point(size=6)
         + p9.geom_line()
         + p9.theme_classic()
         + p9.ylab("Count")
         + p9.xlab("Day")
         + p9.theme(text=p9.element_text(size=18))

         + p9.scale_y_log10()
         + p9.facet_wrap("~barcode")
         + p9.geom_hline(inoculum, p9.aes(yintercept='cnt', color='dnaid'), linetype="dashed", size=1)
         )
    return g