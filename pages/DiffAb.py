import streamlit as st
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.express as px
from itertools import cycle
import dash_bio as dashbio
import yaml

def load_data(design_file, count_file):
    sampleData = pd.read_csv(design_file, index_col=0)
    countData = pd.read_csv((count_file), index_col=0)
    return sampleData, countData


def app(datadir):
    #DATADIR = Path('/Users/ansintsova/git_repos/tnseq_app/data/ath_rnaseq')

    clrs = px.colors.qualitative.Plotly
    st.write('## Differentical Expression Results')

    fdf = pd.read_csv(list(datadir.glob("*unfiltered*results*kegg.csv"))[0])
    c1, c2, c3, c4 = st.columns(4)

    with open(datadir / "pages.yaml") as fh:
        config = yaml.safe_load(fh)
    annotation_cols = config['annotation']
    sampleID = config['sampleID'][0]
    pval = config['pval_col'][0]
    lfc = config['lfc_col'][0]
    gene_name = st.radio('Choose gene annotation', annotation_cols, key='ann')
    contrast = 'contrast'

    fdf['log10FDR'] = -10 * np.log10(fdf[pval])
    if contrast == 'All':
        df = fdf.copy()
    else:
        contrasts = fdf[contrast].unique()
        contrast_to_show = st.selectbox('Select a contrast', contrasts)
        df = fdf[fdf.contrast == contrast_to_show].copy()

    with st.expander('Show Volcano Plot'):
        c1, c2 = st.columns(2)
        fdr = c1.number_input('FDR cutoff', value=0.05)
        lfc_th = c2.number_input('Log FC cutoff (absolute)', value=1)
        df['hit'] = ((abs(df[lfc]) > lfc_th) & (df[pval] < fdr))
        fig = px.scatter(df, x=lfc, y='log10FDR', color='hit',
                         height=700,
                         color_discrete_map={
                             True: clrs[1],
                             False: clrs[0]},
                         hover_name=df.index, hover_data=[lfc, pval])
        fig.add_vline(x=lfc_th, line_width=2, line_dash="dash", line_color="grey")
        fig.add_vline(x=-lfc_th, line_width=2, line_dash="dash", line_color="grey")
        fig.add_hline(y=-10*np.log10(fdr), line_width=2, line_dash="dash", line_color="grey")
        fig.update_layout(autosize=True, font=dict(size=18), paper_bgcolor='rgba(0,0,0,0)',
                          )
        fig.update_traces(marker=dict(size=8,
                                      line=dict(width=1,
                                                color='DarkSlateGrey')),
                          selector=dict(mode='markers'))
        st.plotly_chart(fig, use_container_width=True)

    with st.expander('LFC rankings by Pathway'):
        contrast_col, lfc_col, fdr_col, lfc_lib_col = st.columns(4)
        fdr = fdr_col.number_input('FDR cutoff', value=0.05, key='kegg_pval')
        lfc_th = lfc_col.number_input('Log FC cutoff (absolute)', min_value=0.0, step=0.5, value=1.0,key='kegg_lfc')
        df['hit'] = ((abs(df[lfc]) > lfc_th) & (df[pval] < fdr))
        show_kegg = st.selectbox('Show KEGG Pathway', ['All'] + list(df.KEGG_Pathway.unique()))
        if show_kegg != 'All':
            df = df[df.KEGG_Pathway == show_kegg]
        df = df.sort_values(lfc).reset_index().reset_index().rename({'level_0': 'ranking'}, axis=1)
        fig = px.scatter(df, x='ranking', y=lfc, color='hit',
                         height=700,
                         color_discrete_map={
                             True: clrs[1],
                             False: clrs[0]},
                         hover_name=gene_name,
                         title=f"{contrast_to_show} - {show_kegg}",
                         hover_data={lfc: True,
                                     'log10FDR': False,
                                    'ranking': False,
                                     pval: True},
                         labels={"ranking": '', lfc: 'Log2 FC'}
                         )
        fig.add_hline(y=0, line_width=2, line_dash="dash", line_color="grey")
        fig.update_xaxes(showticklabels=False)
        fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                          font=dict(size=18))
        fig.update_traces(marker=dict(size=14,
                                      line=dict(width=2,
                                                color='DarkSlateGrey')),
                          selector=dict(mode='markers'))
        st.plotly_chart(fig, use_container_width=True)


    # compContrasts = st.multiselect('Select contrasts to compare', contrasts)
    # if not compContrasts:
    #     st.stop()
    # c1, c2 = st.columns(2)
    # filters = {}
    # for col, contrast in zip(cycle([c1, c2]), compContrasts):
    #    col.write(contrast)
    #    l = col.number_input('LFC cutoff', value=2, key=f'{contrast}_lfc')
    #    f = col.number_input('FDR cutoff', value=0.05, key=f'{contrast}_fdr')
    #    filters[contrast] = (l, f)
    #
    # compDfs = []
    # for key, value in filters.items():
    #     if value[0] > 0:
    #         df = fdf[(fdf.contrast == key) & (fdf.log2FoldChange > value[0])&(fdf.padj < value[1])]
    #     else:
    #         df = fdf[(fdf.contrast == key) & (fdf.log2FoldChange < value[0]) & (fdf.padj < value[1])]
    #     compDfs.append(df)
    #
    # vennDf = pd.concat(compDfs)
    # vennDf = vennDf.reset_index().pivot(index='index', columns='contrast', values='log2FoldChange').dropna()
    #
    # if vennDf.empty:
    #     st.write('No genes matching the filtering criteria found')
    #     st.stop()
    # st.write(vennDf.shape)
    # columns = list(vennDf.columns)
    # rows = list(vennDf.index)
    # # fix error if only one gene or 1 sample
    #
    #
    # cluster = 'all' if len(columns) > 1 else 'row'
    # clustergram = dashbio.Clustergram(
    #     data=vennDf.loc[rows].values,
    #     row_labels=rows,
    #     column_labels=columns,
    #     color_threshold={
    #         'row': 250,
    #         'col': 700
    #     },
    #     height=800,
    #     width=700,
    #     center_values=False,
    #     cluster=cluster
    # )
    #
    # if vennDf.shape[0] < 100:
    #     st.plotly_chart(clustergram, use_container_width=True)
    # else:
    #     st.write('Too many genes to display. Download table?')
