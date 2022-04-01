import streamlit as st
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.express as px
from itertools import cycle
import dash_bio as dashbio


def load_data(design_file, count_file):
    sampleData = pd.read_csv(design_file, index_col=0)
    countData = pd.read_csv((count_file), index_col=0)
    return sampleData, countData


def app(datadir):
    DATADIR = Path('/Users/ansintsova/git_repos/tnseq_app/data/ath_rnaseq')

    clrs = px.colors.qualitative.Plotly
    st.write('## Differentical Expression Results')

    fdf = pd.read_csv(DATADIR/'diffab.kegg.csv')
    c1, c2, c3, c4 = st.columns(4)
    pval = c1.selectbox('Select p-value column', [None] + list(fdf.columns), key='pval')
    lfc = c2.selectbox('Select LFC column', [None] + list(fdf.columns), key='lfc')
    ann = c3.selectbox('Select annotation column', [None] + list(fdf.columns), key='ann')
    contrast = c4.selectbox('Select contrast column', [None, 'All'] + list(fdf.columns), key='contr')
    if not pval or not lfc or not ann or not contrast:
        st.stop()

    fdf['log10FDR'] = -10*np.log10(fdf[pval])
    if contrast == 'All':
        df = fdf.copy()
    else:
        contrasts = fdf[contrast].unique()
        contrast_to_show = st.selectbox('Select a contrast', contrasts)
        df = fdf[fdf.contrast == contrast_to_show].copy()

    fdr = c2.number_input('FDR cutoff', value=0.05)
    lfc_th = c3.number_input('Log FC cutoff (absolute)', value=1)
    df['hit'] = ((abs(df[lfc]) > lfc_th) & (df[pval] < fdr))
    fig = px.scatter(df, x=lfc, y='log10FDR', color='hit',
                     height=700,
                     color_discrete_map={
                         True: clrs[1],
                         False: clrs[0]},
                     hover_name=df.index)
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
    #
    compContrasts = st.multiselect('Select contrasts to compare', contrasts)
    if not compContrasts:
        st.stop()
    c1, c2 = st.columns(2)
    filters = {}
    for col, contrast in zip(cycle([c1, c2]), compContrasts):
       col.write(contrast)
       l = col.number_input('LFC cutoff', value=2, key=f'{contrast}_lfc')
       f = col.number_input('FDR cutoff', value=0.05, key=f'{contrast}_fdr')
       filters[contrast] = (l, f)

    compDfs = []
    for key, value in filters.items():
        if value[0] > 0:
            df = fdf[(fdf.contrast == key) & (fdf.log2FoldChange > value[0])&(fdf.padj < value[1])]
        else:
            df = fdf[(fdf.contrast == key) & (fdf.log2FoldChange < value[0]) & (fdf.padj < value[1])]
        compDfs.append(df)

    vennDf = pd.concat(compDfs)
    vennDf = vennDf.reset_index().pivot(index='index', columns='contrast', values='log2FoldChange').dropna()

    if vennDf.empty:
        st.write('No genes matching the filtering criteria found')
        st.stop()
    st.write(vennDf.shape)
    columns = list(vennDf.columns)
    rows = list(vennDf.index)
    # fix error if only one gene or 1 sample


    cluster = 'all' if len(columns) > 1 else 'row'
    clustergram = dashbio.Clustergram(
        data=vennDf.loc[rows].values,
        row_labels=rows,
        column_labels=columns,
        color_threshold={
            'row': 250,
            'col': 700
        },
        height=800,
        width=700,
        center_values=False,
        cluster=cluster
    )

    if vennDf.shape[0] < 100:
        st.plotly_chart(clustergram, use_container_width=True)
    else:
        st.write('Too many genes to display. Download table?')
