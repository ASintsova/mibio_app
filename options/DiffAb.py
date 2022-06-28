import streamlit as st
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.express as px
from itertools import cycle
import dash_bio as dashbio
import yaml
import requests
from time import sleep



def load_data(design_file, count_file):
    sampleData = pd.read_csv(design_file, index_col=0)
    countData = pd.read_csv((count_file), index_col=0)
    return sampleData, countData

@st.cache
def convert_df(df):
 # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv().encode('utf-8')


def app(datadir):
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
    contrast_col = 'contrast'
    contrasts = fdf[contrast_col].unique()
    contrast_to_show = st.selectbox('Select a contrast', ['All'] + list(contrasts))
    fdf['log10FDR'] = -10 * np.log10(fdf[pval])
    if contrast_to_show == 'All':
        df = fdf.copy()
    else:
        df = fdf[fdf[contrast_col] == contrast_to_show].copy()
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
                         hover_name=df[gene_name], hover_data=[lfc, pval])
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
        lfc_col, fdr_col = st.columns(2)
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

    st.subheader("Protein-protein interactions")
    with st.expander('Tabular Results'):
        vennDf = fdf.copy()
        compContrasts = st.multiselect('Select contrasts to display', contrasts)
        c1, c2 = st.columns(2)
        filters = {}
        for col, con in zip(cycle([c1, c2]), compContrasts):
            col.write(con)
            l = col.number_input('LFC cutoff', value=-1.0, step=0.5, key=f'{con}_lfc')
            f = col.number_input('FDR cutoff', value=0.05, step=0.01, key=f'{con}_fdr')
            filters[con] = (l, f)
        comp_genes = []

        for key, value in filters.items():
            if value[0] > 0:
                genes = set(vennDf[(vennDf[contrast_col] == key) & (vennDf[lfc] > value[0])
                                   & (vennDf[pval] < value[1])][gene_name].values)
            else:
                genes = set(vennDf[(vennDf[contrast_col] == key) & (vennDf[lfc] < value[0])
                                   & (vennDf[pval] < value[1])][gene_name].values)
            comp_genes.append(genes)
        if not comp_genes:
            st.stop()
        intersect_genes = set.intersection(*comp_genes)
        vennDf = vennDf[vennDf[gene_name].isin(intersect_genes)].copy()
        vennDf = vennDf[[gene_name, lfc, pval, contrast_col]].drop_duplicates()
        #vennDf2 = vennDf.pivot(index='Name', columns='contrast_col', values=['LFC', 'fdr'])
        st.write(vennDf.shape)
        string_col, download_col = st.columns(2)
        string_api_url = "https://version-11-5.string-db.org/api"
        output_format = 'tsv-no-header'
        method = 'get_link'
        my_genes = set(vennDf[gene_name].values)
        request_url = "/".join([string_api_url, output_format, method])
        string_col.markdown("### STRING Interaction Network")
        species = string_col.number_input("NCBI species taxid", value=7227, help='Drosophila melanogaster: 7227')

        params = {
            "identifiers": "\r".join(my_genes),  # your protein
            "species": species,  # species NCBI identifier
            "network_flavor": "confidence",  # show confidence links
            "caller_identity": "mibio"  # your app name
        }
    #
        if string_col.button('Get STRING network'):
            results = requests.post(request_url, data=params)
            network_url = results.text.strip()
            st.markdown(f"[Link to STRING network]({network_url})")
            sleep(1)
        download_col.markdown("### Download results as csv")
        fname_default = config['projectName'].replace(' ', '_')
        fname = download_col.text_input("File name", value=fname_default)
        fname = fname+".csv"
        download_col.download_button("Download data as csv file", convert_df(vennDf), file_name=fname)

    # st.subheader("Gene Selector")
    # with st.expander('Gene Selector'):
    #     genes = st.multiselect("Choose gene(s) of interest", fdf[identifier].unique())
    #     if not genes:
    #         st.stop()
    #     c3, c4 = st.columns(2)
    #
    #     for col, gene in zip(cycle([c3, c4]), genes):
    #         gene_df = fdf[fdf[identifier] == gene].copy()
    #         gene_df = gene_df[[identifier, 'library', 'contrast_col', 'LFC', 'fdr']].drop_duplicates()
    #         gene_df = gene_df.sort_values('contrast_col')
    #         fig = px.strip(gene_df, title=gene, x="contrast_col", y='LFC', color='library',
    #                        hover_data=[identifier,"library", "contrast_col", "fdr"], stripmode='overlay')
    #         fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
    #                           font=dict(size=16))
    #         fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='LightGrey')
    #         fig.add_hline(y=0, line_width=2, line_dash="dash", line_color="grey")
    #         col.plotly_chart(fig, use_container_width=True)
    #
    #


