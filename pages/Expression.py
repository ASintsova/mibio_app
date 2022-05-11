import escher
from escher import Builder
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
from itertools import cycle
import plotly.express as px
import yaml

def app(datadir):
    st.subheader('Gene Expression')
    with open(datadir / "pages.yaml") as fh:
        config = yaml.safe_load(fh)

    clrs = px.colors.qualitative.Plotly
    sampleData = pd.read_csv(datadir / "sampleData.csv", index_col=0)
    countData = pd.read_csv(list(datadir.glob("*tpms*.csv"))[0])
    annotation_cols =config['annotation']
    sampleID = config['sampleID'][0]
    gene_name = st.radio('Choose gene annotation', annotation_cols)

    with st.expander('Show Gene Expression'):
        df = countData.set_index(gene_name).copy()
        df = df.apply(lambda x: np.log2(x + 0.5) if np.issubdtype(x.dtype, np.number) else x)
        sampleDataAb = sampleData.reset_index()
        df = df.reset_index()
        c1, c2 = st.columns(2)
        compare_by = c1.selectbox('Compare by', sampleDataAb.columns)
        color_by = c2.selectbox('Color by',  list(sampleDataAb.columns))
        genes = st.multiselect("Choose gene(s) of interest", df[gene_name].unique())

        if not genes:
            st.stop()
        c3, c4 = st.columns(2)
        tpm_label = 'log2 (TPM)'
        for col, gene in zip(cycle([c3, c4]), genes):
            gene_df = df[df[gene_name] == gene]
            samples = list(sampleDataAb[sampleID].values)
            gene_df = gene_df[[gene_name] + samples]

            gene_df = (gene_df.melt(id_vars=[gene_name], value_name=tpm_label, var_name=sampleID)
                       .merge(sampleData, how='left', on=sampleID))
            gene_df = gene_df.sort_values(compare_by)
            fig = px.box(gene_df, title=gene, x=compare_by, y=tpm_label, color=color_by,
                           hover_data=[gene_name] + list(sampleData.columns))
            fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                              font=dict(size=16))
            fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='LightGrey')
            col.plotly_chart(fig, use_container_width=True)
