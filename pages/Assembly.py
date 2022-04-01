import streamlit as st
import pandas as pd
from sklearn.decomposition import PCA
import numpy as np
import plotly.express as px


def app(datadir):
    st.write(datadir)
    st.write('# Assembly Stats')
    st.write('## Statistics for contigs >= 500 bp')
    # Input standart report file
    df1 = pd.read_table("/Users/ansintsova/git_repos/tnseq_app/data/assembly_test/unicycler.report.tsv")
    df2 = pd.read_table("/Users/ansintsova/git_repos/tnseq_app/data/assembly_test/report.tsv")
    df3 = pd.read_table("/Users/ansintsova/git_repos/tnseq_app/data/assembly_test/contigs_report.tsv")
    df = df1.merge(df2, on='Assembly').merge(df3, on='Assembly').set_index("Assembly").T
    stats_to_show = ['# contigs', 'Largest contig', 'Total length', 'GC (%)', 'N50', 'L50']
    stats_to_show = [s for s in stats_to_show if s in df.columns]

    df_all = df[stats_to_show].copy().astype(float).reset_index().rename({'index':'sampleID'}, axis=1)
    c1, c2 = st.columns(2)
    for i, stat in enumerate(stats_to_show):
        f = px.bar(df_all, x='sampleID', y=stat, width=400,
                                         color='sampleID', hover_data=df_all)
        f.update_layout(autosize=True, font=dict(size=18), paper_bgcolor='rgba(0,0,0,0)',
                          )
        f.update_xaxes(showticklabels=False)
        if df_all.sampleID.nunique() < 4:
            if i%2 == 0:
                c2.plotly_chart(f,use_container_width=True)
            else:
                c1.plotly_chart(f, use_container_width=True)
        else:
            st.plotly_chart(f, use_container_width=True)

    if_ref = ['# misassemblies', 'Misassembled contigs length', '# local misassemblies', 'Genome fraction (%)',
              "# N's per 100 kbp", "# mismatches per 100 kbp"]
    if_ref = [s for s in if_ref if s in df.columns]
    if len(if_ref) > 0:
        st.write("## Reference-based stats")
        c3, c4 = st.columns(2)
        df_ref = df[if_ref].copy().astype(float).reset_index().rename({'index': 'sampleID'}, axis=1)
        for i, stat in enumerate(if_ref):
            f = px.bar(df_ref, x='sampleID', y=stat, width=400,
                       color='sampleID', hover_data=df_ref)
            f.update_layout(autosize=True, font=dict(size=18), paper_bgcolor='rgba(0,0,0,0)',
                            )
            f.update_xaxes(showticklabels=False)
            if df_ref.sampleID.nunique() < 4:
                if i % 2 == 0:
                    c3.plotly_chart(f, use_container_width=True)
                else:
                    c4.plotly_chart(f, use_container_width=True)
            else:
                st.plotly_chart(f, use_container_width=True)
    st.write(df.T)

