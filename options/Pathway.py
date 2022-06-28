# import escher
# from escher import Builder
# import streamlit as st
# import streamlit.components.v1 as components
# import pandas as pd
# import numpy as np
# import plotly.express as px
#
# def app(datadir):
#     st.write('Pathway test')
#     clrs = px.colors.qualitative.Plotly
#     fdf = pd.read_csv("/Users/ansintsova/git_repos/fly_rnaseq/data/results_for_app/Test_Control_unfiltered_results_ann.kegg.csv")
#     fdf['logpval'] = -10 * np.log10(fdf.padj)
#     fdf['logpval'] = fdf['logpval'].fillna(0.1)
#     fdf['hits'] = fdf.padj < 0.05
#     pathways = fdf.KEGG_Pathway.dropna().unique()
#     p = st.selectbox('Choose Pathway', pathways)
#     #d = st.selectbox('Choose day', days)
#     subDf = fdf[(fdf.KEGG_Pathway == p)].sort_values('log2FoldChange')
#     fig = px.scatter(subDf, x='Symbol', y='log2FoldChange', size='logpval',
#                      category_orders={"Symbol": subDf.Symbol.values},
#                      color_discrete_map={
#                          True: clrs[1],
#                          False: clrs[0]},
#                      labels={'Name': '', 'z-score': 'Gene LFC'},
#                      height=700,
#                      color='hits')
#     fig.add_hline(y=0, line_width=2, line_dash="dash", line_color="grey",
#                   )
#     fig.update_layout(autosize=True, font=dict(size=18), paper_bgcolor='rgba(0,0,0,0)',
#                       )
#     fig.update_traces(marker=dict(
#                                   line=dict(width=2,
#                                             color='DarkSlateGrey')),
#                       selector=dict(mode='markers'))
#     st.plotly_chart(fig, use_container_width=True)
#
#     # gene_info = pd.read_csv('./data/SL1344_test/Path/S6_RNA-seq_aerobic_to_anaerobic.csv', index_col=0).to_dict()['lfc']
#     # #gene_info = {"xapA": 7.387150687875563, "xapB": 7.5306018918703765, "speF": 8.340221502890024}
#     # builder = Builder(
#     #     map_name= 'iJO1366.Central metabolism',
#     #     #model_name='e_coli_core',
#     #
#     # )
#     # builder.gene_info = gene_info
#     # builder.save_html('./data/example_map.html')
#     # HtmlFile = open("./data/example_map.html", 'r', encoding='utf-8')
#     # source_code = HtmlFile.read()
#     # components.html(source_code, height=800)