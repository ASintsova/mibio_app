import streamlit as st

import pandas as pd
#from diffexpr.py_deseq import py_DESeq2

def load_data(design_file, count_file):
    sampleData = pd.read_csv(design_file, index_col=0)
    countData = pd.read_csv((count_file), index_col=0)
    return sampleData, countData


def app():
    st.write('# RNAseq/RBSeq app prototype')
    st.write('## Loading the data')
    st.write("Please upload experimental design table and count file. No spaces in experimental desing column names,"
             " index of the experimental desing should be same as the columns of the count file")
    design_file = st.file_uploader('Experimental Design File',  accept_multiple_files=False)
    count_file = st.file_uploader('Count File',  accept_multiple_files=False)

    if design_file is None or count_file is None:
        st.stop()

    samples, counts = load_data(design_file, count_file)
    st.write(samples.head())

    #st.write('## Subsetting the data')
    st.write('## Validating the data')
    countData = counts[[c for c in samples.index if c in counts.columns]].copy()
    sampleData = samples.copy()
    st.write(countData.sample(5))
    st.write('## Experimental Design')

    design = []

    design_options = list(sampleData.columns)
    not_include = []
    batch = st.radio(
         "Choose Batch Factor",
         ["None"] + design_options)

    if batch != 'None':
        design.append(batch)
        not_include.append(batch)

    #
    factor1 = st.radio(
         "Choose First Factor",
         [c for c in design_options if c not in not_include])

    not_include.append(factor1)

    factor2 = st.radio(
         "Choose Second Factor",
         ["None"]+[c for c in design_options if c not in not_include])

    if factor2 != "None":
        sampleData['group'] = sampleData[factor1].astype(str) + "_" + sampleData[factor2].astype(str)
        design.append('group')

    else:
        design.append(factor1)
    st.write(sampleData.head())

    design_str = ' + '.join(design)
    #
    #
    manual_design = st.text_input('Enter design manually', value="")
    if manual_design:
        design_str = manual_design
        design = [d.strip() for d in design_str.split("+")]

    st.write(f"Design: {design_str}")

    compRef = st.selectbox(
         'Compare: ',
         sampleData[design[-1]].unique())
    compTreats = st.multiselect('to ',
         [c for c in sampleData[design[-1]].unique() if c != compRef])

    comparisons = [(c, compRef) for c in compTreats]
    st.write(comparisons)



    countData = countData.reset_index()
    #st.write(countData.head())
    dds = py_DESeq2(count_matrix=countData,
                    design_matrix=sampleData,
                    design_formula='~ group',
                    gene_column='Geneid')  # <- telling DESeq2 this should be the gene ID column

    dds.run_deseq()

    dds.get_deseq_result(contrast=['group', 'Ileum_PBS', 'Ileum_Lactulose'])
    res = dds.deseq_result
    st.write(res.head())

    st.write('## Transforming the data')

# @st.cache
# def load_files():
#     results_dir = "data/counts"
#     files = [f for f in Path(results_dir).iterdir()]
#     controls = pd.read_table('data/controls.txt', header=None, names=['DN', 'barcode', 'phenotype', 'conc'])
#     return pd.concat([pd.read_csv(f, index_col=0) for f in files]), controls
#
#
# @st.cache
# def subset_experiment(df,  exp, dnaid, col1='experiment', col2='dnaid'):
#     '''
#     example query string : '(exp=="TV5490A") & (dnaid == "dnaid2023")'
#     '''
#     query_string = f'({col1} == "{exp}") & ({col2} == "{dnaid}")'
#     return df.copy().query(query_string)
#
# @st.cache
# def filter_inoculum(exp_df, filter_below=0):
#     filt_df = (exp_df.copy()
#                .drop(['ShortName', 'locus_tag'], axis=1)
#                .drop_duplicates()
#                .pivot(index='barcode', columns='sampleID', values='cnt'))
#     columns_to_filter = [f for f in filt_df.columns if 'inoculum' in f]
#     filt_df = filt_df[(filt_df[columns_to_filter] >= filter_below).all(1)].reset_index()
#     fd = filt_df.melt(id_vars=['barcode'], var_name='sampleID', value_name='cnt').fillna(0)
#     fd = fd.merge(exp_df[['barcode', 'sampleID', 'mouse', 'day']], on=['barcode', 'sampleID'])
#     return fd
#
#
# @st.cache
# def filter_all_exps(df, to_filter=1000):
#     filt_dfs = []
#     for i, g in df.groupby(['dnaid', 'experiment']):
#         fildf = filter_inoculum(g, to_filter).reset_index().rename({'index': 'barcode'}, axis=1)
#         fildf = fildf.melt(id_vars=['barcode'], value_name='cnt').assign(dnaid=i[0], experiment=i[1])
#         filt_dfs.append(fildf)
#     fdf = pd.concat(filt_dfs)
#     tdf = fdf.merge(
#         df[['barcode', 'ShortName', 'locus_tag', 'dnaid', 'experiment', 'sampleID', 'mouse', 'day', 'organ']],
#         on=['dnaid', 'experiment', 'sampleID', 'barcode'], how='left')
#     return tdf
#
#
# def get_table_download_link(df):
#     """Generates a link allowing the data in a given panda dataframe to be downloaded
#     in:  dataframe
#     out: href string
#     """
#     csv = df.to_csv(index=False)
#     b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
#     href = f'<a href="data:file/csv;base64,{b64}"> Download results </a>'
#     return href
#
# def main():
#     st.title("TNSeq Data Exploration and Analysis")
#     df, controls = load_files()
#
#     # Give some info
#     st.write(f"Loaded data for {', '.join(sorted(df.dnaid.unique()))} ")
#     to_display = ['barcode', 'ShortName', 'locus_tag', 'cnt', 'sampleID', 'dnaid', 'experiment']
#     st.write("Sample of the loaded data: ")
#     st.dataframe(df[to_display].sample(10))
#
#     st.write("# Step 1: Choose Experiment")
#     col1, col2 = st.beta_columns(2)
#     with col1:
#         dnaid = st.selectbox("dnaid: ", tuple(['Choose a dnaid'] + list(df.dnaid.unique())))
#     with col2:
#         experiment = st.selectbox("Experiment: ",
#                                   tuple(['Choose an experiment'] + list(df[df.dnaid == dnaid].experiment.unique())))  # Example 'TV5490A'
#
#     if (experiment == 'Choose an experiment') | (dnaid == 'Choose a dnaid'):
#         st.write('No experiment selected')
#         st.stop()
#     else:
#         st.write(f'You selected {experiment}, {dnaid}')
#
#     exp_df = subset_experiment(df, experiment, dnaid)
#
#     st.write("# Data Exploration")
#
#     st.markdown("## Looking at WITS")
#     corr_df, good_samples = calculate_correlation(exp_df, "data/controls.txt")
#
#     phenotype = st.selectbox("Phenotype: ",  tuple(['Choose phenotype'] + list(controls.phenotype.unique())))
#     st.write("You can graph all samples (All), just the inoculum (Inoculum), or choose one or more mouse tags")
#     what_to_plot = st.multiselect('Samples', (["All", "Inoculum"] + list(exp_df.mouse.unique())))
#     if not what_to_plot:
#         st.stop()
#     elif 'All' in what_to_plot:
#         samples = 'all'
#     elif 'Inoculum' in what_to_plot:
#         samples = 'inoculum'
#     else:
#         samples = what_to_plot
#
#     g = graph_wits(exp_df, controls, corr_df, phenotype, what=samples, file_path='')
#     st.pyplot(p9.ggplot.draw(g))
#
#     #st.markdown("## Look at the count data distribution")
#     #st.write("Filter barcodes with low counts in the inoculum?")
#     #to_filter = st.number_input("Filter barcodes with counts below:", min_value=0, value=0)
#
#
#     st.write("# Analyze data")
#
#     # Choose correlation cutoff
#
#     st.markdown("## Choose Mice")
#
#     cutoff = st.slider('Cutoff', 0.0, 1.0, 0.9, 0.01)
#     wits_corr_df, good_mice = calculate_correlation(exp_df, "data/controls.txt", cutoff=cutoff)
#
#
#     mice = (exp_df[(exp_df.sampleID).isin(good_mice)][['day', 'mouse']]
#             .drop_duplicates()
#             .sort_values('day'))
#
#     bad_mice = (exp_df[~exp_df.sampleID.isin(good_mice)][['day', 'mouse']]
#             .drop_duplicates()
#             .sort_values('day').groupby('day').agg({'mouse': ['unique']}).reset_index())
#     if bad_mice.empty:
#         mice_sum = mice.groupby('day').agg({'mouse': ['count', 'unique']}).reset_index()
#         mice_sum.columns = ['day', 'Number of samples', 'Passing Samples']
#         mice_sum['Passing Samples'] = mice_sum['Passing Samples'].apply(lambda x: ", ".join(x))
#         mice_sum['Failed Samples'] = 'N/A'
#     else:
#         mice_sum = mice.groupby('day').agg({'mouse':['count', 'unique']}).reset_index().merge(bad_mice, on='day')
#         mice_sum.columns = ['day', 'Number of samples', 'Passing Samples', 'Failed Samples']
#         mice_sum['Passing Samples'] = mice_sum['Passing Samples'].apply(lambda x: ", ".join(x))
#         mice_sum['Failed Samples'] = mice_sum['Failed Samples'].apply(lambda x: ", ".join(x))
#     st.table(mice_sum)
#
#     st.write(" ## Filter barcodes with low counts in the inoculum?")
#     filter_below = st.number_input('Filter barcodes with counts below:', min_value=0, value=1000,)
#     filt_df = filter_inoculum(exp_df, filter_below=filter_below)
#     filt_d0 = filt_df[filt_df.day == 'd0']
#     st.table((filt_d0.groupby('mouse').barcode.count()
#               .reset_index()
#               .rename({'mouse': 'sample', 'barcode': 'Number of barcodes'}, axis=1))
#              .set_index('sample'))
#     g2 = view_inoculum_counts(filt_d0)
#     st.pyplot(p9.ggplot.draw(g2))
#
#     st.write('## Analysis Parameters:')
#     st.markdown(f"### Parameters:\n")
#     st.markdown(f"dnaid: {dnaid}\n")
#     st.markdown(f"experiment: {experiment}\n")
#     st.markdown(f"samples: {', '.join(good_mice)}\n")
#     st.markdown(f"filtering below: {filter_below}")
#
#
#     results = analysis.analyze_experiment(exp_df, dnaid, experiment,
#                                     good_mice, "data/controls.txt", filter_below, "data/analysis")
#     results = results.reset_index().rename({'index':'barcode'}, axis=1)
#
#     st.markdown(get_table_download_link(results), unsafe_allow_html=True)
#     days = sorted(exp_df.day.unique())
#     days.remove('d0')
#     st.markdown("## Show hits (padj < 0.05)")
#     hitday = st.selectbox('Choose Day', tuple(days))
#
#     hits = results[results[f'{hitday}_padj'] < 0.05]
#     hits = hits[['gene', 'locus', 'num_barcodes', 'library'] + [c for c in hits.columns if hitday in c]]
#     st.write(hits)
#
#
#     st.write("## Show Specific Genes")
#
#     gene = st.text_input("Gene Name or Locus tag")
#     if gene in results.gene.values:
#         st.write(results[results.gene == gene])
#     elif gene in results.locus.values:
#         st.write(results[results.locus == gene])
#     else:
#         st.write(f'{gene} not found')
#         st.stop()
#     st.write("### View Barcode count data")
#     data_choice = st.radio('Show all data?', ('This experiment', 'All data'))
#     filter_choice = st.radio("Filtered?", ('All Barcodes', 'Filtered Barcodes'))
#
#     if data_choice == 'This experiment':
#         to_graph = exp_df
#     elif data_choice == 'All data':
#         to_graph = df
#     else:
#         st.stop()
#
#     if filter_choice == 'All Barcodes':
#         graph_filter = 0
#     elif filter_choice == 'Filtered Barcodes':
#         graph_filter = filter_below
#     else:
#         st.stop()
#
#     g3 = view_barcodes(to_graph, gene, to_filter=graph_filter)
#     st.pyplot(p9.ggplot.draw(g3))
#     st.stop()
#
#
#
#         #st.write('Choose columns to explore')
#
#         #subset = st.multiselect('', tuple(results.columns))
#         #res = results[subset]
#         #st.write(res.sample(5))
#
#         # st.stop()
#         #
#         # st.write(len(set(results.gene.values)))
#         # st.dataframe(results.sample(25))
#         #     # fit, final = process_results(exp_df, controls, good_mice, filter_below=filter_below, cntrl_type='wt')
#         #     # final['experiment'] = experiment
#         #     # st.dataframe(final.sample(10))
#         #
#

# if __name__ == "__main__":
#     main()