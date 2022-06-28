import streamlit as st
import yaml

def app(datadir):

    descriptions = {'EDA': ['EDA', 'Explore patterns and identify outliers'],
                    'DiffAb': ['Differential Gene Expression', 'Explore diffrential expresion results produced by DESeq2'],
                    'Expression': ['Explore expression levels for any gene of interest']}

    with open(datadir / "pages.yaml") as fh:
        config = yaml.safe_load(fh)

    st.header(config['projectName'])
    for page in config['pages']:
        if page in descriptions.keys():
            st.subheader(descriptions[page][0])
            st.write(descriptions[page][1])




