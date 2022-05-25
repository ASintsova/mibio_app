import streamlit as st
import tarfile
from pathlib import Path
import yaml

# Custom imports
from multipage import MultiPage
from pages import EDA, DiffAb, Summary, Pathway, Assembly, Expression, Home
st.set_page_config(page_title="NCCR Microbiomes ETHZ", layout='wide')

#DATADIR = Path('/Users/ansintsova/git_repos/tnseq_app/data/ath_rnaseq')
#DATADIR = Path('/Users/ansintsova/git_repos/tnseq_app/data/results_for_app')

DATADIR = Path(".")
fname = st.sidebar.file_uploader('Project Results', accept_multiple_files=False)
if not fname:
    st.stop()
tar = tarfile.open(fileobj=fname, mode="r:gz")
tar.extractall()



with open("pages.yaml") as fh:
    config = yaml.safe_load(fh)
st.stop()
st.title("NCCR Microbiomes ETHZ")

# Create an instance of the app
app = MultiPage()
project_sites = config['pages']

# results = {
#            'Summary': ('Summary', Summary.app, DATADIR),
#            'Path': ('Pathway Analysis',  Pathway.app, DATADIR),
#            'Assembly': ('Assembly', Assembly.app, DATADIR)}

available_pages = {'Home': ('Home', Home.app, DATADIR),
           'EDA': ('Exploratory Data Analysis', EDA.app, DATADIR),
           'DiffAb': ('Differential Expression/Abundance', DiffAb.app, DATADIR),
           'Expression': ('Gene Expression', Expression.app, DATADIR)}

for page_name, page in available_pages.items():
    if page_name not in project_sites:
        continue
    app.add_page(page[0], page[1], page[2])

# The main app
app.run()