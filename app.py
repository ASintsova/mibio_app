import streamlit as st
import tarfile
from pathlib import Path
import yaml


# Custom imports
from multipage import MultiPage
from pages import PCA, DiffAb, Summary, Pathway, Assembly, Expression, Home
st.set_page_config(page_title="NCCR Microbiomes ETHZ", layout='wide')

DATADIR = Path('/Users/ansintsova/git_repos/tnseq_app/data/ath_rnaseq')

#DATADIR = Path('/Users/ansintsova/git_repos/tnseq_app/data/assembly_test')
# fname = st.sidebar.file_uploader('Project Results', accept_multiple_files=False)
# if not fname:
#     st.stop()
# tar = tarfile.open(fileobj=fname, mode="r:gz")
# tar.extractall(DATADIR)


with open(DATADIR/"pages.yaml") as fh:
    config = yaml.safe_load(fh)

st.title("NCCR Microbiomes ETHZ")


# Create an instance of the app
app = MultiPage()
project_sites = config['pages']

results = {
           'Summary': ('Summary', Summary.app, DATADIR),
           'PCA': ('PCA', PCA.app, DATADIR),
           'DiffAb': ('Differential Expression/Abundance', DiffAb.app, DATADIR),
           'Path': ('Pathway Analysis',  Pathway.app, DATADIR),
           'Expression': ('Expression', Expression.app, DATADIR),
           'Assembly': ('Assembly', Assembly.app, DATADIR)}

results = {'Home': ('Home', Home.app, DATADIR),
           'PCA': ('PCA', PCA.app, DATADIR),
           'DiffAb': ('Differential Expression/Abundance', DiffAb.app, DATADIR)}
for page_name, page in results.items():
    if page_name not in project_sites:
        continue
    app.add_page(page[0], page[1], page[2])

#app.add_page('Assembly', Assembly.app)
#app.add_page('Pathway Analysis', Pathway.app)
#app.add_page('Differential Expression', DiffAb.app)




#


#app.add_page("Analyse featureCounts Data", featureCounts_analysis.app)
#app.add_page("Explore DE analysis results", explore_results.app)

# The main app
app.run()