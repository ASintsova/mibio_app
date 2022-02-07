import streamlit as st
import tarfile
from pathlib import Path

# Custom imports
from multipage import MultiPage
from pages import PCA, DiffAb, Summary
st.set_page_config(page_title="mBARq App", layout='wide')

DATADIR = Path('./data/test')
# Create an instance of the app
app = MultiPage()
# Title of the main page
st.title("NCCR Microbiomes ETHZ prototype")
app.add_page("PCA", PCA.app)
# fname = st.sidebar.file_uploader('Project Results', accept_multiple_files=False)
# if not fname:
#     st.stop()
# tar = tarfile.open(fileobj=fname, mode="r:gz")
# tar.extractall(DATADIR)
# project_sites = [s.name for s in DATADIR.glob('**/*') if not s.is_file()]
# project_name = project_sites[0]
#
# st.header(f'Project Name: {project_name}')
#
# results = {'Summary': ('Summary', Summary.app),
#            'PCA': ('PCA', PCA.app),
#            'DiffAb': ('Differential Expression/Abundance', DiffAb.app)}
#
# for page_name, page in results.items():
#     if page_name not in project_sites:
#         continue
#     app.add_page(page[0], page[1])

#app.add_page("Analyse featureCounts Data", featureCounts_analysis.app)
#app.add_page("Explore DE analysis results", explore_results.app)

# The main app
app.run()