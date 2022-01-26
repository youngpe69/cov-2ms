from re import M
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import seaborn as sns
import json as jn

# page configuration
st.set_page_config(
    page_title="CoVs-2 Mutations",
    layout='wide',
    menu_items={
        'Get Help': 'https://twitter.com/youngvintage69',
        'Report a bug': 'https://github.com/youngpe69/cov-2ms/issues',
        'About': 'Simple streamlit app for exploring SARS-CoV-2 mutations'
    }
)

# session initialisation
if 'mutSelection' not in st.session_state:
    st.session_state['mutSelection'] = set()
if 'mutByDate' not in st.session_state:
    st.session_state['mutByDate'] = pd.DataFrame()

# initialisation
DATA_URL = ('https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv')
DATA_FILEPATH = ('../corona/mut_data/cog_metadata.csv')
MUTATION_FILEPATH = ('./mutations.json')
genes = {'S', 'M', 'N', 'E', 'Orf3a', 'Orf7a', 'Orf8'}

variants = set()
mutations = dict()
currentGene = 'S'

# main area and sidebar 
st.title('CoVs-2 Mutations')
st.subheader('Plot Mutation Frequency')
st.text('Plot mutation frequency across all the UK SARS-CoV-2 sequences held in MRC-CLIMB using the filters below.')

# mutation json
def load_mutationTable(variants, mutations):
    with open(MUTATION_FILEPATH) as f:
        mutation_table = jn.load(f)
        print(type(mutation_table))
        for i in mutation_table:
            variants.add(i['Variant'])
            mutations[i['Variant']] = i
        print(variants)
    return mutation_table

# COG-Metadata loader
@st.cache
def load_data():
    cogDf = pd.read_csv(DATA_URL)
    cogDf = cogDf.rename(columns={"sample_date": "date"})
    cogDf["date"] = pd.to_datetime(cogDf["date"])
    cogDf = cogDf.sort_values(by='date')
    return cogDf

# search function
def search_data(muts, cogDf):
    print(muts)
    mutByDate = pd.DataFrame()
    seriesList = list()
    for mut in muts:
        idSeries = cogDf['mutations'].str.contains(mut, na=False)

        m = cogDf[idSeries]
        countByDate = m.groupby('date').count()
        countByDate = countByDate.rename(columns={"sequence_name": mut})
        series = countByDate[mut]
        seriesList.append(series)
    mutByDate = pd.concat(seriesList, axis=1)
    return mutByDate

# load the mutation table
mutation_table = load_mutationTable(variants, mutations)

# columns
col1, col2, col3 = st.columns(3)

#  variant
with col1:
    selected_variant = st.selectbox('Select a variant', variants)

# gene
with col2:
    option = st.selectbox('Select a gene', (genes))
    currentGene = option

# mutation
with col3:
    if currentGene in mutations[selected_variant]:
        mutSelection = st.multiselect('Select mutations', mutations[selected_variant][currentGene])
        st.session_state['mutSelection'].update(mutSelection)

    
# total selected mutations
if st.session_state['mutSelection']:
    mutTotalSelection = st.multiselect('Selected mutations', st.session_state['mutSelection'], default=st.session_state['mutSelection'])
    st.session_state['mutSelection'].clear()
    st.session_state['mutSelection'].update(mutTotalSelection)


# run a search
if st.button('Search') :
    cogDf = load_data()
    st.session_state['mutByDate'] = search_data(st.session_state['mutSelection'], cogDf)

st.line_chart(data=st.session_state['mutByDate'])
st.caption("Data source: MRC-CLIMB via https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv")
st.caption("Contact| twitter: @youngvintage69")





