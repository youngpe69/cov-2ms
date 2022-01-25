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

# session initialisation
if 'mutSelection' not in st.session_state:
    st.session_state['mutSelection'] = set()

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
st.text('Plot mutation frequency across all the UK SARS-CoV-2 sequences held in MRC-CLIMB\nusing the filters on the left.')
st.sidebar.title("Filters")

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
    cogDf = pd.read_csv(DATA_FILEPATH)
    cogDf = cogDf.rename(columns={"sample_date": "date"})
    cogDf["date"] = pd.to_datetime(cogDf["date"])
    cogDf = cogDf.sort_values(by='date')
    return cogDf

# search function
def search_data(muts, cogDf):
    print(muts)
    mutByDate = pd.DataFrame()
    for mut in muts:
        idSeries = cogDf['mutations'].str.contains(mut, na=False)

        m = cogDf[idSeries]
        countByDate = m.groupby('date').count()
        countByDate = countByDate.rename(columns={"sequence_name": mut})
        series = countByDate[mut]
        mutByDate[mut] = series
    return mutByDate

# data display
def searchDisplay(mutByDate, fromDate, toDate):
    plotTypePY(mutByDate, fromDate, toDate)

def plotTypePY(mutByDate, fromDate, toDate):
    fig2, ax2 = plt.subplots(num=None, figsize=(25, 20), facecolor="w", edgecolor="k")

    fig2.text(
            0.51,
            0.05,
            f"Date: {str(dt.date.today())} | Mutation Frequency Patterns | Data source: MRC-CLIMB via https://cog-uk.s3.climb.ac.uk/phylogenetics/latest/cog_metadata.csv | @youngvintage69",
            size=16,
            va="bottom",
            ha="center",
        )

    sns.lineplot(data=mutByDate, palette="tab10", linewidth=2.5)

    plt.legend(loc=2, prop={'size': 20})
    locator = mdates.AutoDateLocator()
    formatter = mdates.ConciseDateFormatter(locator)
    #ax2.set_ylim(0,100)
    ax2.set_xlim([fromDate, toDate])
    ax2.tick_params(axis='both', which='major', labelsize=20)
    ax2.tick_params(axis='both', which='minor', labelsize=10) 
    ax2.xaxis.set_major_locator(locator)
    ax2.xaxis.set_major_formatter(formatter)
    ax2.xaxis.set_minor_locator(ticker.MultipleLocator(base=1.0))
    ax2.set_xlabel("Sample Date")
    ax2.xaxis.label.set_size(22)
    ax2.set_ylabel("Count")
    ax2.yaxis.label.set_size(22)
    ax2.grid(True, which="major", linewidth=0.25)
    ax2.grid(True, which="minor", linewidth=0.1)
    ax2.set_axisbelow(True)
    plt.show()

    st.pyplot(fig2)

def plotTypeST(mutByDate):
    st.line_chart(data=mutByDate)


# load the mutation table
mutation_table = load_mutationTable(variants, mutations)

# sidebar controls
#
#  variant
selected_variant = st.sidebar.selectbox('Select a variant', variants)

# gene
option = st.sidebar.selectbox(
    'Select a gene',
     (genes))
currentGene = option

# mutation
mutSelection = st.sidebar.multiselect('Select mutations', mutations[selected_variant][currentGene])

# total selected mutations
totalSelection = st.session_state['mutSelection']
totalSelection.update(mutSelection)
mutTotalSelection = st.sidebar.multiselect('Selected mutations', st.session_state['mutSelection'], default=st.session_state['mutSelection'])
totalSelection.clear()
totalSelection.update(mutTotalSelection)

# date controls
fromDate = st.sidebar.date_input(
     "From",
     dt.date(2020, 4, 1))
toDate = st.sidebar.date_input("To", dt.date(2022, 1, 22))

# run a search
if st.sidebar.button('Search') :
    cogDf = load_data()
    mutByDate = search_data(st.session_state['mutSelection'], cogDf)
    searchDisplay(mutByDate, fromDate, toDate)





