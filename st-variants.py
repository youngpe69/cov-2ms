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
import networkx as nx

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
if 'initialise_plot' not in st.session_state:
    st.session_state['initialise_plot'] = True
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

#
# Functions
# 
#
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

#@st.cache
def load_graph():
    G = nx.read_graphml(GRAPH_FILEPATH)
    return G

# search function
def search_data(muts, cogDf):
    mutByDate = pd.DataFrame()
    seriesList = list()
    for mut in muts:
        idSeries = cogDf['mutations'].str.contains(mut, na=False)

        m = cogDf[idSeries]
        countByDate = m.groupby('date').count()
        countByDate = countByDate.rename(columns={"sequence_name": mut})
        series = countByDate[mut]
        seriesList.append(series)
        #mutByDate[mut] = series
    mutByDate = pd.concat(seriesList, axis=1)
    return mutByDate
     

# load the mutation table
mutation_table = load_mutationTable(variants, mutations)

# get initial data to display in plot
if st.session_state['initialise_plot']:
    cogDf = load_data()
    #st.session_state['mutByDate'] = search_data(st.session_state['mutSelection'], cogDf)
    st.session_state['initialise_plot'] = False

#
# Actual page controls
# 
# main area and sidebar 
st.title('CoVs-2 Mutations')
st.subheader('Variant')

# columns for spacing
col1a, col1b, col1c = st.columns(3)
with col1a:
    selected_variant = st.selectbox('Select a variant', variants)

st.header(" ")
st.subheader('Mutation Frequency Plot')
st.text("Select a gene and add mutations to plot their observed frequency over time")

    

# columns for mutation selection
col2a, col2b = st.columns(2)
    
# gene selector
with col2a:
    option = st.selectbox('Select a gene', (genes))
    currentGene = option

# mutation selector
with col2b:
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


# Mutation relationships graph
st.header(" ")
st.subheader("Explore Mutation Relationships")

# columns for graph filter
col3a, col3b = st.columns(2)

with col3a:
    selected_threshold = st.slider(label='Select Threshold', min_value=0.1, max_value=1.0, value=0.7)

with col3b:
    select_edge_labels = st.checkbox('Show Weights', value=True)

# Now create the equivalent Node and Edge lists
G = load_graph()

# Create net for weight threshold
eligible_edges = [(from_node,to_node,edge_attributes) for from_node,to_node,edge_attributes in G.edges(data=True) if edge_attributes['weight'] > selected_threshold]
gmut = nx.Graph()
gmut.add_edges_from(eligible_edges)

# Create objects for graph UI
nodes = [Node(id=i, label=str(i), size=300) for i in gmut.nodes]
edges = [Edge(source=i, target=j, type="STRAIGHT", label=f"{gmut[i][j]['weight']:.1f}") for (i,j) in gmut.edges
        if i in gmut.nodes and j in gmut.nodes]


config = Config(width=1800, 
                height=900, 
                directed=False,
                nodeHighlightBehavior=True, 
                highlightColor="#F7A7A6", # or "blue"
                collapsible=False,
                node={'labelProperty':'label'},
                link={'labelProperty': 'label', 'renderLabel': select_edge_labels},
                staticGraphWithDragAndDrop=True,
                maxZoom=1,
                minZoom=1,
                initialZoom = 1,
                disableLinkForce = True
                # **kwargs e.g. node_size=1000 or node_color="blue"
                ) 

return_value = agraph(nodes=nodes, 
                      edges=edges, 
                      config=config)

st.caption("Contact| twitter: @youngvintage69")





