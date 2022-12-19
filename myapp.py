import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



st.title('Methylation Portal')

BRCA = '../Datasets/GeneTableBRCAPVals.csv'

@st.cache
def load_data():
    data = pd.read_csv(BRCA)
    return data

@st.cache
def load_gene_data(dataset):
    gene_data = dataset.loc[dataset['Gene'] == text_input]
    return gene_data

# Create a text element and let the reader know the data is loading.
data_load_state = st.text('Loading data...')
# Load 10,000 rows of data into the dataframe.
BRCA = load_data()
# Notify the reader that the data was successfully loaded.
data_load_state.text("Done! (using st.cache)")

st.subheader('Enter Gene')
text_input = st.text_input("Enter Gene: ")
if text_input:
    st.write("You entered: ", text_input)

st.subheader('Raw data')
st.write(BRCA)

BRCA_gene_data = load_gene_data(BRCA)
st.write(BRCA_gene_data)

fig,ax = plt.subplots()
ax.hist = plt.hist(BRCA_gene_data[['MedianBvalTumor']],bins = np.arange(0,1.2,0.2))
st.pyplot(fig)