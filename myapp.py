import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



st.title('Methylation Portal')
BLCA = './TCGADatasets/BLCA/GeneTableBLCAPVals.csv'
BRCA = './TCGADatasets/BRCA/GeneTableBRCAPVals.csv'
CESC = './TCGADatasets/CESC/GeneTableCESCPVals.csv'
COAD = './TCGADatasets/COAD/GeneTableCOADPVals.csv'
ESCA = './TCGADatasets/ESCA/GeneTableESCAPVals.csv'
GBM = './TCGADatasets/GBM/GeneTableGBMPVals.csv'
HNSC = './TCGADatasets/HNSC/GeneTableHNSCPVals.csv'
KIRC = './TCGADatasets/KIRC/GeneTableKIRCPVals.csv'
KIRP = './TCGADatasets/KIRP/GeneTableKIRPPVals.csv'
LGG = './TCGADatasets/LGG/GeneTableLGGPVals.csv'
LIHC = './TCGADatasets/LIHC/GeneTableLIHCPVals.csv'
LUAD = './TCGADatasets/LUAD/GeneTableLUADPVals.csv'
LUSC = './TCGADatasets/LUSC/GeneTableLUSCPVals.csv'
OV = './TCGADatasets/OV/GeneTableOVPVals.csv'
PAAD = './TCGADatasets/PAAD/GeneTablePAADPVals.csv'
PCPG = './TCGADatasets/PCPG/GeneTablePCPGPVals.csv'
PRAD = './TCGADatasets/PRAD/GeneTablePRADPVals.csv'
READ = './TCGADatasets/READ/GeneTableREADPVals.csv'
SARC = './TCGADatasets/SARC/GeneTableSARCPVals.csv'
SKCM = './TCGADatasets/SKCM/GeneTableSKCMPVals.csv'
STAD = './TCGADatasets/STAD/GeneTableSTADPVals.csv'
TCGT = './TCGADatasets/TCGT/GeneTableTCGTPVals.csv'
THCA = './TCGADatasets/THCA/GeneTableTHCAPVals.csv'
THYM = './TCGADatasets/THYM/GeneTableTHYMPVals.csv'
UCEC = './TCGADatasets/UCEC/GeneTableUCECPVals.csv'

all_datasets = [BLCA,BRCA,CESC,COAD,ESCA,GBM,HNSC,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TCGT,THCA,THYM,UCEC]

@st.cache
def load_data(dataset):
    data = pd.read_csv(dataset)
    return data

@st.cache
def load_gene_data(dataset):
    gene_data = dataset.loc[dataset['Gene'] == text_input]
    return gene_data

# Create a text element and let the reader know the data is loading.
data_load_state = st.text('Loading data...')
# Load 10,000 rows of data into the dataframe.
for dataset in all_datasets:
    dataset = load_data(dataset)
# Notify the reader that the data was successfully loaded.
data_load_state.text("Done! (using st.cache)")

st.subheader('Enter Gene')
text_input = st.text_input("Enter Gene: ")
if text_input:
    st.write("You entered: ", text_input)

st.subheader('Raw data')
st.write(BLCA)

BLCA_gene_data = load_gene_data(BLCA)
st.write(BLCA_gene_data)

fig,ax = plt.subplots()
ax.hist = plt.hist(BRCA_gene_data[['MedianBvalTumor']],bins = np.arange(0,1.2,0.2))
st.pyplot(fig)