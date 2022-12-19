import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



st.title('Methylation Portal')
BLCAdata = './TCGADatasets/BLCA/GeneTableBLCAPVals.csv'
BRCAdata = './TCGADatasets/BRCA/GeneTableBRCAPVals.csv'
CESCdata = './TCGADatasets/CESC/GeneTableCESCPVals.csv'
COADdata = './TCGADatasets/COAD/GeneTableCOADPVals.csv'
ESCAdata = './TCGADatasets/ESCA/GeneTableESCAPVals.csv'
GBMdata = './TCGADatasets/GBM/GeneTableGBMPVals.csv'
HNSCdata = './TCGADatasets/HNSC/GeneTableHNSCPVals.csv'
KIRCdata = './TCGADatasets/KIRC/GeneTableKIRCPVals.csv'
KIRPdata = './TCGADatasets/KIRP/GeneTableKIRPPVals.csv'
LGGdata = './TCGADatasets/LGG/GeneTableLGGPVals.csv'
LIHCdata = './TCGADatasets/LIHC/GeneTableLIHCPVals.csv'
LUADdata = './TCGADatasets/LUAD/GeneTableLUADPVals.csv'
LUSCdata = './TCGADatasets/LUSC/GeneTableLUSCPVals.csv'
OVdata = './TCGADatasets/OV/GeneTableOVPVals.csv'
PAADdata = './TCGADatasets/PAAD/GeneTablePAADPVals.csv'
PCPGdata = './TCGADatasets/PCPG/GeneTablePCPGPVals.csv'
PRADdata = './TCGADatasets/PRAD/GeneTablePRADPVals.csv'
READdata = './TCGADatasets/READ/GeneTableREADPVals.csv'
SARCdata = './TCGADatasets/SARC/GeneTableSARCPVals.csv'
SKCMdata = './TCGADatasets/SKCM/GeneTableSKCMPVals.csv'
STADdata = './TCGADatasets/STAD/GeneTableSTADPVals.csv'
TGCTdata = './TCGADatasets/TGCT/GeneTableTGCTPVals.csv'
THCAdata = './TCGADatasets/THCA/GeneTableTHCAPVals.csv'
THYMdata = './TCGADatasets/THYM/GeneTableTHYMPVals.csv'
UCECdata = './TCGADatasets/UCEC/GeneTableUCECPVals.csv'

all_data = [BLCAdata,BRCAdata,CESCdata,COADdata,ESCAdata,GBMdata,HNSCdata,KIRCdata,KIRPdata,LGGdata,LIHCdata,LUADdata,LUSCdata,OVdata,PAADdata,PCPGdata,PRADdata,READdata,SARCdata,SKCMdata,STADdata,TGCTdata,THCAdata,THYMdata,UCECdata]
#BLCA,BRCA,CESC,COAD,ESCA,GBM,HNSC,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC
#all_datasets = [BLCA,BRCA,CESC,COAD,ESCA,GBM,HNSC,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC]

@st.cache
def load_data(dataset):
    data = pd.read_csv(dataset)
    return data

def load_gene_data(dataset):
    gene_data = dataset.loc[dataset['Gene'] == text_input]
    st.write(gene_data)
    return gene_data

# Create a text element and let the reader know the data is loading.
data_load_state = st.text('Loading data...')
# Load 10,000 rows of data into the dataframe.
for i,data in enumerate(all_data):
    all_data[i] = load_data(data)
# Notify the reader that the data was successfully loaded.
data_load_state.text("Done! (using st.cache)")

st.subheader('Enter Gene')
text_input = st.text_input("Enter Gene: ")
if text_input:
    st.write("You entered: ", text_input)

st.subheader('Raw data')
st.write(all_data[0])

for data in all_data:
    load_gene_data(data)

fig,ax = plt.subplots()
# ax.hist = plt.hist(BLCA_gene_data[['MedianBvalTumor']],bins = np.arange(0,1.2,0.2))
st.pyplot(fig)