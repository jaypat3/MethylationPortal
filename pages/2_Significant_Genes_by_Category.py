import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from math import log

st.set_page_config(page_title="Significant Genes by Category")

st.title('Significant Genes by Category')

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
TGCT = './TCGADatasets/TGCT/GeneTableTGCTPVals.csv'
THCA = './TCGADatasets/THCA/GeneTableTHCAPVals.csv'
THYM = './TCGADatasets/THYM/GeneTableTHYMPVals.csv'
UCEC = './TCGADatasets/UCEC/GeneTableUCECPVals.csv'
ProbeLink = './TCGADatasets/ProbeMapping/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy.csv'
GeneLink = './TCGADatasets/ProbeMapping/probeMap_hugo_gencode_good_hg19_V24lift37_probemap.csv'

all_data = [BLCA,BRCA,CESC,COAD,ESCA,GBM,HNSC,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC]
#BLCA,BRCA,CESC,COAD,ESCA,GBM,HNSC,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC
#all_datasets = [BLCA,BRCA,CESC,COAD,ESCA,GBM,HNSC,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC]
#OV REMOVED FOR NOW
gene_data_array = []

@st.cache
def load_data(dataset):
    data = pd.read_csv(dataset)
    return data

def switch(data):
    if data == "BLCA":
        return all_data[0]
    elif data == "BRCA":
        return all_data[1] 
    elif data == "CESC":
        return all_data[2]
    elif data == "COAD":
        return all_data[3]
    elif data == "ESCA":
        return all_data[4] 
    elif data == "GBM":
        return all_data[5]
    elif data == "HNSC":
        return all_data[6]
    elif data == "KIRC":
        return all_data[7]
    elif data == "KIRP":
        return all_data[8]
    elif data == "LGG":
        return all_data[9]
    elif data == "LIHC":
        return all_data[10]
    elif data == "LUAD":
        return all_data[11]
    elif data == "LUSC":
        return all_data[12]
    elif data == "OV":
        return all_data[13]  
    elif data == "PAAD":
        return all_data[14]
    elif data == "PCPG":
        return all_data[15]
    elif data == "PRAD":
        return all_data[16]
    elif data == "READ":
        return all_data[17]
    elif data == "SARC":
        return all_data[18]
    elif data == "SKCM":
        return all_data[19]
    elif data == "STAD":
        return all_data[20]
    elif data == "TCGT":
        return all_data[21]
    elif data == "THCA":
        return all_data[22]
    elif data == "THYM":
        return all_data[23]
    elif data == "UCEC":
        return all_data[24]
    else:
        st.write("INVALID INPUT")
        return

def unique(list1):
 
    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    unique_list = (list(list_set))
    return unique_list

@st.experimental_memo
def convert_df(df):
   return df.to_csv(index=False).encode('utf-8')

for i,data in enumerate(all_data):
    all_data[i] = load_data(data)
    

which_dataset = st.selectbox("Enter Dataset: ",("BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC"))
# which_category = st.text_input("Enter Category: ")
which_category = st.selectbox("Enter Category: ",("GTMT","GTMB","GBMT","GBMB","GTCG","GTCD","GBCG","GBCD","GTMut","GBMut"))
which_threshold = st.number_input("Enter threshold: ")

if which_dataset and which_category and which_threshold:
    main_dataset = switch(which_dataset)

    tempdf = main_dataset.loc[:,["Dataset","Gene","Location","Probe","MedianBvalTumor","MedianBvalNormal","Corr",which_category]]
    tempdf = tempdf[tempdf[which_category] <= which_threshold]
    st.write("The entire dataframe filtered by category and threshold")
    st.write(tempdf)
    important_genes = tempdf["Gene"]
    important_genes = unique(important_genes)
    st.write("The significant genes")
    st.write(important_genes)

    check_gene = st.text_input("Check if gene is in this list: ")
    if check_gene:
        result = important_genes.count(check_gene)
        if result > 0:
            st.write("This gene is in the list!")
        else:
            st.write("This gene is NOT in the list!")

    important_genes = pd.DataFrame(important_genes)
    important_genes.columns = ['Gene']
    csv1 = convert_df(tempdf)
    st.download_button("Press to download full table filtered by category and threshold",csv1,"significant_genes_by_category_dataframe.csv","text/csv")
    csv2 = convert_df(important_genes)
    st.download_button("Press to download the list of significant genes",csv2,"significant_genes_by_category_genelist.csv","text/csv")
