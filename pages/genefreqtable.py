import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from math import log


st.title('Frequency Table Across Genes')

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
all_categories = ["GTMT","GTMB","GBMT","GBMB","GTCG","GTCD","GBCG","GBCD","GTMut","GBMut"]

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

def get_count_binary(gene,data):
    arr = []
    data = data.loc[data['Gene'] == gene]
    for category in all_categories:
        temp = data[data.loc[:,category] <= 0.05]
        if category == "GTMT":
            temp = temp[temp["GTMTnSig"] >= 0.05]
        if category == "GTMB":
            temp = temp[temp["GTMBnSig"] >= 0.05]
        if category == "GBMT":
            temp = temp[temp["GBMTnSig"] >= 0.05]
        if category == "GBMB":
            temp = temp[temp["GBMBnSig"] >= 0.05]


        if temp.empty:
            arr.append(0)
        else:
            arr.append(1)
    return arr
    

for i,data in enumerate(all_data):
    all_data[i] = load_data(data)

which_dataset = st.text_input("Enter Dataset: ")
input = st.text_input("Enter list of genes (space separated)")

if which_dataset and input:
    main_dataset = switch(which_dataset)
    genes_list = input.split()
    L = []
    for gene in genes_list:
        to_add = []
        to_add.append(gene)
        category_nums = get_count_binary(gene,main_dataset)
        to_add.extend(category_nums)
        L.append(pd.DataFrame(to_add).T)
    
    df = pd.concat(L,ignore_index=True)
    x = ['Dataset','GTMT','GTMB','GBMT','GBMB','GTCG','GTCD','GBCG','GBCD','GTMut','GBMut']
    df.columns = x
    df = df.set_index('Dataset')
    st.write(df)
    st.write("Note: The above is a binary matrix; 1 for has significant probes, 0 for doesn't.")
    st.write("The idea of the plot below is to get an idea of how many genes are represented in a category, not how many probes, as that may skew the plots since some genes may have many more probes than others.")

    x.pop(0)
    y = df.sum(axis=0)
    plt.figure()
    plt.title("Frequency Counts of Category in Gene List")
    plt.bar(x,y)
    plt.show()
    st.pyplot(plt)
    
