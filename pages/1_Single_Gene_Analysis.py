import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from math import log
from math import ceil

st.set_page_config(page_title="Single Gene Analysis")

st.title('Single Gene Analysis')

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
ProbeLink = './TCGADatasets/ProbeMapping/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy.csv'
GeneLink = './TCGADatasets/ProbeMapping/probeMap_hugo_gencode_good_hg19_V24lift37_probemap.csv'

all_data = [BLCAdata,BRCAdata,CESCdata,COADdata,ESCAdata,GBMdata,HNSCdata,KIRCdata,KIRPdata,LGGdata,LIHCdata,LUADdata,LUSCdata,OVdata,PAADdata,PCPGdata,PRADdata,READdata,SARCdata,SKCMdata,STADdata,TGCTdata,THCAdata,THYMdata,UCECdata]
#BLCA,BRCA,CESC,COAD,ESCA,GBM,HNSC,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC
#all_datasets = [BLCA,BRCA,CESC,COAD,ESCA,GBM,HNSC,KIRC,KIRP,LGG,LIHC,LUAD,LUSC,OV,PAAD,PCPG,PRAD,READ,SARC,SKCM,STAD,TGCT,THCA,THYM,UCEC]
#OV REMOVED FOR NOW
gene_data_array = []

@st.cache
def load_data(dataset):
    data = pd.read_csv(dataset)
    return data

def load_gene_data(dataset):
    gene_data = dataset.loc[dataset['Gene'] == text_input]
    if not gene_data.empty:
        # st.write(gene_data)
        # st.write("num rows:", gene_data.shape[0])
        gene_data_array.append(gene_data)
    return gene_data

@st.cache
def get_significant_probe_counts(dataset):
    cols = ["GTMT","GTMB","GBMT","GBMB","GTCG","GTCD","GBCG","GBCD","GTMut","GBMut"]
    to_add = [0]*(len(cols)+1)
    for i,col in enumerate(cols):
        tempdf = dataset[dataset[col] <= 0.05]
        if col == "GTMT":
            tempdf = tempdf[tempdf["GTMTnSig"] > 0.05]
        if col == "GTMB":
            tempdf = tempdf[tempdf["GTMBnSig"] > 0.05]
        if col == "GBMT":
            tempdf = tempdf[tempdf["GBMTnSig"] > 0.05]
        if col == "GBMB":
            tempdf = tempdf[tempdf["GBMBnSig"] > 0.05]    
        to_add[i+1] = tempdf.shape[0]

    if not all(v == 0 for v in to_add):
        to_add[0] = dataset.iat[0,0]
        return to_add
    else:
        return []

@st.cache
def get_significant_probes(dataset,col):
    tempdf = dataset[dataset[col] < 0.05]
    if col == "GTMT":
        tempdf = tempdf[tempdf["GTMTnSig"] >= 0.05]
    if col == "GTMB":
        tempdf = tempdf[tempdf["GTMBnSig"] >= 0.05]
    if col == "GBMT":
        tempdf = tempdf[tempdf["GBMTnSig"] >= 0.05]
    if col == "GBMB":
        tempdf = tempdf[tempdf["GBMBnSig"] >= 0.05]
    
    if not tempdf.empty:
        to_return = tempdf["Probe"].tolist()
        return to_return
    else:
        return []

def unique(list1):
 
    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    unique_list = (list(list_set))
    return unique_list

def combine_p_vals(dataset,col):
    tempdf = dataset.loc[:,col]
    # st.write(tempdf)
    to_return = stats.combine_pvalues(tempdf,method="fisher",weights=None)[1]
    # st.write(to_return)
    return to_return

def negln(x):
    return -log(x)

#LOADING DATA
for i,data in enumerate(all_data):
    all_data[i] = load_data(data)

probe_info = load_data(ProbeLink)
gene_info = load_data(GeneLink)


text_input = st.text_input("Enter Gene: ")
st.text("Note 1: BRCA1 and BRCA1/NBR2 are 2 separate genes")
st.text("Note 2: OV location data is incorrect due to hg18/19 error")
if text_input:
    st.write("You entered: ", text_input)

# st.subheader("Here's a dataset for testing")
# st.write(all_data[0])
if text_input:

    tab1,tab2,tab3 = st.tabs(["Visualizing Probes","Plotting p values","Plotting Probe Locations"])

    for data in all_data:
        load_gene_data(data)

    probe_count_df = pd.DataFrame()
    probe_count_matrix = []
    for data in gene_data_array:
        to_add = get_significant_probe_counts(data)
        if to_add:
            to_add = pd.DataFrame(to_add).T
            probe_count_matrix.append(to_add)

    if probe_count_matrix:
        probe_count_df = pd.concat(probe_count_matrix,ignore_index = True) 
        probe_count_df.columns = ['Dataset','GTMT','GTMB','GBMT','GBMB','GTCG','GTCD','GBCG','GBCD','GTMut','GBMut']
        probe_count_df = probe_count_df.set_index('Dataset')

    with tab1:
        if not probe_count_df.empty:
            st.subheader("Section 1: Visualizing significant Probes")
            probe_count_df.plot(kind='bar',stacked=True,title='Significant Probes by Dataset')
            plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
            plt.show()
            st.pyplot(plt)
            st.write("Note: there is (almost always) overlap between the classification groups GTMT,GTMB,etc.")
            st.write("This is because a probe's significant expression may influence or be influenced by anothers, or they may be entirely independent but both present.")
            st.write("Not just that, but the same probe may even show up as significant across multiple tumor types.")
            st.write("\n\n")

    important_probes_gene = []
    significant_probes_matrix = []
    labels = ['GTMT','GTMB','GBMT','GBMB','GTCG','GTCD','GBCG','GBCD','GTMut','GBMut']
    for label in labels:
        all_probes = []
        to_add = []
        for data in gene_data_array:
            to_add = get_significant_probes(data,label)
            if to_add:
                all_probes.extend(to_add)
        
        if all_probes:
            important_probes_gene = unique(all_probes)
            significant_probes_matrix.append(len(unique(all_probes)))
        else:
            significant_probes_matrix.append(0)

    with tab1:
        if not probe_count_df.empty:
            plt.figure()
            plt.bar(labels,significant_probes_matrix)
            plt.title('Significant Probes by Characteristic (unique)')
            plt.show()
            st.pyplot(plt)
            st.write("This plot gives counts of how many probes are significant in a certain category. Again, note that a probe may show up as significant across multiple categories.")


    p_value_df = pd.DataFrame()
    p_value_matrix = []
    order = []
    for data in gene_data_array:
        dataset_p_vals = []
        dataset_p_vals.append(data.iloc[0,0])
        order.append(data.iloc[0,0])
        for label in labels:
            p_val = combine_p_vals(data,label)
            if p_val:
                p_val = negln(p_val)
                dataset_p_vals.append(p_val)
            else:
                dataset_p_vals.append(-1)
        p_value_matrix.append(pd.DataFrame(dataset_p_vals).T)

    with tab2:
        if p_value_matrix:
            p_value_df = pd.concat(p_value_matrix,ignore_index=True)
            p_value_df.columns = ['Dataset','GTMT','GTMB','GBMT','GBMB','GTCG','GTCD','GBCG','GBCD','GTMut','GBMut']
            p_value_df = p_value_df.set_index('Dataset')
            p_value_df = p_value_df.T
            p_value_df_copy = p_value_df.copy()
            p_value_df['Categories'] = ['GTMT','GTMB','GBMT','GBMB','GTCG','GTCD','GBCG','GBCD','GTMut','GBMut']
            st.subheader("Section 2: Visualizing p values per dataset")
            st.write("Depending on the number of datasets this gene appears in, there may be many plots showing which categories a gene is significant in.\n We first provide a table showing a summary of what will be plotted")
            st.write(p_value_df_copy)
        
        for item in order:
            plt.figure()
            p_value_df.plot.scatter('Categories',item)
            plt.axhline(y = negln(0.05), color = 'r', linestyle = '-')
            plt.title(item + ' overall p value (negative log)')
            plt.show()
            st.pyplot(plt)

    with tab3:
        st.subheader("Section 3: Visualizing the location of probes inside a gene")
        probe_info = probe_info[probe_info["gene"] == text_input]
        gene_info = gene_info[gene_info["gene"] == text_input]
        gene_start = gene_info.iloc[0,3]
        gene_end = gene_info.iloc[0,4]
        probe_info.sort_values('chromStart',ascending=True,inplace=True)
        x = probe_info.loc[:,"#id"]
        y = probe_info.loc[:,"chromStart"]
        intersections = list(set(x).intersection(important_probes_gene))
        diff_color_int = np.where(x.isin(intersections),'r','b')
        plt.figure()
        plt.scatter(x,y,c=diff_color_int)
        plt.axhline(y = gene_start, color = 'r', linestyle = '-')
        plt.axhline(y = gene_end, color = 'r', linestyle = '-')
        plt.suptitle(text_input + ' Probe Location map')
        plt.title('Red probes are significant in some dataset')
        if len(x) > 10:
            plt.xticks([])
        plt.show()
        st.pyplot(plt)

