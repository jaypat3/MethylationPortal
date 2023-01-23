import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from math import log
from math import ceil
import io

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
all_data_names = ["BLCA","BRCA","CESC","COAD","ESCA","GBM","HNSC","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC"]
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

@st.experimental_memo
def convert_df(df):
   return df.to_csv(index=False).encode('utf-8')

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
    else:
        st.error("Invalid input; note that inputs are case sensitive! (Datasets are BLCA, BRCA, CESC,etc and genes are BRCA1,NXN,ERBB2, etc)")
        st.stop()
    with tab1:
        if not probe_count_df.empty:
            st.subheader("Section 1: Visualizing significant Probes")
            plt.figure()
            probe_count_df.plot(kind='bar',stacked=True,title='Significant Probes by Dataset')
            plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
            plt.ylabel('Frequency')
            plt.show()
            st.pyplot(plt)
            st.write("Note: there is (almost always) overlap between the classification groups GTMT,GTMB,etc.")
            st.write("This is because a probe's significant expression may influence or be influenced by anothers, or they may be entirely independent but both present.")
            st.write("Not just that, but the same probe may even show up as significant across multiple tumor types.")
            st.write("\n\n")

            
            x = probe_count_df.index
            plt.figure()
            plt.suptitle("Probe Frequencies across Categorical Pairs")
            plt.subplot(2,2,1)
            y1 = probe_count_df.loc[:,"GTMT"]
            y2 = probe_count_df.loc[:,"GTMB"]
            plt.bar(x,y1,color="r")
            plt.bar(x,y2,bottom=y1,color="b")
            plt.legend(["GTMT","GTMB"],bbox_to_anchor=(1.05, 1.0),loc = 'upper left')
            plt.xlabel('Dataset')
            plt.ylabel('Frequency')

            plt.subplot(2,2,2)
            y1 = probe_count_df.loc[:,"GBMT"]
            y2 = probe_count_df.loc[:,"GBMB"]
            plt.bar(x,y1,color="r")
            plt.bar(x,y2,bottom=y1,color="b")
            plt.legend(["GBMT","GBMB"],bbox_to_anchor=(1.05, 1.0),loc='upper left')
            plt.xlabel('Dataset')
            plt.ylabel('Frequency')

            plt.subplot(2,2,3)
            y1 = probe_count_df.loc[:,"GTCG"]
            y2 = probe_count_df.loc[:,"GTCD"]
            plt.bar(x,y1,color="r")
            plt.bar(x,y2,bottom=y1,color="b")
            plt.legend(["GTCG","GTCD"],bbox_to_anchor=(1.05, 1.0),loc='upper left')
            plt.xlabel('Dataset')
            plt.ylabel('Frequency')

            plt.subplot(2,2,4)
            y1 = probe_count_df.loc[:,"GBCG"]
            y2 = probe_count_df.loc[:,"GBCD"]
            plt.bar(x,y1,color="r")
            plt.bar(x,y2,bottom=y1,color="b")
            plt.legend(["GBCG","GBCD"],bbox_to_anchor=(1.05, 1.0),loc='upper left')
            plt.xlabel('Dataset')
            plt.ylabel('Frequency')

            plt.subplots_adjust(wspace=1,hspace=1)
            plt.show()
            st.pyplot(plt)

            csv = convert_df(probe_count_df)
            st.download_button("Press to download significant probe counts table",csv,"significant_probes_dataset.csv","text/csv")


    important_probes_gene = []
    significant_probes_matrix_counts = []
    important_probes_gene_methylation = []
    significant_probes_df = pd.DataFrame()
    labels = ['GTMT','GTMB','GBMT','GBMB','GTCG','GTCD','GBCG','GBCD','GTMut','GBMut']
    for label in labels:
        all_probes = []
        to_add = []
        for data in gene_data_array:
            to_add = get_significant_probes(data,label)
            if to_add:
                all_probes.extend(to_add)
            if to_add and (label == 'GTMT' or label == 'GTMB' or label == 'GBMT' or label == 'GBMB'):
                important_probes_gene_methylation.extend(to_add)

        if all_probes:
            important_probes_gene = unique(all_probes)
            significant_probes_df = pd.concat([significant_probes_df,pd.Series(important_probes_gene,name=label)],axis=1)
            significant_probes_matrix_counts.append(len(unique(all_probes)))
        else:
            significant_probes_matrix_counts.append(0)

    important_probes_gene_methylation = unique(important_probes_gene_methylation)
    with tab1:
        if not probe_count_df.empty:
            plt.figure()
            plt.bar(labels,significant_probes_matrix_counts)
            plt.title('Significant Probes by Characteristic (unique)')
            plt.ylabel('Frequency')
            plt.show()
            st.pyplot(plt)
            st.write("This plot gives counts of how many probes are significant in a certain category. Again, note that a probe may show up as significant across multiple categories.")
        
        csv = convert_df(significant_probes_df)
        st.download_button("Press to download unique probes category table",csv,"significant_probes_category.csv","text/csv")




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
            categories = ['GTMT','GTMB','GBMT','GBMB','GTCG','GTCD','GBCG','GBCD','GTMut','GBMut']
            p_value_df['Categories'] = categories
            st.subheader("Section 2: Visualizing p values per dataset")

        plt.figure(figsize=(25,25))
        plt.suptitle("P Values per dataset",y=0.95,size=32)
        x = categories
        for index,item in enumerate(all_data_names,1):
            plt.subplot(5,5,index,)
            plt.title(item)
            default_x_ticks = range(len(x))
            plt.xticks(default_x_ticks, x)
            plt.xlabel('Category')
            plt.ylabel('p value (-ln)')
            if item in order:
                y = p_value_df.loc[:,item]
                plt.scatter(x,y)
                plt.axhline(y = negln(0.05),color='r',linestyle='-')
        
        
        plt.subplots_adjust(wspace=0.25,hspace=0.25)
        plt.show()
        st.pyplot(plt) 

        csv = convert_df(p_value_df)
        st.download_button("Press to download p value table",csv,"p_value_table.csv","text/csv")

        fn = 'p_values.png'
        img = io.BytesIO()
        plt.savefig(img, format='png')
 
        btn = st.download_button(
            label="Download image",
            data=img,
            file_name=fn,
            mime="image/png"
        )

        # for item in order:
            # plt.figure()
            # p_value_df.plot.scatter('Categories',item)
            # plt.axhline(y = negln(0.05), color = 'r', linestyle = '-')
            # plt.title(item + ' overall p value (negative log)')
            # plt.show()
            # st.pyplot(plt)

    with tab3:
        st.subheader("Section 3: Visualizing the location of probes inside a gene")
        probe_info = probe_info[probe_info["gene"] == text_input]
        gene_info = gene_info[gene_info["gene"] == text_input]
        gene_start = gene_info.iloc[0,3]
        gene_end = gene_info.iloc[0,4]
        probe_info.sort_values('chromStart',ascending=True,inplace=True)
        x = probe_info.loc[:,"#id"]
        y = probe_info.loc[:,"chromStart"]
        intersections = list(set(x).intersection(important_probes_gene_methylation))
        diff_color_int = np.where(x.isin(intersections),'r','b')
        plt.figure()
        plt.scatter(y,x,c=diff_color_int)
        plt.axvline(x = gene_start, color = 'r', linestyle = '-')
        plt.axvline(x = gene_end, color = 'r', linestyle = '-')
        plt.suptitle(text_input + ' Probe Location map')
        plt.xlabel('Location (BP)')
        plt.ylabel('Probe')
        plt.title('Red probes are significant (methylation-wise) in a dataset')
        if len(y) > 10:
            plt.yticks([])
        plt.show()
        st.pyplot(plt)

