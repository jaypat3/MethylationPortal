import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from math import log

st.set_page_config(page_title="Home")
st.title('TCGA Integrative Analysis Portal')

st.subheader("Introduction:")
st.write("Welcome! This website aims to provide analysis of TCGA (The Cancer Genome Atlas) Datasets with regards to associations between gene expression and 3 epigenetic markers of cancer: Methylation, Copy Number Variation, and Mutation.")
st.write("We utilize tables which determine the significance of association between a probe's gene expression and its epigenetic markers, and we call probes with a p value of < 0.05 in a specific category to be 'significant'.")
st.write("All data is provided in the data section of the website.")

st.header("A brief guide on using this website:")
st.subheader("Single Gene Analysis:")
st.write("Enter a single gene, and we will give you 3 things:")
st.write("1. A visualization of the number of significant probes per dataset and category.")
st.write("2. A plot for each dataset plotting the category against its p value for that gene, where > 3 as a y axis label is significant (using negative log)")
st.write("3. A plot showing the location of the probes wtih respect to the gene, with red probes being significant in some dataset.")

st.subheader("Significant Genes by Category:")
st.write("Enter a dataset, a category, and a threshold for p values, and we will return a list of genes that are significant in that category. Furthermore, you can search to see whether a particular gene is in the list or not.")

st.subheader("Multiple Gene Frequency:")
st.write("Enter a dataset and a list of space separated or comma separated genes, and we will return a frequency table of each category and how many genes of the list appear in each category.")

st.subheader("Categories:")
st.write("We test associations between 2 groups; the group breakdown is as follows:")
st.write("GT/GB: High Gene Expression/Low Gene Expression")
st.write("MT/MB: Hypermethylation/Hypomethylation")
st.write("CG/CD: Copy Number Gain/Copy Number Deletion")
st.write("Mut: Mutation")

