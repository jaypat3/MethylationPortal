import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from math import log

st.title('Methylation Portal')

st.write("A brief guide on using this website:")
st.subheader("singlegene:")
st.write("Enter a single gene, and we will give you 3 things:")
st.write("1. A description of all the significant probes for the gene based on category and dataset.")
st.write("2. A plot for each dataset plotting the category against its p value for that gene, where > 3 is significant")
st.write("3. A plot showing the location of the probes wtih respect to the gene, with red probes being significant in some dataset.")

st.subheader("significantgenes:")
st.write("Enter a dataset, a category, and a threshold for p values, and we will return a list of genes that are significant in that category.")
st.write("Furthermore, you can search to see whether a particular gene is in the list or not.")

st.subheader("genefreqtable:")
st.write("Enter a dataset and a list of space separated genes (will change this to comma separated if that's more convenient, let me know).")
st.write("We will return a frequency table of each category and how many genes of the list appear in each category.")

st.subheader("Categories:")
st.write("We test associations between 2 groups; the group breakdown is as follows:")
st.write("GT/GB: High Gene Expression/Low Gene Expression")
st.write("MT/MB: Hypermethylation/Hypomethylation")
st.write("CG/CD: Copy Number Gain/Copy Number Deletion")
st.write("Mut: Mutation")

