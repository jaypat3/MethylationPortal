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
st.write("We utilize tables which determine the significance of association between a probe's gene expression and its epigenetic markers, and we call probes with a p value of < 0.05 in a specific category to be 'significant'. All tables used to create figures as well as the main association tables themselves are downloadable.")
st.write("A more extensive overview of what the website does and aims to accomplish is given in the about section.")
st.write("Also, guides on using the website are posted in their appropriate tab as well.")
st.write("Note: All inputs are case sensitive!")

