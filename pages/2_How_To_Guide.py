import streamlit as st
from PIL import Image

st.header("A brief guide on using this website:")
st.write("Click on a tab to see a tutorial on using a feature.")

tab1,tab2,tab3 = st.tabs(["Single Gene Analysis","Significant Genes by Category","Multiple Gene Frequency"])
with tab1:
    st.subheader("Single Gene Analysis:")
    st.write("Enter a single gene, and we will give you 3 things:")
    st.write("1. A visualization of the number of significant probes per dataset and category.")
    st.write("2. A plot for each dataset plotting the category against its p value for that gene, where > 3 as a y axis label is significant (using negative log)")
    st.write("3. A plot showing the location of the probes wtih respect to the gene, with red probes being significant in some dataset.")

    st.write("Example: ERBB2: In the categorical pairs table, note that many Ovarian probes show a positive correlation between gene expression and copy number. In general, most significant probes come from CNV and increased gene expression and hypermethylation. The p value tables tell us that for some datasets like BRCA and OV, copy number seems to be the driving cause of gene expression, while for GBM and THYM, increased gene expression with hypermethylation seems to be the main cause. Finally, we see that most significant probes with respect to methylation come from a uniform region of the gene in the middle, and probes towards the end of the gene don't often cause changes in gene expression.")

with tab2:
    st.subheader("Significant Genes by Category:")
    st.write("Enter a dataset, a category, and a threshold for p values, and we will return a list of genes that are significant in that category. Furthermore, you can search to see whether a particular gene is in the list or not.")

    st.write("Example: BLCA, GTMT, 0.05: We obtain a list of genes in the BLCA dataset whose association between high gene expression and hypermethylation is significant, as they are all less than 0.05. Querying 'NXN' gives that it is a significant gene in this category, while 'BRCA1' is not.")

with tab3:
    st.subheader("Multiple Gene Frequency:")
    st.write("Enter a dataset and a list of space separated genes, comma separated genes, or a CSV file of genes, and we will return a frequency table of each category and how many genes of the list appear in each category.")
    st.write("If you enter a CSV file, please have it in this format: 1 column, column name, gene list: ")
    image1 = Image.open('./pictures/sample_CSV_input.jpg')
    st.image(image1, caption="sample CSV file to input")
    st.write("Example: All genes from Significant Genes by Category with the following selection: BLCA, GTMT, 0.05. This gives us a binary table of which genes are in each category as well as its plot. Note that all genes in this list are in the GTMT category, but these genes are also enriched in other categories, namely GBMB: its symettric counterpart. Thus, the positive relationship between the two variables is visibly present here.")
    image2 = Image.open('./pictures/frequency_counts_example.jpg')
    st.image(image2, caption="frequency table output for the given list of genes")

st.subheader("Categories:")
st.write("We test associations between 2 groups; the group breakdown is as follows:")
st.write("GT/GB: High Gene Expression/Low Gene Expression")
st.write("MT/MB: Hypermethylation/Hypomethylation")
st.write("CG/CD: Copy Number Gain/Copy Number Deletion")
st.write("Mut: Mutation")