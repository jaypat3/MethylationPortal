import streamlit as st

TCGA_data_URL = "https://xenabrowser.net/datapages/"
Association_Table_Data_URL = "https://github.com/jaypat3/MethylationPortal/tree/main/TCGADatasets"

st.header("Data:")
st.subheader("All the files used on this website can be found here:")
st.write("TCGA data: [https://xenabrowser.net/datapages/](%s)" % TCGA_data_URL)
st.write("Association Table data: [https://github.com/jaypat3/MethylationPortal/tree/main/TCGADatasets](%s)" % Association_Table_Data_URL)
