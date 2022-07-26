# import mols2grid
import pandas as pd 
import streamlit as st 
import streamlit.components.v1 as components
from rdkit import Chem
# from rdkit.Chem.Draw import rdDepictor, rdMolDraw2D
from rdkit.Chem.Descriptors import ExactMolWt, MolLogP, NumHDonors, NumHAcceptors

st.title("Filter FDA-Approved Drugs by Lipinski's Rule of Five")

st.markdown("What is Lipinski's Rule of Five?")

# Use cache to allow for faster reloads after each update
@st.cache(allow_output_mutation=True)

# Read data into pandas df
def download_dataset():
    # Load data once, and then cache for future loads
    df = pd.read_csv(
        # Read data separated by tabs, and drop any missing values
        # File contains list of FDA-approved drops and SMILES notations
        "https://www.cureffi.org/wp-content/uploads/2013/10/drugs.txt", sep="\t"
    ).dropna()
    return df

'''Descriptors used for Lipinski's Rule of Five
   Each function will accept a SMILES string, which is a text-representation of a molecular structure
   Using rdkit
'''
# Calculate molecular weight
def calculate_molWeight(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    return ExactMolWt(mol)

# Calculate logP (solubility)
def calculate_logP(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    return MolLogP(mol)

# Calculate number of hydrogen donors
def calculate_hydroDonors(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    return NumHDonors(mol)

# Calculate number of hydrogen acceptors
def calculate_hydroAcceptors(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    return NumHAcceptors(mol)

# Copy the data so that changes aren't made to the original cached version
df = download_dataset().copy()

# Perform caclulation on each drug and display df columns
# Lambda function takes in the dataset, and extracts the "smiles" column, and applies the function
df["MW"] = df.apply(lambda x: calculate_molWeight(x["smiles"]), axis=1)
df["LogP"] = df.apply(lambda x: calculate_logP(x["smiles"]), axis=1)
df["NumHDonors"] = df.apply(lambda x: calculate_hydroDonors(x["smiles"]), axis=1)
df["NumHAcceptors"] = df.apply(lambda x: calculate_hydroAcceptors(x["smiles"]), axis=1)

# Create Sidebar panel 
st.sidebar.header('Set Parameters')
st.sidebar.write('*Note: Display compounds having values less than the following thresholds*')

# Molecular Weight Slider
molWeight_cutoff = st.sidebar.slider(
    label="Molecular Weight",
    min_value=0,
    max_value=1000,
    value=500, # Default Value
    step=10,
)

# LogP Slider
logP_cutoff = st.sidebar.slider(
    label="LogP",
    min_value=-10,
    max_value=10,
    value=5,
    step=1,
)

# Hydrogen Donors Slider
hydrogenDonors_cutoff = st.sidebar.slider(
    label="Number of Hydrogen Donors",
    min_value=0,
    max_value=15,
    value=5,
    step=1,
)

# Hydrogen Acceptors Slider
hydrogenAcceptors_cutoff = st.sidebar.slider(
    label="Number of Hydrogen Acceptors",
    min_value=0,
    max_value=20,
    value=10,
    step=1,
)

# Apply the filters/cutoff values from sidebar to the dataframe, and update the df accordingly
# Take the first filter, and shift it to the following filter in order to apply all 4 filters simultaneously
# The number of compounds will continuously decrease as filters narrow the search
df_result = df[df["MW"] < molWeight_cutoff]
df_result2 = df_result[df_result["LogP"] < logP_cutoff]
df_result3 = df_result2[df_result2["NumHDonors"] < hydrogenDonors_cutoff]
df_result4 = df_result3[df_result3["NumHAcceptors"] < hydrogenAcceptors_cutoff]

# Print out the number of drugs at the end of the search
st.write(df_result4.shape)
# Display the dataframe
st.write(df_result4)


# Display the compounds using mol2grid by embedded html inside a window
# Using the final filtered df (df_result4)
raw_html = mols2grid.display(df_result4,
    subset=["img", "Name", "MW", "LogP", "NumHDonors", "NumHAcceptors"],
    # Name of columns in dataframe
    mapping = {"smiles": "SMILES", "generic_name": "Name"})._repr_html_()

components.html(raw_html, width=900, height=1100, scrolling=False)