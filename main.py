import mols2grid
import pandas as pd 
import streamlit as st 
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt, MolLogP, NumHDonors, NumHAcceptors

st.title("Filter FDA-Approved Drugs by Lipinski's Rule of Five")

st.markdown("What is Lipinski's Rule of Five?")

# Use cache to allow for faster reloads after each update
@st.cache(allow_output_mutation=True)
def download_dataset():
    # Load data once, and then cache for future loads
    df = pd.read_csv(
        # Read data separated by tabs, and drop any missing values
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
    return ExtractMolWt(mol)

# Calculate logP (solubility)
def calculate_logP(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    return calculate_logP

# Calculate number of hydrogen donors
def calculate_hydroDonors(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    return NumHDonors(mol)

# Calculate number of hydrogen acceptors
def calculate_hydroAcceptors(smiles_string):
    mol = Chem.MoleFromSmiles(smile_string)
    return NumHAcceptors(mol)
