# Drug-Discovery-Filter
Drug discovery tool to filter and display compounds based on pharmacokinetic attributes 
Deployed at: https://drugdiscoveryfilter.herokuapp.com/

## Description
This web application uses a list of FDA-approved drugs (https://www.cureffi.org/wp-content/uploads/2013/10/drugs.txt) and allows the users
to filter based on a set of criteria: Molecular Weight, LogP (solubility), Number of Hydrogen Donors/Acceptors. These criteria are outlined by
Lipinski's Rule of Five, which are a set of guidelines that determine if an orally administered drug is more likely to be active in the human body.
The user can adjust these filters by using the sliders located in the sidebar. The resulting compounds are displayed using a Pandas dataframe. In addition,
images of the compounds are displayed using the mols2grid library.

## How to Run the Project
1. Clone the repository with `git clone`
2. Install dependencies using `pip install -r requirements.txt`
3. Start app with `streamlit run main.py`

## Technologies
- RDKit
- Mols2Grid
- Pandas
- Streamlit

# Preview
![alt text](https://github.com/lincolnwu/Drug-Discovery-Filter/blob/master/dsf_homepage1.png)
![alt text](https://github.com/lincolnwu/Drug-Discovery-Filter/blob/master/dsf_homepage2.png)
