# param.py

# PubChem API configuration
CID_LOOKUP_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/JSON"
SMILES_LOOKUP_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/property/CanonicalSMILES/JSON"

#  timeout settings for API calls
API_TIMEOUT = 10
API_SLEEP_SECONDS = 0.3

# SMILES cleaning options
DROP_INCOMPLETE_ENTRIES = True

# Output filenames
CLEANED_SMILES_CSV_FILENAME = "cleaned_drugbank_smiles_mapping.csv"

# Column names for drug interaction TSV file
INTERACTION_COLUMNS = ["drug1", "drug2"]

# Data settings (optional for consistency)
EXPECTED_INTERACTION_FILE_FORMAT = "tsv"
