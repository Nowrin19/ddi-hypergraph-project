import os

# Root directory of the project
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

# Data directories
DATA_DIR = os.path.join(ROOT_DIR, 'data')
INPUT_DIR = os.path.join(DATA_DIR, 'input')

# Input files
INTERACTION_TSV = os.path.join(INPUT_DIR, 'ChCh-Miner_durgbank-chem-chem.tsv.gz')
CLEANED_SMILES_CSV = os.path.join(INPUT_DIR, 'cleaned_drugbank_smiles_mapping.csv')

# Output files
GRAPH_PT = os.path.join(DATA_DIR, 'hypergraph_with_smiles_features.pt')
GRAPH_HTML = os.path.join(DATA_DIR, 'drug_graph_interactive.html')
GRAPH_PNG = os.path.join(DATA_DIR, 'full_drug_graph.png')
GRAPH_GRAPHML = os.path.join(DATA_DIR, 'drug_drug_graph.graphml')

