import os
import pandas as pd
import torch
import requests
import time
import sys
from sklearn.feature_extraction.text import CountVectorizer
from torch_geometric.data import Data
from torch_geometric.utils import to_undirected
import networkx as nx

import config
import param


def get_cid_from_drugbank(drugbank_id):
    url = f"{param.CID_LOOKUP_URL}/{drugbank_id}/cids/JSON"
    try:
        response = requests.get(url, timeout=param.API_TIMEOUT)
        if response.status_code == 200:
            return response.json()['IdentifierList']['CID'][0]
    except Exception as e:
        print(f"[CID Error] {drugbank_id}: {e}")
    return None


def get_smiles_from_cid(cid):
    url = f"{param.SMILES_LOOKUP_URL}/{cid}/property/CanonicalSMILES/JSON"
    try:
        response = requests.get(url, timeout=param.API_TIMEOUT)
        if response.status_code == 200:
            return response.json()['PropertyTable']['Properties'][0]['CanonicalSMILES']
    except Exception as e:
        print(f"[SMILES Error] CID {cid}: {e}")
    return None


def extract_and_save_smiles():
    df = pd.read_csv(config.INTERACTION_TSV, sep='\t', header=None, names=['drug1', 'drug2'])
    unique_ids = pd.unique(df[['drug1', 'drug2']].values.ravel())

    results = []
    for i, drug_id in enumerate(unique_ids):
        print(f"🔎 {i+1}/{len(unique_ids)}: {drug_id}")
        cid = get_cid_from_drugbank(drug_id)
        smiles = get_smiles_from_cid(cid) if cid else None
        results.append({'DrugBank_ID': drug_id, 'PubChem_CID': cid, 'SMILES': smiles})
        time.sleep(param.API_DELAY)

    smiles_df = pd.DataFrame(results).dropna(subset=['PubChem_CID', 'SMILES'])
    smiles_df.to_csv(config.CLEANED_SMILES_CSV, index=False)
    print(f" Saved cleaned SMILES to: {config.CLEANED_SMILES_CSV}")
    return df, smiles_df


def build_and_save_hypergraph(interactions_df, smiles_df):
    # Step 1: Vectorize SMILES
    vectorizer = CountVectorizer(analyzer='char', ngram_range=(2, 4))
    X = vectorizer.fit_transform(smiles_df['SMILES'])
    x = torch.tensor(X.toarray(), dtype=torch.float)

    # Step 2: Create edge list
    drug2idx = {drug: idx for idx, drug in enumerate(smiles_df['DrugBank_ID'])}
    edges = [
        [drug2idx[d1], drug2idx[d2]]
        for d1, d2 in zip(interactions_df['drug1'], interactions_df['drug2'])
        if d1 in drug2idx and d2 in drug2idx
    ]
    edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()
    edge_index = to_undirected(edge_index)

    # Step 3: Save graph
    data = Data(x=x, edge_index=edge_index)
    torch.save(data, config.GRAPH_PT)
    print(f"Hypergraph saved to: {config.GRAPH_PT}")
    return data


def export_graphml(interactions_df, valid_drugs):
    edge_list = [
        (row['drug1'], row['drug2'])
        for _, row in interactions_df.iterrows()
        if row['drug1'] in valid_drugs and row['drug2'] in valid_drugs
    ]
    G = nx.Graph()
    G.add_edges_from(edge_list)
    nx.write_graphml(G, config.GRAPH_GRAPHML)
    print(f" Exported GraphML to: {config.GRAPH_GRAPHML}")


def main():
    print("🚀 Starting DDI Hypergraph Pipeline")
    interactions_df, smiles_df = extract_and_save_smiles()
    data = build_and_save_hypergraph(interactions_df, smiles_df)
    export_graphml(interactions_df, set(smiles_df['DrugBank_ID']))
    print("🎉 Pipeline completed successfully!")


if __name__ == "__main__":
    main()
