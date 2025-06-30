#!/usr/bin/env python
# coding: utf-8

# Step 1: Build Hypergraph with `manager.py`

# In[1]:


get_ipython().system('python manager.py')


#  Step 2: Load Graph & SMILES

# In[2]:


import torch
import pandas as pd
import config


data = torch.load(config.GRAPH_PT, weights_only=False)
smiles_df = pd.read_csv(config.CLEANED_SMILES_CSV)
drug_ids = smiles_df['DrugBank_ID'].tolist()

print(f"Loaded graph with {data.num_nodes} nodes and {data.num_edges} edges.")


# In[3]:


import networkx as nx
from torch_geometric.utils import to_networkx


G = to_networkx(data, to_undirected=True)
G = nx.relabel_nodes(G, {i: drug_ids[i] for i in range(len(drug_ids))})
print(" NetworkX graph created with DrugBank ID labels.")


# Step 4: Cluster & Visualize

# In[9]:


from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import networkx as nx

# Clustering
kmeans = KMeans(n_clusters=5, random_state=42, n_init='auto')
clusters = kmeans.fit_predict(data.x.numpy())


color_map = {i: int(cluster) for i, cluster in enumerate(clusters)}


plt.figure(figsize=(14, 14))
pos = nx.spring_layout(G, seed=42)


colors = [color_map[n] for n in G.nodes()]

nx.draw(G, pos, node_color=colors, cmap=plt.cm.Set3, node_size=40, with_labels=True, font_size=6)
plt.title("Drug-Drug Graph (Clustered)")
plt.axis('off')
plt.show()




# 

# Step 5: Export to GraphML

# In[18]:


import networkx as nx

valid_edges = [
    (row['drug1'], row['drug2'])
    for _, row in pd.read_csv(config.INTERACTION_TSV, sep='\t', names=['drug1', 'drug2']).iterrows()
    if row['drug1'] in drug_ids and row['drug2'] in drug_ids
]

G_export = nx.Graph()
G_export.add_edges_from(valid_edges)
nx.write_graphml(G_export, config.GRAPH_GRAPHML)
print(f"Exported GraphML to: {config.GRAPH_GRAPHML}")


# In[ ]:




