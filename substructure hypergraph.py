#!/usr/bin/env python
# coding: utf-8

# In[1]:


# In[1]: imports & paths
import os, random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print("Device:", device)

LOCAL_PATH = r"C:\Users\tasmi\Downloads\cleaned_drugbank_smiles_mapping.csv"  # <-- change if needed
OUT_DIR = r"C:\Users\tasmi\Downloads"
os.makedirs(OUT_DIR, exist_ok=True)


# In[2]:


# In[2]: load dataset
df = pd.read_csv(LOCAL_PATH)

# normalize column names
df.columns = [c.lower() for c in df.columns]
if "drugbank_id" in df.columns:
    df = df.rename(columns={"drugbank_id": "drug_id"})
elif "drug_id" not in df.columns:
    raise ValueError("No drug_id column found!")

if "smiles" not in df.columns:
    raise ValueError("No smiles column found!")

# drop duplicates/NA
df = df.dropna(subset=["drug_id", "smiles"]).drop_duplicates(subset=["drug_id"])
print(f"Loaded {len(df)} drugs")
df.head()


# In[3]:


# In[3]: SMILES -> k-mers
from collections import Counter, defaultdict

KMER_K = 3
MIN_FREQ = 5

def smiles_to_kmers(smiles, k=3):
    return [smiles[i:i+k] for i in range(len(smiles)-k+1)]

# build vocabulary
vocab_counter = Counter()
for smi in df["smiles"]:
    vocab_counter.update(smiles_to_kmers(smi, k=KMER_K))

sub_vocab = [s for s, c in vocab_counter.items() if c >= MIN_FREQ]
sub2idx = {s:i for i,s in enumerate(sub_vocab)}
idx2sub = {i:s for s,i in sub2idx.items()}

print(f"Substructure vocab size: {len(sub_vocab)}")

# map drugs to substructures
drug2idx = {d:i for i,d in enumerate(df["drug_id"])}
idx2drug = {i:d for d,i in drug2idx.items()}

incidence = []
for d, smi in zip(df["drug_id"], df["smiles"]):
    j = drug2idx[d]
    for sub in smiles_to_kmers(smi, k=KMER_K):
        if sub in sub2idx:
            i = sub2idx[sub]
            incidence.append((i,j))

print(f"Positive incidences: {len(incidence)}")


# In[4]:


# In[4]: Build incidence tensor H_idx
H_idx = torch.tensor(incidence, dtype=torch.long).t().contiguous()  # shape (2, num_pos)
V, E = len(sub_vocab), len(df)
print(f"#substructures (V)={V}, #drugs (E)={E}, incidences={H_idx.shape[1]}")


# In[5]:


# In[5]: hyper-parameters
D_NODE_IN = 128
D_EDGE_IN = 128
D_HIDDEN  = 128  # same as D_OUT
D_OUT     = 128

EPOCHS = 5
BATCHES_PER_EPOCH = 200
BATCH_SIZE = 4096
LR = 1e-3
WEIGHT_DECAY = 1e-4


# In[6]:


# In[6]: HyGNN encoder (corrected)
class HyGNNEncoder(nn.Module):
    def __init__(self, d_node_in, d_edge_in, d_hidden, d_out, dropout=0.1):
        super().__init__()
        # edge -> node
        self.W1 = nn.Linear(d_edge_in, d_hidden, bias=False)
        self.W2 = nn.Linear(d_edge_in, d_hidden, bias=False)
        self.W3 = nn.Linear(d_node_in, d_hidden, bias=False)
        # node -> edge
        self.W4 = nn.Linear(d_hidden, d_hidden, bias=False)
        self.W5 = nn.Linear(d_hidden, d_hidden, bias=False)
        self.W6 = nn.Linear(d_edge_in, d_hidden, bias=False)

        self.node_proj = nn.Linear(d_hidden, d_hidden)
        self.edge_proj = nn.Linear(d_hidden, d_out)
        self.drop = nn.Dropout(dropout)
        self.leaky = nn.LeakyReLU(0.2)
        self.act = nn.ReLU()

    def forward(self, P, Q, H_idx, num_nodes, num_edges):
        node_idx, edge_idx = H_idx[0], H_idx[1]

        # edge -> node
        QW1, QW2, PW3 = self.W1(Q), self.W2(Q), self.W3(P)
        s_ij = self.leaky((QW2[edge_idx] * PW3[node_idx]).sum(-1))
        alpha = segment_softmax(s_ij, node_idx, num_nodes)
        msg_node = torch.zeros((num_nodes, QW1.size(-1)), device=P.device)
        msg_node.index_add_(0, node_idx, alpha.unsqueeze(-1)*QW1[edge_idx])
        P_new = self.act(self.node_proj(self.drop(msg_node)))

        # node -> edge
        PW4, PW5, QW6 = self.W4(P_new), self.W5(P_new), self.W6(Q)
        t_ji = self.leaky((PW5[node_idx]*QW6[edge_idx]).sum(-1))
        beta = segment_softmax(t_ji, edge_idx, num_edges)
        msg_edge = torch.zeros((num_edges, PW4.size(-1)), device=Q.device)
        msg_edge.index_add_(0, edge_idx, beta.unsqueeze(-1)*PW4[node_idx])
        Q_out = self.edge_proj(self.drop(msg_edge))
        return P_new, Q_out

def segment_softmax(scores, group_index, num_segments, eps=1e-12):
    device = scores.device
    max_per = torch.full((num_segments,), -1e30, device=device)
    for g in torch.unique(group_index):
        m = (group_index == g)
        max_per[g] = scores[m].max()
    shifted = scores - max_per[group_index]
    ex = shifted.exp()
    sum_per = torch.zeros(num_segments, device=device)
    for g in torch.unique(group_index):
        m = (group_index == g)
        sum_per[g] = ex[m].sum()
    return ex / (sum_per[group_index] + eps)


# In[7]:


# In[7]: init model/optimizer
torch.manual_seed(0); np.random.seed(0); random.seed(0)

H_idx_dev = H_idx.to(device)
node_embed = nn.Embedding(V, D_NODE_IN).to(device)
edge_embed = nn.Embedding(E, D_EDGE_IN).to(device)
encoder = HyGNNEncoder(D_NODE_IN, D_EDGE_IN, D_HIDDEN, D_OUT, dropout=0.1).to(device)

params = list(encoder.parameters()) + list(node_embed.parameters()) + list(edge_embed.parameters())
optimizer = torch.optim.Adam(params, lr=LR, weight_decay=WEIGHT_DECAY)

pos_node, pos_edge = H_idx_dev[0], H_idx_dev[1]
num_pos = pos_node.numel()

# sanity check forward
with torch.no_grad():
    P0, Q0 = node_embed.weight, edge_embed.weight
    P_test, Q_test = encoder(P0, Q0, H_idx_dev, V, E)
print("P dim:", P_test.shape, "Q dim:", Q_test.shape)


# In[8]:


# In[8]: safer negative sampling
bce = nn.BCEWithLogitsLoss()
pos_pairs = set(map(tuple, H_idx_dev.t().cpu().tolist()))

def sample_batch(batch_size=BATCH_SIZE, max_retries=3):
    idx = torch.randint(0, num_pos, (batch_size,), device=device)
    i_pos, j_pos = pos_node[idx], pos_edge[idx]
    i_neg = i_pos.clone()
    j_neg = torch.randint(0, E, (batch_size,), device=device)

    # reject true positives
    for _ in range(max_retries):
        bad = torch.tensor(
            [ (int(i_neg[t]), int(j_neg[t])) in pos_pairs for t in range(batch_size) ],
            device=device, dtype=torch.bool
        )
        if not bad.any():
            break
        j_neg[bad] = torch.randint(0, E, (bad.sum().item(),), device=device)

    return i_pos, j_pos, i_neg, j_neg


# In[9]:


# In[9]: training loop
for epoch in range(EPOCHS):
    encoder.train()
    total = 0.0
    for _ in range(BATCHES_PER_EPOCH):
        i_pos, j_pos, i_neg, j_neg = sample_batch()

        P0, Q0 = node_embed.weight, edge_embed.weight
        P_new, Q_out = encoder(P0, Q0, H_idx_dev, V, E)

        s_pos = (P_new[i_pos]*Q_out[j_pos]).sum(-1)
        s_neg = (P_new[i_neg]*Q_out[j_neg]).sum(-1)
        y = torch.cat([torch.ones_like(s_pos), torch.zeros_like(s_neg)])
        s = torch.cat([s_pos, s_neg])

        loss = bce(s, y)
        optimizer.zero_grad(); loss.backward(); optimizer.step()
        total += loss.item()
    print(f"Epoch {epoch+1}/{EPOCHS} loss={total/BATCHES_PER_EPOCH:.4f}")


# In[10]:


# In[10]: save embeddings
encoder.eval()
with torch.no_grad():
    P0, Q0 = node_embed.weight, edge_embed.weight
    _, Q_out = encoder(P0, Q0, H_idx_dev, V, E)

Z = Q_out.cpu().numpy()
npy_path = os.path.join(OUT_DIR, f"drug_embeddings_SUB_k{KMER_K}_minfreq{MIN_FREQ}_dim{D_OUT}.npy")
csv_path = os.path.join(OUT_DIR, f"drug_embeddings_SUB_k{KMER_K}_minfreq{MIN_FREQ}_dim{D_OUT}.csv")
np.save(npy_path, Z)
pd.DataFrame(Z, index=df["drug_id"]).rename_axis("drug_id").to_csv(csv_path)

print("Saved:", npy_path, "\nSaved:", csv_path)


# In[ ]:




