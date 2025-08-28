# protein_embedding_clustering.py

import torch
from transformers import AutoTokenizer, AutoModel
import pandas as pd
import numpy as np
import umap
import hdbscan
from sklearn.cluster import KMeans, DBSCAN
import matplotlib.pyplot as plt

# Step 1: Load protein sequences from .faa file
def load_faa(filename): 
    sequences = []
    headers = []
    with open(filename, 'r') as f:
        seq = ''
        header = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq:
                    sequences.append(seq)
                    headers.append(header)
                    seq = ''
                header = line[1:]
            else:
                seq += line
        if seq:
            sequences.append(seq)
            headers.append(header)
    return headers, sequences

# Step 2: Generate embeddings using ESM-2
def generate_embeddings(sequences, batch_size=8, device='cpu'):
    tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
    model = AutoModel.from_pretrained("facebook/esm2_t6_8M_UR50D").to(device)
    model.eval()
    embeddings = []

    with torch.no_grad():
        for i in range(0, len(sequences), batch_size):
            batch = sequences[i:i+batch_size]
            encoded = tokenizer(batch, return_tensors="pt", padding=True, truncation=True, max_length=1024)
            input_ids = encoded['input_ids'].to(device)
            attention_mask = encoded['attention_mask'].to(device)
            outputs = model(input_ids, attention_mask=attention_mask)
            batch_embeddings = outputs.last_hidden_state.mean(dim=1).cpu().numpy()
            embeddings.append(batch_embeddings)
    embeddings = np.vstack(embeddings)
    return embeddings

# Step 3: Dimensionality reduction with UMAP
def reduce_dimensions(embeddings, n_neighbors=15, min_dist=0.1, n_components=2):
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=n_components, random_state=42)
    reduced = reducer.fit_transform(embeddings)
    return reduced

# Step 4: Clustering

def run_clustering(embeddings):
    clusters = {}

    # HDBSCAN
    hdb = hdbscan.HDBSCAN(min_cluster_size=5)
    clusters['HDBSCAN'] = hdb.fit_predict(embeddings)

    # K-Means (k=10 for example)
    kmeans = KMeans(n_clusters=10, random_state=42)
    clusters['KMeans'] = kmeans.fit_predict(embeddings)

    # DBSCAN
    db = DBSCAN(eps=5, min_samples=5)
    clusters['DBSCAN'] = db.fit_predict(embeddings)

    return clusters

# Step 5: Visualization
def plot_clusters(reduced, cluster_labels, method_name):
    plt.figure(figsize=(8,6))
    plt.scatter(reduced[:,0], reduced[:,1], c=cluster_labels, cmap='tab20', s=5)
    plt.title(f"Clusters ({method_name})")
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.show()

# Main driver
def main():
    # faa_file = 'downloads/MGYS00006491/MGYA00679207/ERZ17499708_FASTA_predicted_cds.faa.gz'
    faa_file = '../results/MGYA00679207/genes/proteins.faa'
    print(f"Loading protein sequences from {faa_file}")
    headers, sequences = load_faa(faa_file)

    print(f"Loaded {len(sequences)} protein sequences")

    print("Generating embeddings...")
    embeddings = generate_embeddings(sequences)
    print(f"Generated embeddings of shape: {embeddings.shape}")

    print("Reducing dimensions...")
    reduced = reduce_dimensions(embeddings)
    print("Reduced embeddings to 2D for visualization")

    print("Clustering embeddings...")
    cluster_results = run_clustering(reduced)

    for method, labels in cluster_results.items():
        print(f"Plotting clusters for {method}")
        plot_clusters(reduced, labels, method)

if __name__ == '__main__':
    main()
