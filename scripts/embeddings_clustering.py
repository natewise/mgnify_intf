# protein_embedding_clustering.py

import torch
from transformers import AutoTokenizer, AutoModel
import pandas as pd
import numpy as np
import umap
import hdbscan
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
import matplotlib.pyplot as plt
import os

STUDY_DIR = f"{os.path.dirname(os.path.abspath(__file__))}/../results/MGYA00679207"
# check if directory exists, otherwise throw error
if not os.path.exists(STUDY_DIR):
    raise ValueError(f"STUDY_DIR {STUDY_DIR} does not exist")

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
            print(f"Processing batch {i//batch_size + 1} / {(len(sequences) + batch_size - 1)//batch_size} {100 * (i//batch_size + 1) // ((len(sequences) + batch_size - 1)//batch_size)}%")
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

# Step 5: Visualization (save to file)
def plot_clusters(reduced, cluster_labels, method_name, output_dir='plots'):
    os.makedirs(output_dir, exist_ok=True)
    plt.figure(figsize=(8,6))
    plt.scatter(reduced[:,0], reduced[:,1], c=cluster_labels, cmap='tab20', s=5)
    plt.title(f"Clusters ({method_name})")
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.tight_layout()
    output_file = os.path.join(output_dir, f'{method_name}_clusters.png')
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Saved cluster plot to {output_file}")

# Step 6: Evaluate clustering
def evaluate_clustering(embeddings, cluster_results):
    metrics = {}
    for method, labels in cluster_results.items():
        # Keep only points that are not noise (-1)
        valid_idx = labels != -1
        y = labels[valid_idx]
        X = embeddings[valid_idx]

        # Check if there are at least 2 clusters
        if len(np.unique(y)) < 2:
            metrics[method] = {'silhouette': np.nan,
                               'davies_bouldin': np.nan,
                               'calinski_harabasz': np.nan}
            continue

        metrics[method] = {
            'silhouette': silhouette_score(X, y),
            'davies_bouldin': davies_bouldin_score(X, y),
            'calinski_harabasz': calinski_harabasz_score(X, y)
        }
    return metrics


# Main driver
def main():
    embeddings_dir = f"{STUDY_DIR}/embeddings"
    os.makedirs(embeddings_dir, exist_ok=True)
    embeddings_file = f"{embeddings_dir}/protein_embeddings.npz"

    if os.path.exists(embeddings_file):
        print(f"Loading embeddings from {embeddings_file}")
        embeddings = np.load(embeddings_file)['embeddings']
    else:
        faa_file = f"{STUDY_DIR}/genes/proteins.faa"
        if not os.path.exists(faa_file):
            raise ValueError(f"FASTA file {faa_file} does not exist")

        print(f"Loading protein sequences from {faa_file}")
        headers, sequences = load_faa(faa_file)
        sequences = sequences[:1000]  # Limit for testing
        print(f"Loaded {len(sequences)} protein sequences")

        print("Generating embeddings...")
        embeddings = generate_embeddings(sequences)
        np.savez_compressed(embeddings_file, embeddings=embeddings)
        print(f"Saved embeddings to {embeddings_file}")

    print(f"Embeddings shape: {embeddings.shape}")

    print("Reducing dimensions...")
    reduced = reduce_dimensions(embeddings)
    print("Reduced embeddings to 2D for visualization")

    print("Clustering embeddings...")
    plots_dir = f"{STUDY_DIR}/plots"
    os.makedirs(plots_dir, exist_ok=True)
    cluster_results = run_clustering(reduced)

    for method, labels in cluster_results.items():
        print(f"Plotting clusters for {method}")
        plot_clusters(reduced, labels, method, plots_dir)

    print("Evaluating clustering quality...")
    metrics = evaluate_clustering(embeddings, cluster_results)
    metrics_df = pd.DataFrame(metrics).T
    metrics_file = f"{STUDY_DIR}/report/cluster_metrics.csv"
    metrics_df.to_csv(metrics_file)
    print(f"Saved clustering evaluation metrics to {metrics_file}")
    print(metrics_df)

if __name__ == '__main__':
    main()
