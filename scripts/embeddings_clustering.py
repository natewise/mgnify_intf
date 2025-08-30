import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import umap
import hdbscan
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
from transformers import AutoTokenizer, AutoModel
from typing import List, Tuple, Dict, Any

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the script."""
    parser = argparse.ArgumentParser(description="Embed, reduce, and cluster protein sequences.")
    parser.add_argument("--study_dir", type=str, required=True, help="Directory containing the study data (e.g., '../results/MGYA00679207').")
    parser.add_argument("--model_name", type=str, default="facebook/esm2_t6_8M_UR50D", help="Name of the Hugging Face ESM model to use.")
    parser.add_argument("--batch_size", type=int, default=8, help="Batch size for embedding generation.")
    parser.add_argument("--n_sequences", type=int, default=None, help="Limit processing to the first N sequences (for testing).")
    parser.add_argument("--force_regenerate", action="store_true", help="Force regeneration of embeddings even if a cached file exists.")
    parser.add_argument("--k_means_clusters", type=int, default=10, help="Number of clusters for KMeans.")
    return parser.parse_args()

def load_faa(filename: str) -> Tuple[List[str], List[str]]:
    """
    Loads protein headers and sequences from a FASTA (.faa) file.

    Args:
        filename (str): Path to the FASTA file.

    Returns:
        Tuple[List[str], List[str]]: A tuple containing a list of headers and a list of sequences.
    """
    sequences, headers = [], []
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append("".join(current_seq))
                    current_seq = []
                headers.append(line[1:])
            else:
                current_seq.append(line)
        if current_seq:
            sequences.append("".join(current_seq))
            
    return headers, sequences

def generate_embeddings(sequences: List[str], model_name: str, batch_size: int, device: str) -> np.ndarray:
    """
    Generates protein embeddings using a specified ESM model.

    Args:
        sequences (List[str]): A list of protein amino acid sequences.
        model_name (str): The name of the pre-trained model from Hugging Face.
        batch_size (int): The number of sequences to process in each batch.
        device (str): The device to run the model on ('cpu' or 'cuda').

    Returns:
        np.ndarray: A numpy array of protein embeddings.
    """
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name).to(device)
    model.eval()
    
    all_embeddings = []
    with torch.no_grad():
        for i in range(0, len(sequences), batch_size):
            batch = sequences[i:i+batch_size]
            num_batches = (len(sequences) + batch_size - 1) // batch_size
            print(f"Processing batch {i//batch_size + 1}/{num_batches} {100 * (i//batch_size + 1)//num_batches}%")
            
            encoded = tokenizer(batch, return_tensors="pt", padding=True, truncation=True, max_length=1022).to(device)
            outputs = model(**encoded)
            
            # Use mean pooling of the last hidden state
            batch_embeddings = outputs.last_hidden_state.mean(dim=1).cpu().numpy()
            all_embeddings.append(batch_embeddings)
            
    return np.vstack(all_embeddings)

def reduce_dimensions(embeddings: np.ndarray, n_components: int = 2, n_neighbors: int = 15, min_dist: float = 0.1) -> np.ndarray:
    """
    Reduces the dimensionality of embeddings using UMAP.

    Args:
        embeddings (np.ndarray): The high-dimensional embedding data.
        n_components (int): The number of dimensions to reduce to.
        n_neighbors (int): UMAP n_neighbors parameter.
        min_dist (float): UMAP min_dist parameter.

    Returns:
        np.ndarray: The low-dimensional representation of the data.
    """
    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        random_state=42
    )
    return reducer.fit_transform(embeddings)

def run_clustering(embeddings: np.ndarray, k_means_clusters: int) -> Dict[str, np.ndarray]:
    """
    Performs clustering on the embeddings using multiple algorithms.

    Args:
        embeddings (np.ndarray): The embeddings to cluster (can be high or low dimensional).
        k_means_clusters (int): The number of clusters to use for KMeans.

    Returns:
        Dict[str, np.ndarray]: A dictionary mapping algorithm name to cluster labels.
    """
    clusters = {}
    
    print("Running HDBSCAN...")
    hdb = hdbscan.HDBSCAN(min_cluster_size=5, gen_min_span_tree=True)
    clusters['HDBSCAN'] = hdb.fit_predict(embeddings)

    print(f"Running KMeans with k={k_means_clusters}...")
    kmeans = KMeans(n_clusters=k_means_clusters, random_state=42, n_init='auto')
    clusters['KMeans'] = kmeans.fit_predict(embeddings)

    print("Running DBSCAN...")
    # Note: DBSCAN's eps parameter is highly sensitive and may need tuning based on the data scale.
    db = DBSCAN(eps=0.7, min_samples=5) 
    clusters['DBSCAN'] = db.fit_predict(embeddings)

    return clusters

def evaluate_clustering(embeddings: np.ndarray, cluster_results: Dict[str, np.ndarray]) -> pd.DataFrame:
    """
    Calculates clustering evaluation metrics.

    Args:
        embeddings (np.ndarray): The embeddings used for evaluation.
        cluster_results (Dict[str, np.ndarray]): The clustering results.

    Returns:
        pd.DataFrame: A DataFrame containing the evaluation metrics for each method.
    """
    metrics = {}
    for method, labels in cluster_results.items():
        # Filter out noise points (labeled -1) for fair evaluation
        valid_indices = labels != -1
        if np.sum(valid_indices) == 0: continue
        
        X_filtered = embeddings[valid_indices]
        labels_filtered = labels[valid_indices]
        
        # Metrics require at least 2 clusters to be calculated
        if len(np.unique(labels_filtered)) < 2:
            metrics[method] = {'silhouette': np.nan, 'davies_bouldin': np.nan, 'calinski_harabasz': np.nan}
            continue
            
        metrics[method] = {
            'silhouette': silhouette_score(X_filtered, labels_filtered),
            'davies_bouldin': davies_bouldin_score(X_filtered, labels_filtered),
            'calinski_harabasz': calinski_harabasz_score(X_filtered, labels_filtered)
        }
    return pd.DataFrame(metrics).T

def plot_clusters(reduced_embeddings: np.ndarray, cluster_labels: np.ndarray, method_name: str, output_dir: str):
    """
    Generates and saves a scatter plot of the clusters.

    Args:
        reduced_embeddings (np.ndarray): The 2D embeddings for plotting.
        cluster_labels (np.ndarray): The cluster labels for coloring points.
        method_name (str): The name of the clustering method for the title.
        output_dir (str): The directory to save the plot in.
    """
    plt.figure(figsize=(10, 8))
    
    # Differentiate noise points by plotting them separately in gray
    noise_mask = cluster_labels == -1
    clusters_mask = ~noise_mask
    
    plt.scatter(reduced_embeddings[clusters_mask, 0], reduced_embeddings[clusters_mask, 1], c=cluster_labels[clusters_mask], cmap='viridis', s=10, alpha=0.7)
    plt.scatter(reduced_embeddings[noise_mask, 0], reduced_embeddings[noise_mask, 1], c='lightgray', s=5, alpha=0.5, label='Noise')
    
    plt.title(f"Protein Clusters ({method_name})")
    plt.xlabel('UMAP Dimension 1')
    plt.ylabel('UMAP Dimension 2')
    plt.colorbar(label='Cluster ID')
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, f'{method_name}_clusters.png')
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Saved cluster plot to {output_file}")

def main():
    """Main execution function."""
    args = parse_args()

    if not os.path.isdir(args.study_dir):
        raise FileNotFoundError(f"Study directory not found: {args.study_dir}")

    # --- Setup Directories and Paths ---
    embeddings_dir = os.path.join(args.study_dir, "embeddings")
    plots_dir = os.path.join(args.study_dir, "plots")
    report_dir = os.path.join(args.study_dir, "report")
    os.makedirs(embeddings_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(report_dir, exist_ok=True)
    
    embeddings_file = os.path.join(embeddings_dir, f"protein_embeddings_{os.path.basename(args.model_name)}.npz")
    faa_file = os.path.join(args.study_dir, "genes/proteins.faa")

    # --- 1. Load or Generate Embeddings ---
    if os.path.exists(embeddings_file) and not args.force_regenerate:
        print(f"Loading cached embeddings from {embeddings_file}")
        data = np.load(embeddings_file, allow_pickle=True)
        embeddings = data['embeddings']
        headers = data['headers']
    else:
        if not os.path.exists(faa_file):
            raise FileNotFoundError(f"FASTA file not found: {faa_file}")
        
        print(f"Loading protein sequences from {faa_file}")
        headers, sequences = load_faa(faa_file)
        
        if args.n_sequences:
            headers = headers[:args.n_sequences]
            sequences = sequences[:args.n_sequences]
        print(f"Loaded {len(sequences)} protein sequences.")

        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        print(f"Using device: {device}")
        
        embeddings = generate_embeddings(sequences, args.model_name, args.batch_size, device)
        np.savez_compressed(embeddings_file, embeddings=embeddings, headers=np.array(headers, dtype=object))
        print(f"Saved embeddings to {embeddings_file}")

    print(f"Embeddings shape: {embeddings.shape}")

    # --- 2. Dimensionality Reduction ---
    print("Reducing dimensions with UMAP...")
    reduced_embeddings = reduce_dimensions(embeddings)

    # --- 3. Clustering ---
    print("Clustering reduced embeddings...")
    cluster_results = run_clustering(reduced_embeddings, args.k_means_clusters)

    # --- 4. Plotting and Saving Results ---
    results_df = pd.DataFrame({
        'header': headers,
        'umap1': reduced_embeddings[:, 0],
        'umap2': reduced_embeddings[:, 1],
    })

    for method, labels in cluster_results.items():
        print(f"Plotting clusters for {method}")
        plot_clusters(reduced_embeddings, labels, method, plots_dir)
        results_df[f'{method}_label'] = labels
    
    results_file = os.path.join(report_dir, "clustering_results.csv")
    results_df.to_csv(results_file, index=False)
    print(f"Saved detailed clustering results to {results_file}")

    # --- 5. Evaluation ---
    print("Evaluating clustering quality...")
    # Evaluate based on the space where clustering was performed (reduced embeddings)
    metrics_reduced = evaluate_clustering(reduced_embeddings, cluster_results)
    metrics_reduced.to_csv(os.path.join(report_dir, "cluster_metrics_reduced_space.csv"))
    print("\nMetrics based on 2D UMAP space:")
    print(metrics_reduced)

    # Optionally, evaluate based on the original high-dimensional space
    metrics_original = evaluate_clustering(embeddings, cluster_results)
    metrics_original.to_csv(os.path.join(report_dir, "cluster_metrics_original_space.csv"))
    print("\nMetrics based on original embedding space:")
    print(metrics_original)

if __name__ == '__main__':
    main()
