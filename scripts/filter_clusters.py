import pandas as pd
import numpy as np
import argparse
import os
from typing import Dict, List, Tuple

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the script."""
    parser = argparse.ArgumentParser(
        description="Select representative proteins from a specific cluster based on UMAP coordinates."
    )
    parser.add_argument(
        "--cluster_results",
        type=str,
        required=True,
        help="Path to the clustering_results.csv file."
    )
    parser.add_argument(
        "--fasta_file",
        type=str,
        required=True,
        help="Path to the original .faa file containing all protein sequences."
    )
    parser.add_argument(
        "--cluster_id",
        type=int,
        required=True,
        help="The ID of the cluster to analyze (e.g., 77)."
    )
    parser.add_argument(
        "--output_file",
        type=str,
        default=None,
        help="Path to save the output FASTA file. Defaults to 'cluster_[id]_representatives.faa'."
    )
    return parser.parse_args()

def load_fasta(filename: str) -> Dict[str, str]:
    """
    Loads a FASTA file into a dictionary mapping headers to sequences.
    Cleans the header to match the format in the results CSV.
    """
    sequences = {}
    current_seq_header = None
    current_seq = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq_header:
                    sequences[current_seq_header] = "".join(current_seq)
                
                # Clean the header by removing the '>' and splitting at the first space
                current_seq_header = line[1:].split(' ')[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_seq_header:
            sequences[current_seq_header] = "".join(current_seq)
            
    return sequences

def select_representatives(cluster_df: pd.DataFrame) -> List[str]:
    """
    Selects up to 5 representative proteins from a cluster:
    1. The one closest to the centroid.
    2. The four on the extreme edges (top, bottom, left, right).
    """
    if cluster_df.empty:
        return []

    # If the cluster is small, just return all members
    if len(cluster_df) <= 5:
        return cluster_df['header'].tolist()
        
    # --- 1. Find the Centroid Protein ---
    centroid = cluster_df[['umap1', 'umap2']].mean().values
    distances = np.sqrt(
        (cluster_df['umap1'] - centroid[0])**2 + (cluster_df['umap2'] - centroid[1])**2
    )
    centroid_protein_header = cluster_df.loc[distances.idxmin()]['header']

    # --- 2. Find the Edge Proteins ---
    top_protein_header = cluster_df.loc[cluster_df['umap2'].idxmax()]['header']
    bottom_protein_header = cluster_df.loc[cluster_df['umap2'].idxmin()]['header']
    left_protein_header = cluster_df.loc[cluster_df['umap1'].idxmin()]['header']
    right_protein_header = cluster_df.loc[cluster_df['umap1'].idxmax()]['header']
    
    # --- 3. Combine and ensure uniqueness ---
    # Using a dictionary automatically handles duplicates if one protein occupies multiple roles
    representative_headers = {
        'centroid': centroid_protein_header,
        'top': top_protein_header,
        'bottom': bottom_protein_header,
        'left': left_protein_header,
        'right': right_protein_header
    }
    
    return list(representative_headers.values())


def main():
    """Main execution function."""
    args = parse_args()

    # Set default output file name if not provided
    if args.output_file:
        output_file = args.output_file
    else:
        # By default, place the output file in the same directory as the input FASTA file
        fasta_dir = os.path.dirname(args.fasta_file)
        file_name = f"cluster_{args.cluster_id}_representatives.faa"
        output_file = os.path.join(fasta_dir, file_name)

    print(f"Loading cluster results from: {args.cluster_results}")
    all_clusters_df = pd.read_csv(args.cluster_results)
    
    # Clean the header column to ensure it matches the FASTA headers
    all_clusters_df['header'] = all_clusters_df['header'].str.split(' ').str[0]

    print(f"Filtering for cluster ID: {args.cluster_id}")
    target_cluster_df = all_clusters_df[all_clusters_df['HDBSCAN_label'] == args.cluster_id]

    if target_cluster_df.empty:
        print(f"Error: Cluster {args.cluster_id} not found in the results file. Exiting.")
        return

    print(f"Found {len(target_cluster_df)} proteins in cluster {args.cluster_id}.")
    
    print("Selecting representative proteins...")
    representative_headers = select_representatives(target_cluster_df)
    print(f"Selected {len(representative_headers)} unique representative proteins.")

    print(f"Loading sequences from: {args.fasta_file}")
    all_sequences = load_fasta(args.fasta_file)

    print(f"Writing representative sequences to: {output_file}")
    with open(output_file, 'w') as f:
        for i, header in enumerate(representative_headers):
            sequence = all_sequences.get(header)
            if sequence:
                f.write(f">{header}\n")
                # Write sequence in lines of 60 characters for standard FASTA format
                for j in range(0, len(sequence), 60):
                    f.write(sequence[j:j+60] + "\n")
            else:
                print(f"  Warning: Could not find sequence for header '{header}'")
    
    print("Script finished successfully.")

if __name__ == '__main__':
    main()
