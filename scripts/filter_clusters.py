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
                # Remove any spaces within the sequence line itself
                current_seq.append(line.replace(" ", ""))
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

    if len(cluster_df) <= 5:
        return cluster_df['header'].tolist()
        
    centroid = cluster_df[['umap1', 'umap2']].mean().values
    distances = np.sqrt(
        (cluster_df['umap1'] - centroid[0])**2 + (cluster_df['umap2'] - centroid[1])**2
    )
    centroid_protein_header = cluster_df.loc[distances.idxmin()]['header']

    top_protein_header = cluster_df.loc[cluster_df['umap2'].idxmax()]['header']
    bottom_protein_header = cluster_df.loc[cluster_df['umap2'].idxmin()]['header']
    left_protein_header = cluster_df.loc[cluster_df['umap1'].idxmin()]['header']
    right_protein_header = cluster_df.loc[cluster_df['umap1'].idxmax()]['header']
    
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

    # Define output file paths based on the input FASTA file's location
    fasta_dir = os.path.dirname(os.path.abspath(args.fasta_file))
    faa_output_file = os.path.join(fasta_dir, f"cluster_{args.cluster_id}_representatives.faa")
    txt_output_file = os.path.join(fasta_dir, f"cluster_{args.cluster_id}_colabfold_input.txt")

    print(f"Loading cluster results from: {args.cluster_results}")
    all_clusters_df = pd.read_csv(args.cluster_results)
    
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

    # Prepare lists for both output formats
    sanitized_representatives_for_faa = []
    sequences_for_colabfold = []
    for i, header in enumerate(representative_headers):
        sequence = all_sequences.get(header)
        if sequence:
            sanitized_sequence = sequence.replace("*", "").replace("X", "").replace(" ", "")
            simple_header = f"cluster_{args.cluster_id}_protein_{i+1}"
            
            sanitized_representatives_for_faa.append((simple_header, sanitized_sequence))
            sequences_for_colabfold.append(sanitized_sequence)
        else:
            print(f"  Warning: Could not find sequence for header '{header}'")
    
    # --- Write the standard .faa file ---
    print(f"Writing standard FASTA file to: {faa_output_file}")
    with open(faa_output_file, 'w') as f:
        for simple_header, sanitized_sequence in sanitized_representatives_for_faa:
            f.write(f">{simple_header}\n")
            # Write sanitized sequence in lines of 60 characters
            for j in range(0, len(sanitized_sequence), 60):
                f.write(sanitized_sequence[j:j+60] + "\n")

    # --- Write the .txt file for ColabFold ---
    final_colabfold_string = ":".join(sequences_for_colabfold)
    print(f"Writing colon-separated sequences to: {txt_output_file}")
    with open(txt_output_file, 'w') as f:
        f.write(final_colabfold_string)
    
    print("\n--- Instructions ---")
    print(f"Two files have been created in '{fasta_dir}':")
    print(f"1. {os.path.basename(faa_output_file)}: A standard FASTA file for your records.")
    print(f"2. {os.path.basename(txt_output_file)}: A text file with sequences for ColabFold.")
    print("\nTo run ColabFold:")
    print(f"- Open '{os.path.basename(txt_output_file)}'.")
    print("- Copy the ENTIRE single line of text.")
    print("- Paste it into the 'query_sequence' box in the ColabFold notebook.")
    print("--------------------")

    print("\nScript finished successfully.")

if __name__ == '__main__':
    main()

# Example usage:
# python3 filter_clusters.py --cluster_results "../results/MGYA00679207/report/mgnify/clustering_results.csv" --fasta_file "../results/MGYA00679207/genes/mgnify/proteins.faa" --cluster_id 77

# Output is used as input for AlphaFold, which has a google colab notebook for running it for free:
# https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb?authuser=1

# Output files from AlphaFold can be analyzed with FoldSeek and ChimeraX:
# https://search.foldseek.com/search
# https://www.cgl.ucsf.edu/chimerax/download.html
# Once downloaded, load all the pdb files (File > Open) run a few commands:
# rainbow
# matchmaker #2,3,4,5 to #1 (adjust numbers depending on number of proteins)
# hide; show cartoon