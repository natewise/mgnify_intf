import pandas as pd
import argparse
import os

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the script."""
    parser = argparse.ArgumentParser(description="Generate a detailed report for a specific protein cluster.")
    parser.add_argument(
        "--cluster_results", 
        type=str, 
        required=True, 
        help="Path to the clustering_results.csv file."
    )
    parser.add_argument(
        "--interpro_annotations", 
        type=str, 
        required=True, 
        help="Path to the InterPro matches .tsv file."
    )
    parser.add_argument(
        "--cluster_id", 
        type=int, 
        required=True, 
        help="The ID of the cluster to analyze (e.g., 31)."
    )
    parser.add_argument(
        "--output_file", 
        type=str, 
        help="Optional path to save the report. Defaults to the same directory as cluster_results."
    )
    return parser.parse_args()

def generate_cluster_report(clusters_df: pd.DataFrame, annotations_df: pd.DataFrame, cluster_id: int):
    """
    Filters for a specific cluster, merges with annotations, and returns a detailed DataFrame.
    """
    # Filter for the target cluster
    cluster_members_df = clusters_df[clusters_df['HDBSCAN_label'] == cluster_id].copy()

    # Perform a LEFT merge. This is crucial as it keeps ALL proteins from the cluster,
    # even those without a match in the annotations file.
    merged_df = pd.merge(
        cluster_members_df, 
        annotations_df, 
        left_on='header', 
        right_on='protein_id', 
        how='left'
    )
    
    # Fill in missing functional descriptions for our report
    merged_df['signature_description'].fillna('--- Unannotated (Potential Novel Candidate) ---', inplace=True)
    
    # We might have multiple annotations per protein, so we'll group by protein header
    report_data = merged_df.groupby('header')['signature_description'].apply(lambda x: ' | '.join(x.unique())).reset_index()
    
    return report_data

def main():
    """Main execution function."""
    args = parse_args()

    print(f"Loading cluster results from: {args.cluster_results}")
    clusters_df = pd.read_csv(args.cluster_results)
    
    # Clean the protein IDs in the cluster results file
    clusters_df['header'] = clusters_df['header'].str.split(' ').str[0]

    print(f"Loading InterPro annotations from: {args.interpro_annotations}")
    try:
        annotations_df = pd.read_csv(
            args.interpro_annotations, 
            sep='\t', 
            header=0,
            usecols=['protein_accession', 'signature_description']
        )
        annotations_df.rename(columns={'protein_accession': 'protein_id'}, inplace=True)
    except Exception as e:
        print(f"Error loading annotation file: {e}")
        return

    print(f"Generating detailed report for Cluster {args.cluster_id}...")
    report_df = generate_cluster_report(clusters_df, annotations_df, args.cluster_id)

    # --- Determine output file path ---
    if args.output_file:
        output_path = args.output_file
    else:
        # Default to the same directory as the cluster_results file
        output_dir = os.path.dirname(args.cluster_results)
        output_filename = f"cluster_{args.cluster_id}_deep_dive.txt"
        output_path = os.path.join(output_dir, output_filename)

    # --- Generate the final report string ---
    report_string = f"Deep Dive Report for Cluster {args.cluster_id}\n"
    report_string += "=" * 40 + "\n"
    report_string += f"Total proteins in cluster: {len(report_df)}\n\n"
    
    for _, row in report_df.iterrows():
        report_string += f"Protein: {row['header']}\n"
        report_string += f"  Function(s): {row['signature_description']}\n\n"
        
    print(f"Saving report to: {output_path}")
    # Ensure the directory exists before writing
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(report_string)

    print("Analysis complete.")

if __name__ == '__main__':
    main()

# Usage
# python deep_dive_cluster.py \
#     --cluster_results "../results/MGYA00679207/report/mgnify/clustering_results.csv" \
#     --interpro_annotations "../results/MGYA00679207/genes/mgnify/interPro.tsv" \
#     --cluster_id 31
# This will create a new file named `cluster_31_deep_dive.txt`. This file will list every one of the 30 proteins in the cluster and clearly show which ones are our known `ArsR` regulators and which ones are our unannotated, high-priority novel candidates.
