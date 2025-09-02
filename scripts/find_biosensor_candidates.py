import pandas as pd
import argparse
import os
import re
import sys

# --- Keyword Dictionary ---
# This is the "brain" of the script. You can add new analytes and their associated
# biological keywords here to expand the script's capabilities in the future.
ANALYTE_KEYWORDS = {
    'uranium': {
        'Phosphate Sensing (Top Tier)': [
            'phosphate', 'phoR', 'phoB', 'pstS', 'phosphate-binding'
        ],
        'Heavy Metal Efflux (Good Candidate)': [
            'heavy metal', 'efflux', 'P-type ATPase', 'CzcA', 'arsC', 'cadA', 'merT'
        ],
        'General Metal Regulation (Relevant)': [
            'ArsR', 'MerR', 'metal-binding', 'metal resistance', 'Fur family'
        ]
    },
    'pfas': {
        'Dehalogenation Systems (Top Tier)': [
            'dehalogenase', 'haloacid', 'halogenated'
        ],
        'Toxin Efflux (Good Candidate)': [
            'efflux', 'multidrug resistance', 'AcrB', 'CupA', 'MFS transporter'
        ],
        'General Stress Response (Relevant)': [
            'stress protein', 'oxidative stress', 'chaperone'
        ]
    }
}

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the script."""
    parser = argparse.ArgumentParser(description="Find potential biosensor candidates for a specified analyte.")
    parser.add_argument(
        "--analyte",
        type=str,
        required=True,
        choices=ANALYTE_KEYWORDS.keys(),
        help=f"The target analyte to search for. Choices: {list(ANALYTE_KEYWORDS.keys())}"
    )
    parser.add_argument(
        "--interpro_annotations", 
        type=str, 
        required=True, 
        help="Path to the InterPro matches .tsv file from MGnify."
    )
    parser.add_argument(
        "--cluster_results", 
        type=str, 
        help="Optional path to the clustering_results.csv file to add cluster info to the report."
    )
    parser.add_argument(
        "--output_file", 
        type=str, 
        help="Optional path to save the report. Defaults to the same directory as the clustering_results.csv file."
    )
    return parser.parse_args()

def find_sensor_candidates(annotations_df: pd.DataFrame, analyte: str) -> pd.DataFrame:
    """
    Searches the annotations DataFrame for keywords related to a specified analyte.
    """
    sensor_keywords = ANALYTE_KEYWORDS.get(analyte, {})
    if not sensor_keywords:
        print(f"Error: No keywords defined for analyte '{analyte}'.")
        return pd.DataFrame()

    found_proteins = []

    for category, keywords in sensor_keywords.items():
        pattern = re.compile('|'.join(keywords), re.IGNORECASE)
        matches = annotations_df[annotations_df['signature_description'].str.contains(pattern, na=False)].copy()
        
        if not matches.empty:
            matches['Candidate_System'] = category
            found_proteins.append(matches)
    
    if not found_proteins:
        return pd.DataFrame()
        
    result_df = pd.concat(found_proteins, ignore_index=True)
    result_df = result_df[['protein_accession', 'signature_description', 'Candidate_System']].drop_duplicates()
    
    return result_df

def main():
    """Main execution function."""
    args = parse_args()

    print(f"Loading InterPro annotations from: {args.interpro_annotations}")
    try:
        annotations_df = pd.read_csv(
            args.interpro_annotations, 
            sep='\t', 
            header=0,
            usecols=['protein_accession', 'signature_description']
        )
    except Exception as e:
        print(f"Error loading annotation file: {e}")
        return

    print(f"Searching for {args.analyte.upper()} biosensor candidate signatures...")
    candidates_df = find_sensor_candidates(annotations_df, args.analyte)

    if candidates_df.empty:
        print(f"No potential {args.analyte} biosensor candidates were found based on the defined keywords.")
        return

    if args.cluster_results:
        print(f"Loading cluster results from: {args.cluster_results} to add cluster IDs.")
        clusters_df = pd.read_csv(args.cluster_results)
        clusters_df['header'] = clusters_df['header'].str.split(' ').str[0]
        
        candidates_df = pd.merge(
            candidates_df,
            clusters_df[['header', 'HDBSCAN_label']],
            left_on='protein_accession',
            right_on='header',
            how='left'
        )
        candidates_df['HDBSCAN_label'].fillna('N/A', inplace=True)
        del candidates_df['header']

    # --- Determine output file path ---
    if args.output_file:
        output_path = args.output_file
    else:
        if not args.cluster_results:
            print("Error: --output_file must be specified if --cluster_results is not provided for default path.", file=sys.stderr)
            sys.exit(1)
        output_dir = os.path.dirname(args.cluster_results)
        output_filename = f"{args.analyte}_sensor_candidates.txt"
        output_path = os.path.join(output_dir, output_filename)

    report_string = f"{args.analyte.upper()} Biosensor Candidate Report\n"
    report_string += "=" * 40 + "\n\n"
    
    for system, group in candidates_df.groupby('Candidate_System'):
        report_string += f"--- {system} ---\n"
        for _, row in group.iterrows():
            cluster_info = f"(Cluster: {int(row['HDBSCAN_label'])})" if 'HDBSCAN_label' in row and row['HDBSCAN_label'] != 'N/A' else ""
            report_string += f"  Protein: {row['protein_accession']} {cluster_info}\n"
            report_string += f"    Evidence: {row['signature_description']}\n\n"
        report_string += "\n"

    print(f"Saving report to: {output_path}")
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(report_string)

    print("Search complete.")

if __name__ == '__main__':
    main()

# clear; python3 find_biosensor_candidates.py --analyte "uranium" --cluster_results "../studies/MGYS00006491/MGYA00679207/report/mgnify/clustering_results.csv" --interpro_annotations "../studies/MGYS00006491/MGYA00679207/downloads/interPro.tsv"
