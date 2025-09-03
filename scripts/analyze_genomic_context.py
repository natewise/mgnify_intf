import pandas as pd
import argparse
import os
import re
from typing import Dict, List, Optional

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Analyze the genomic context of biosensor candidates.")
    parser.add_argument(
        "--candidates_report",
        type=str,
        required=True,
        help="Path to the biosensor candidate report file (e.g., uranium_sensor_candidates.txt)."
    )
    parser.add_argument(
        "--contigs_fasta",
        type=str,
        required=True,
        help="Path to the FASTA file containing the assembled contigs (e.g., 'contigs.fasta')."
    )
    parser.add_argument(
        "--interpro_annotations",
        type=str,
        required=True,
        help="Path to the InterPro matches .tsv file for annotating all proteins on the contig."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        help="Directory to save reports and FASTA files. Defaults to the same directory as the candidates_report."
    )
    parser.add_argument(
        "--max_candidates",
        type=int,
        required=False,
        default=None,
        help="Optional: Limit the analysis to the top N candidates found in the report. Defaults to all."
    )
    return parser.parse_args()

def parse_candidate_report(report_file: str, max_candidates: Optional[int] = None) -> List[str]:
    """Extracts protein IDs from the candidate report, with an optional limit."""
    protein_ids = []
    with open(report_file, 'r') as f:
        for line in f:
            if line.strip().startswith("Protein:"):
                protein_ids.append(line.strip().split()[1])
    
    unique_ids = list(set(protein_ids))
    
    if max_candidates is not None and len(unique_ids) > max_candidates:
        print(f"Found {len(unique_ids)} candidates. Limiting analysis to the first {max_candidates}.")
        return unique_ids[:max_candidates]
        
    return unique_ids

def get_contig_id_from_protein(protein_id: str) -> str:
    """Parses the contig ID from a Prodigal-formatted protein ID."""
    match = re.search(r'-(NODE-[^_\s]+)', protein_id)
    if match:
        return match.group(1)
    return None

def find_genes_on_contig(all_protein_ids: List[str], target_contig_id: str) -> List[str]:
    """Finds all protein IDs that belong to a specific contig."""
    genes_on_contig = []
    contig_pattern = f"-{target_contig_id}_"
    for prot_id in all_protein_ids:
        if contig_pattern in prot_id:
            genes_on_contig.append(prot_id)
    return sorted(genes_on_contig)

def get_annotations(protein_ids: List[str], annotations_df: pd.DataFrame) -> Dict[str, str]:
    """Gets annotations for a list of proteins, safely handling missing values."""
    annotations = {}
    subset_df = annotations_df[annotations_df['protein_accession'].isin(protein_ids)]
    for protein_id, group in subset_df.groupby('protein_accession'):
        unique_descriptions = group['signature_description'].dropna().astype(str).unique()
        annotations[protein_id] = ' | '.join(unique_descriptions)
    return annotations
    
def extract_contig_fasta(contigs_fasta: str, contig_id: str, output_path: str):
    """Extracts a single contig sequence to a new FASTA file."""
    with open(contigs_fasta, 'r') as infile, open(output_path, 'w') as outfile:
        found = False
        for line in infile:
            if line.startswith('>') and contig_id in line:
                found = True
                outfile.write(line)
            elif found and line.startswith('>'):
                break
            elif found:
                outfile.write(line)

def main():
    """Main execution function."""
    args = parse_args()

    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = os.path.dirname(os.path.abspath(args.candidates_report))
    
    os.makedirs(output_dir, exist_ok=True)
    
    fasta_subdir = os.path.join(output_dir, "contig_fastas_for_blast")
    os.makedirs(fasta_subdir, exist_ok=True)

    print("Step 1: Parsing candidate report...")
    candidate_proteins = parse_candidate_report(args.candidates_report, args.max_candidates)
    if not candidate_proteins:
        print("No protein IDs found in the report. Exiting.")
        return

    print(f"Analyzing {len(candidate_proteins)} unique candidates.")

    print("Loading all InterPro annotations...")
    try:
        annotations_df = pd.read_csv(
            args.interpro_annotations, sep='\t', header=0,
            usecols=['protein_accession', 'signature_description']
        )
        all_protein_ids = annotations_df['protein_accession'].unique().tolist()
    except Exception as e:
        print(f"Error reading annotations file: {e}")
        return

    all_reports = []
    processed_contigs = set()

    for candidate in candidate_proteins:
        contig_id = get_contig_id_from_protein(candidate)
        if not contig_id or contig_id in processed_contigs:
            continue
        
        print(f"\nAnalyzing context for contig: {contig_id}")
        processed_contigs.add(contig_id)

        contig_gene_ids = find_genes_on_contig(all_protein_ids, contig_id)
        contig_annotations = get_annotations(contig_gene_ids, annotations_df)
        
        report = f"Genomic Neighborhood for Contig: {contig_id}\n"
        report += "=" * (30 + len(contig_id)) + "\n"
        report += f"(Found {len(contig_gene_ids)} total proteins on this contig)\n\n"
        
        for prot_id in contig_gene_ids:
            annotation = contig_annotations.get(prot_id, "--- Unannotated ---")
            marker = "<- YOUR CANDIDATE" if prot_id == candidate else ""
            report += f"  - Protein: {prot_id} {marker}\n"
            report += f"    Function: {annotation}\n"
        all_reports.append(report)

        contig_fasta_path = os.path.join(fasta_subdir, f"{contig_id}_for_blast.fasta")
        extract_contig_fasta(args.contigs_fasta, contig_id, contig_fasta_path)
        print(f"  - Extracted contig sequence to: {contig_fasta_path}")

    summary_path = os.path.join(output_dir, "genomic_context_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("\n\n".join(all_reports))
    print(f"\nSaved combined neighborhood analysis to: {summary_path}")
    print("\nAnalysis complete.")

if __name__ == '__main__':
    main()

# python3 analyze_genomic_context.py \
#     --candidates_report "../studies/MGYS00006491/MGYA00679207/report/mgnify/uranium_sensor_candidates.txt" \
#     --contigs_fasta "../studies/MGYS00006491/MGYA00679207/downloads/contigs.fasta" \
#     --interpro_annotations "../studies/MGYS00006491/MGYA00679207/downloads/interPro.tsv"
#     --max_candidates 50

# **Step 2: Analyze the "Genomic Neighborhood" Report**

# The script will create a file named `genomic_context_summary.txt`. This is your primary result. Look for patterns:
# * Is your candidate protein next to a **P-type ATPase**?
# * Is it next to a known **phosphate-binding protein (PstS)**?
# * Is it surrounded by other **unannotated** proteins, suggesting a novel operon?

# **Step 3: Identify the Microbe via BLAST**

# The script will also create a `genomic_context` directory containing several `.fasta` files, one for each contig that had a candidate. Now you perform the final step manually:

# 1.  Go to the [NCBI Nucleotide BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=MegaBlast) website.
# 2.  Click the "Choose File" button and upload one of your generated `.fasta` files (e.g., `NODE-2603..._for_blast.fasta`).
# 3.  Keep the default database (`nt`) and click the "BLAST" button at the bottom.
# 4.  The results page will show you which known bacteria have the most similar DNA sequence, giving you the likely taxonomic identity of your contig.

# This workflow seamlessly connects our cluster analysis to a deep, gene-by-gene investigation, bringing you to the very brink of confirming the function and origin of your novel biosensor candidates.