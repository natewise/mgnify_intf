# Connection to MGnify API
from jsonapi_client import Session as APISession
# Data download
from urllib.request import urlretrieve
# OS operations
import os
# Dataframes and display
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)

def fetch_study_downloads(study_accession):
    endpoint = f"studies/{study_accession}/analyses"
    # Make sure the downloads directory exists
    if not os.path.exists("downloads"):
        os.makedirs("downloads")
        print("Created 'downloads' directory.")
    # Create 'study_accession' directory if it doesn't exist
    if not os.path.exists(f"downloads/{study_accession}"):
        os.makedirs(f"downloads/{study_accession}")
        print(f"Created 'downloads/{study_accession}' directory.")

    with APISession("https://www.ebi.ac.uk/metagenomics/api/v1") as mgnify:
        analyses = map(lambda r: r.json, mgnify.iterate(endpoint))
        analyses = pd.json_normalize(analyses)
        for analysisId in analyses.head(2)['attributes.accession']:
            # Create 'analysisId' directory if it doesn't exist
            downloads_dir = f"downloads/{study_accession}/{analysisId}"
            if not os.path.exists(downloads_dir):
                os.makedirs(downloads_dir)
                print(f"Created '{downloads_dir}' directory.")
            print(f"\nStarting downloads for analysis {analysisId}")
            for download in mgnify.iterate(f"analyses/{analysisId}/downloads"):
                # print(f"{download.alias}: {download.description.label}")
                print(f"Downloading {download.alias}...")
                try:
                    # Check if file already exists
                    if os.path.exists(f"{downloads_dir}/{download.alias}"):
                        print(f"{download.alias} already exists, skipping download.")
                        continue
                    urlretrieve(download.links.self.url, f"{downloads_dir}/{download.alias}")
                except Exception as e:
                    print(f"Failed to download {download.alias}: {e}")

    print("Download completed.")
        
def parse_fasta_file(filename):
    """Parse a FASTA file and return a list of sequences."""
    sequences = []
    with open(filename, 'r') as file:
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence)
    print(f"Parsed {len(sequences)} sequences from {filename}.")

if __name__ == "__main__":
    # Study focusing on microbial communities in heavy metal contaminated soils
    # https://www.ebi.ac.uk/metagenomics/studies/MGYS00001589#overview
    fetch_study_downloads("MGYS00001589")

    parse_fasta_file("downloads/MGYS00001589/MGYA00103632/ERR788946_MERGED_FASTQ_16SrRNA.fasta")
