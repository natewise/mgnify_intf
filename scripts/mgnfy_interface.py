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
        
# Extract from gz files if needed
import gzip
def extract_gz_file(gz_filename):
    """Extract a .gz file and return the name of the extracted file."""
    extracted_filename = gz_filename[:-3]  # Remove the .gz extension
    with gzip.open(gz_filename, 'rb') as f_in:
        with open(extracted_filename, 'wb') as f_out:
            f_out.write(f_in.read())
    print(f"Extracted {gz_filename} to {extracted_filename}.")
    return extracted_filename

def parse_fasta_file(filename):
    # Check if the file exists
    if not os.path.exists(filename):
        print(f"File {filename} does not exist.")
        return
    
    # Check if the file is gzipped
    if filename.endswith('.gz'):
        print(f"File {filename} is gzipped, extracting...")
        filename = extract_gz_file(filename)

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
    # fetch_study_downloads("MGYS00001589")

    # Parse all the files in the downloads directory using the parse_fasta_file function
    downloads_dir = "downloads/MGYS00001589/MGYA00103632"
    for root, dirs, files in os.walk(downloads_dir):
        for file in files:
            if file.endswith('.fasta') or file.endswith('.fa'):
                file_path = os.path.join(root, file)
                print(f"Parsing {file_path}...")
                parse_fasta_file(file_path)
            elif file.endswith('.gz'):
                file_path = os.path.join(root, file)
                print(f"Parsing gzipped file {file_path}...")
                parse_fasta_file(file_path)
