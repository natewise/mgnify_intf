import os
import subprocess
import shutil
from scripts.mgnfy_interface import fetch_study_downloads

STUDY_ACCESSION = "MGYS00006491"
STUDY_DIR = f"downloads/{STUDY_ACCESSION}"
RESULTS = "results"
HMMs = ["hmm/PF01022.hmm", "hmm/PF00376.hmm"]  # ArsR/SmtB, MerR

os.makedirs(RESULTS, exist_ok=True)

############################################
# Helper functions
############################################

def list_analyses():
    """Return a list of analysis folder names under STUDY_DIR."""
    return [d for d in os.listdir(STUDY_DIR) if os.path.isdir(os.path.join(STUDY_DIR, d))]

def detect_format(path):
    if path.endswith((".fastq.gz", ".fq.gz")):
        return "fastq"
    elif path.endswith((".fasta.gz", ".fa.gz", ".fna.gz")):
        return "fasta"
    else:
        raise ValueError(f"Unknown file format: {path}")

def get_fasta_files_from_dir(dir):
    files = sorted([os.path.join(dir,f) for f in os.listdir(dir) if f.endswith(".fasta.gz") or f.endswith(".fasta")])
    if not files:
        raise ValueError(f"No fasta files found in {dir}")
    return files

def get_MGnfy_single_end_fasta_sample_reads(analysis_dir):
    fasta_files = get_fasta_files_from_dir(analysis_dir)
    # extract single file with 'MERGED_FASTQ.fasta.gz' ending
    files = [f for f in fasta_files if f.endswith(".fasta.gz")]
    if len(files) > 1 or len(files) == 0:
        raise ValueError(f"Expecting only 1 file with 'MERGED_FASTQ.fasta.gz' ending in {analysis_dir}, found {len(files)}")
    return files

############################################
# Pipeline steps
############################################

def try_run_fastp(reads, outdir):
    os.makedirs(outdir, exist_ok=True)
    fmt = detect_format(reads[0])
    json_out = os.path.join(outdir, "fastp.json")
    html_out = os.path.join(outdir, "fastp.html")
    R1_out = None
    R2_out = None

    if fmt == "fastq":
        R1_out = os.path.join(outdir, "clean_R1.fastq.gz")
        R2_out = os.path.join(outdir, "clean_R2.fastq.gz")
        if len(reads) > 1:
            cmd = ["fastp", "-i", reads[0], "-I", reads[1], "-o", R1_out, "-O", R2_out, "-j", json_out, "-h", html_out]
        else:
            cmd = ["fastp", "-i", reads[0], "-o", R1_out, "-j", json_out, "-h", html_out]
            if not os.path.exists(R2_out):
                os.symlink(R1_out, R2_out)
        subprocess.run(cmd, check=True)
    elif fmt == "fasta":
        R1_out = os.path.join(outdir, "clean_R1.fasta.gz")
        R2_out = os.path.join(outdir, "clean_R2.fasta.gz")
        # Only create R2 if there are two reads (paired-end), otherwise set R2 to None
        os.makedirs(os.path.dirname(R1_out), exist_ok=True)
        shutil.copy2(reads[0], R1_out)
        if len(reads) > 1:
            os.makedirs(os.path.dirname(R2_out), exist_ok=True)
            shutil.copy2(reads[1], R2_out)
        else:
            R2_out = None
        with open(json_out, "w") as f: f.write("{}")
        with open(html_out, "w") as f: f.write("<html><body>FASTA input — no fastp QC</body></html>")
    else:
        raise ValueError("Unsupported format")
    return R1_out, R2_out

def run_megahit(R1, R2, outdir):
    if os.path.exists(outdir):
        print(f"Assembly folder {outdir} exists — removing it")
        shutil.rmtree(outdir)

    # Only use -2 if R2 is a valid file
    if not R2 or not os.path.exists(R2) or R2 == R1:
        cmd = ["megahit", "-r", R1, "-o", outdir, "--min-contig-len", "1000"]
    else:
        cmd = ["megahit", "-1", R1, "-2", R2, "-o", outdir, "--min-contig-len", "1000"]

    print("Running MEGAHIT: " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    return os.path.join(outdir, "final.contigs.fa")


def run_prodigal(contigs, outdir):
    os.makedirs(outdir, exist_ok=True)
    faa = os.path.join(outdir, "proteins.faa")
    ffn = os.path.join(outdir, "genes.ffn")
    gff = os.path.join(outdir, "prodigal.gff")
    cmd = ["prodigal", "-i", contigs, "-a", faa, "-d", ffn, "-o", gff, "-p", "meta"]
    print("Running prodigal: " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    return faa

def run_hmmsearch(faa, hmm, outdir):
    os.makedirs(outdir, exist_ok=True)
    tbl = os.path.join(outdir, os.path.basename(hmm).replace(".hmm", ".tblout"))
    cmd = ["hmmsearch", "--cpu", "4", "--domE", "1e-5", "--tblout", tbl, hmm, faa]
    print("Running hmmsearch: " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    return tbl

def parse_hits(faa, tbls, outcsv):
    os.makedirs(os.path.dirname(outcsv), exist_ok=True)
    with open(outcsv, "w") as out:
        out.write("family\ttarget\tevalue\n")
        for fam, tbl in tbls.items():
            with open(tbl) as f:
                for line in f:
                    if line.startswith("#"): continue
                    parts = line.strip().split()
                    out.write(f"{fam}\t{parts[0]}\t{parts[4]}\n")

############################################
# Main loop (sequential)
############################################

if __name__ == "__main__":
  # Download the study files if not already available
  if not os.path.exists(STUDY_DIR):
    fetch_study_downloads(STUDY_ACCESSION)

  # Start processing metagenomic pipeline
  analyses = list_analyses()
  for analysis in analyses:
      print(f"Processing {analysis}")
      analysis_dir = os.path.join(STUDY_DIR, analysis)
      cleaned_dir = os.path.join(RESULTS, analysis, "cleaned")
      assembly_dir = os.path.join(RESULTS, analysis, "assembly")
      genes_dir = os.path.join(RESULTS, analysis, "genes")
      hmm_dir = os.path.join(RESULTS, analysis, "hmm")
      report_dir = os.path.join(RESULTS, analysis, "report")

      # Paired-end data is supported too, but MGnify appears to provide single-end merged reads
      reads = get_MGnfy_single_end_fasta_sample_reads(analysis_dir)
      print(f"Reads: {reads}")
      if( len(reads) > 2 ):
          raise ValueError(f"Getting {len(reads)} reads from {analysis_dir}, expecting single-end or paired-end reads.")

      # NOTE: Can skip all these steps if input is already contigs or proteins
      # Only needed for fastq input, will just copy fasta files
      R1, R2 = try_run_fastp(reads, cleaned_dir)
      contigs = run_megahit(R1, R2, assembly_dir)
      faa = run_prodigal(contigs, genes_dir)

      tbls = {}
      for hmm in HMMs:
          tbl = run_hmmsearch(faa, hmm, hmm_dir)
          tbls[os.path.basename(hmm).replace(".hmm","")] = tbl

      outcsv = os.path.join(report_dir, "hmm_hits.csv")
      parse_hits(faa, tbls, outcsv)
      print(f"Finished {analysis}, report: {outcsv}")
