import os
import subprocess
import shutil

STUDY_DIR = "downloads/MGYS00001589"
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

def sample_reads(analysis_dir):
    files = sorted([os.path.join(analysis_dir,f) for f in os.listdir(analysis_dir) if f.endswith(".gz")])
    if not files:
        raise ValueError(f"No .gz files found in {analysis_dir}")
    return files

############################################
# Pipeline steps
############################################

def run_fastp(reads, outdir):
    os.makedirs(outdir, exist_ok=True)
    fmt = detect_format(reads[0])
    R1_out = os.path.join(outdir, "clean_R1.fq.gz")
    R2_out = os.path.join(outdir, "clean_R2.fq.gz")
    json_out = os.path.join(outdir, "fastp.json")
    html_out = os.path.join(outdir, "fastp.html")

    if fmt == "fastq":
        if len(reads) > 1:
            cmd = ["fastp", "-i", reads[0], "-I", reads[1], "-o", R1_out, "-O", R2_out, "-j", json_out, "-h", html_out]
        else:
            cmd = ["fastp", "-i", reads[0], "-o", R1_out, "-j", json_out, "-h", html_out]
            if not os.path.exists(R2_out):
                os.symlink(R1_out, R2_out)
        subprocess.run(cmd, check=True)
    elif fmt == "fasta":
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

    subprocess.run(cmd, check=True)
    return os.path.join(outdir, "final.contigs.fa")


def run_prodigal(contigs, outdir):
    os.makedirs(outdir, exist_ok=True)
    faa = os.path.join(outdir, "proteins.faa")
    ffn = os.path.join(outdir, "genes.ffn")
    gff = os.path.join(outdir, "prodigal.gff")
    cmd = ["prodigal", "-i", contigs, "-a", faa, "-d", ffn, "-o", gff, "-p", "meta"]
    subprocess.run(cmd, check=True)
    return faa

def run_hmmsearch(faa, hmm, outdir):
    os.makedirs(outdir, exist_ok=True)
    tbl = os.path.join(outdir, os.path.basename(hmm).replace(".hmm", ".tblout"))
    cmd = ["hmmsearch", "--cpu", "4", "--domE", "1e-5", "--tblout", tbl, hmm, faa]
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

analyses = list_analyses()
for analysis in analyses:
    print(f"Processing {analysis}")
    analysis_dir = os.path.join(STUDY_DIR, analysis)
    cleaned_dir = os.path.join(RESULTS, analysis, "cleaned")
    assembly_dir = os.path.join(RESULTS, analysis, "assembly")
    genes_dir = os.path.join(RESULTS, analysis, "genes")
    hmm_dir = os.path.join(RESULTS, analysis, "hmm")
    report_dir = os.path.join(RESULTS, analysis, "report")

    reads = sample_reads(analysis_dir)
    R1, R2 = run_fastp(reads, cleaned_dir)
    contigs = run_megahit(R1, R2, assembly_dir)
    faa = run_prodigal(contigs, genes_dir)

    tbls = {}
    for hmm in HMMs:
        tbl = run_hmmsearch(faa, hmm, hmm_dir)
        tbls[os.path.basename(hmm).replace(".hmm","")] = tbl

    outcsv = os.path.join(report_dir, "hmm_hits.csv")
    parse_hits(faa, tbls, outcsv)
    print(f"Finished {analysis}, report: {outcsv}")
