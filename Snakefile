# Snakemake pipeline: Screen metagenomes for ArsR/SmtB (PF01022) and MerR (PF00376) regulators
# Expected layout:
#  - Raw reads: downloads/<study>/<analysis>/*_R{1,2}.fasta.gz (or single-end *fastq.gz)
#  - HMMs: hmm/PF01022.hmm, hmm/PF00376.hmm
#  - Outputs under: results/<analysis>/...

import os
from glob import glob

STUDY_DIR = "downloads/MGYS00001589" 
RESULTS = "results"

# Discover analyses automatically
analyses = []
if os.path.exists(STUDY_DIR):
    for root, dirs, files in os.walk(STUDY_DIR):
        # detect FASTQ presence
        fq = [f for f in files if f.endswith((".fasta.gz"))]
        if fq:
            # analysis id = last directory name with fastqs
            analyses.append(root)

analyses = sorted(set(analyses))

HMMs = ["hmm/PF01022.hmm", "hmm/PF00376.hmm"]  # ArsR/SmtB, MerR

rule all:
    input:
        expand("{results}/{analysis}/report/hmm_hits.csv", results=RESULTS, analysis=[os.path.basename(a) for a in analyses])

def sample_reads(analysis_dir):
    """Return tuple (R1, R2 or None) or single-end list for the given analysis_dir."""
    r1 = sorted(glob(os.path.join(analysis_dir, "*_R1*.fasta.gz")) + glob(os.path.join(analysis_dir, "*_1.fasta.gz")))
    r2 = sorted(glob(os.path.join(analysis_dir, "*_R2*.fasta.gz")) + glob(os.path.join(analysis_dir, "*_2.fasta.gz")))
    if r1 and r2:
        return r1[0], r2[0]
    # single-end fallback
    se = sorted(glob(os.path.join(analysis_dir, "*.fasta.gz")))
    if se:
        return (se[0], None)
    raise ValueError(f"No FASTQ files found in {analysis_dir}")

def outdir(analysis):
    return os.path.join(RESULTS, os.path.basename(analysis))

############################################
# Preprocess step (optional fastp)
############################################

rule preprocess_reads:
    input:
        lambda wildcards: sample_reads(os.path.join(STUDY_DIR, wildcards.analysis))
    output:
        R1="results/{analysis}/cleaned/clean_R1.fq.gz",
        R2="results/{analysis}/cleaned/clean_R2.fq.gz",
        json="results/{analysis}/cleaned/fastp.json",
        html="results/{analysis}/cleaned/fastp.html"
    run:
        reads = input
        fmt = detect_format(reads[0])

        shell("mkdir -p {results}/{wildcards.analysis}/cleaned")

        if fmt == "fastq":
            if len(reads) > 1:
                shell(f"""
                    fastp -i {reads[0]} -I {reads[1]} \
                          -o {output.R1} -O {output.R2} \
                          -j {output.json} -h {output.html}
                """)
            else:
                shell(f"""
                    fastp -i {reads[0]} \
                          -o {output.R1} \
                          -j {output.json} -h {output.html}
                    ln -sf $(basename {output.R1}) {output.R2} || true
                """)
        elif fmt == "fasta":
            # Just symlink for consistency
            shell(f"ln -sf {reads[0]} {output.R1}")
            shell(f"ln -sf {reads[0]} {output.R2}")
            # Dummy reports
            shell(f"echo '{{}}' > {output.json}")
            shell(f"echo '<html><body>FASTA input — no fastp QC</body></html>' > {output.html}")

############################################
# MEGAHIT assembly
############################################

rule megahit:
    input:
        R1 = rules.preprocess_reads.output.R1,
        R2 = rules.preprocess_reads.output.R2
    output:
        contigs = "{results}/{analysis}/assembly/final.contigs.fa"
    threads: 8
    shell:
        r"""
        mkdir -p {outdir}
        if [ -s {input.R2} ] && [ "{input.R2}" != "{input.R1}" ]; then
            megahit -1 {input.R1} -2 {input.R2} -o {outdir} -t {threads} --min-contig-len 1000
        else
            megahit -r {input.R1} -o {outdir} -t {threads} --min-contig-len 1000
        fi
        """.replace("{outdir}", r"{results}/{analysis}/assembly")

rule prodigal:
    input:
        contigs = rules.megahit.output.contigs
    output:
        faa = "{results}/{analysis}/genes/proteins.faa",
        ffn = "{results}/{analysis}/genes/genes.ffn",
        gff = "{results}/{analysis}/genes/prodigal.gff"
    shell:
        r"""
        mkdir -p {outdir}
        prodigal -i {input.contigs} -a {output.faa} -d {output.ffn} -o {output.gff} -p meta
        """.replace("{outdir}", r"{results}/{analysis}/genes")

rule hmmsearch:
    input:
        faa = rules.prodigal.output.faa,
        hmm = "hmm/{family}.hmm"
    output:
        tbl = "{results}/{analysis}/hmm/{family}.tblout"
    params:
        domE = "1e-5"
    shell:
        r"""
        mkdir -p {outdir}
        hmmsearch --cpu 4 --domE {params.domE} --tblout {output.tbl} {input.hmm} {input.faa} > /dev/null
        """.replace("{outdir}", r"{results}/{analysis}/hmm")

rule parse_hits:
    input:
        PF01022 = "{results}/{analysis}/hmm/PF01022.tblout",
        PF00376 = "{results}/{analysis}/hmm/PF00376.tblout",
        faa = "{results}/{analysis}/genes/proteins.faa"
    output:
        csv = "{results}/{analysis}/report/hmm_hits.csv"
    shell:
        r"""
        mkdir -p {outdir}
        python scripts/parse_hmm_tblout.py --faa {input.faa} \
            --tbl PF01022:{input.PF01022} PF00376:{input.PF00376} \
            --out {output.csv}
        """.replace("{outdir}", r"{results}/{analysis}/report")