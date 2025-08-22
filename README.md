# Metagenomic Biosensor Screening Pipeline

This Snakemake pipeline screens metagenomic datasets for **metal-sensing transcription factor families** used in biosensors like those in Lux et al. (2021):
- **ArsR/SmtB family** (PF01022) — includes many As/Cd/Zn sensors (e.g., ArsR, CadC, SmtB).
- **MerR family** (PF00376) — includes Hg sensors (MerR-like regulators).

## What it does
1. **QC:** `fastp` on FASTQ reads.
2. **Assembly:** `megahit` (>=1 kb contigs).
3. **Gene calling:** `prodigal -p meta` to predict proteins.
4. **HMM search:** `hmmsearch` against PF01022 and PF00376.
5. **Report:** CSV of hits with scores/e-values per analysis.

## Quick start

```bash
# 1) Create environment
mamba env create -f environment.yml
mamba activate meta-biosensor

# 2) Place your reads under downloads/<study>/<analysis>/*.fastq.gz
#    (You already have MGnify downloads from your Python script.)

# 3) Add HMM profiles:
#    - Put PF01022.hmm (ArsR/SmtB) into downloads/hmm/PF01022.hmm
#    - Put PF00376.hmm (MerR)     into downloads/hmm/PF00376.hmm
#    (Use Pfam-A or curated HMMs; run `hmmpress` on them once.)

# 4) Run the pipeline
snakemake -j 8
