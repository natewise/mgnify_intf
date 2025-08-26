# Create a small, ready-to-run metagenomic screening pipeline scaffold.
# It includes:
# - environment.yml (conda env with common tools)
# - Snakefile (QC -> assembly -> gene calling -> HMM search -> report)
# - scripts/parse_hmm_tblout.py (parses HMMER --tblout to CSV)
# - README.md with instructions
# The pipeline expects FASTQ files under downloads/<study>/<analysis>/*
# and HMM profiles placed in hmm/PF01022.hmm (ArsR/SmtB family) and hmm/PF00376.hmm (MerR family).

import os, textwrap, json, tarfile, io, pathlib

# Get base dir from location of this script
base_dir = pathlib.Path(__file__).parent.resolve()
scripts_dir = f"{base_dir}/scripts"
hmm_dir = f"{base_dir}/hmm"

#!/usr/bin/env python3
import argparse, csv, re
from collections import defaultdict

def read_faa_lengths(faa_path):
    lengths = {}
    with open(faa_path) as fh:
        seq_id = None
        seq = []
        for line in fh:
            if line.startswith(">"):
                if seq_id is not None:
                    lengths[seq_id] = len("".join(seq))
                seq_id = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if seq_id is not None:
            lengths[seq_id] = len("".join(seq))
    return lengths

def parse_tblout(tbl_path):
    hits = []
    with open(tbl_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            # HMMER3 --tblout columns (first 18 are fixed)
            parts = line.strip().split()
            target_name = parts[0]
            query_name  = parts[2]
            evalue      = float(parts[4])
            score       = float(parts[5])
            bias        = float(parts[6])
            hmm_from    = int(parts[15])
            hmm_to      = int(parts[16])
            ali_from    = int(parts[17])
            ali_to      = int(parts[18])
            hits.append(dict(target=target_name, query=query_name, evalue=evalue,
                             score=score, bias=bias, hmm_from=hmm_from, hmm_to=hmm_to,
                             ali_from=ali_from, ali_to=ali_to))
    return hits

def main():
    ap = argparse.ArgumentParser(description="Parse HMMER --tblout files and summarize hits")
    ap.add_argument("--faa", required=True, help="Proteins FASTA used as target (from Prodigal)")
    ap.add_argument("--tbl", nargs="+", required=True, help="List like PF01022:path PF00376:path")
    ap.add_argument("--out", required=True, help="Output CSV")
    args = ap.parse_args()

    lengths = read_faa_lengths(args.faa)

    all_hits = []
    for item in args.tbl:
        fam, path = item.split(":", 1)
        for h in parse_tblout(path):
            h["family"] = fam
            h["prot_len"] = lengths.get(h["target"], "")
            # simple motif flags (optional heuristics)
            # ArsR/SmtB often has Cys pairs for metal binding (Cys-x(n)-Cys)
            seqid = h["target"]
            # motif extraction skipped here to avoid reading sequences; could add if needed
            all_hits.append(h)

    # sort by evalue
    all_hits.sort(key=lambda x: (x["family"], x["evalue"], -x["score"]))

    cols = ["family","target","prot_len","evalue","score","bias","ali_from","ali_to","hmm_from","hmm_to"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for h in all_hits:
            w.writerow({k: h.get(k, "") for k in cols})

if __name__ == "__main__":
    main()
