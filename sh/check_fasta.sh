#!/usr/bin/env bash
# check_fasta.sh
# Usage: ./check_fasta.sh <input.fasta|fastq[.gz]>
# Example: ./sh/check_fasta.sh downloads/MGYS00001589/MGYA00103632/ERR788946_MERGED_FASTQ.fasta

IN="$1"
if [[ -z "$IN" ]]; then
  echo "Usage: $0 <input.fasta|fastq[.gz]>"
  exit 1
fi

# detect compression
if [[ "$IN" == *.gz ]]; then
  CAT="zcat"
else
  CAT="cat"
fi

echo "Checking file: $IN"

# detect format (FASTA vs FASTQ) by first char
firstchar=$($CAT "$IN" | head -n1 | cut -c1)
if [[ "$firstchar" == ">" ]]; then
  fmt="FASTA"
elif [[ "$firstchar" == "@" ]]; then
  fmt="FASTQ"
else
  echo "Could not detect file format (expected > or @ on first line)"
  exit 1
fi
echo "Detected format: $fmt"

# get lengths depending on format
if [[ "$fmt" == "FASTA" ]]; then
  lengths=$($CAT "$IN" \
    | awk 'BEGIN{seq=""} /^>/{if(seq!=""){print length(seq)} seq=""; next} {seq=seq$0} END{if(seq!=""){print length(seq)}}')
else
  lengths=$($CAT "$IN" \
    | awk 'NR%4==2 {print length($0)}')
fi

# compute stats
echo "$lengths" \
  | awk 'NR==1{min=max=$1} {if($1<min)min=$1; if($1>max)max=$1; sum+=$1; n++} 
         END{print "n_reads:",n; print "min_len:",min; print "max_len:",max; print "mean_len:",sum/n}'

# classify heuristically
mean=$(echo "$lengths" | awk '{sum+=$1; n++} END{if(n>0) print sum/n; else print 0}')
max=$(echo "$lengths" | awk 'max<$1{max=$1} END{print max+0}')

if (( $(echo "$mean < 400" | bc -l) )) && (( max < 1000 )); then
  echo "→ Looks like amplicon/16S-style short reads"
else
  echo "→ Looks like shotgun/metagenome-style reads"
fi
