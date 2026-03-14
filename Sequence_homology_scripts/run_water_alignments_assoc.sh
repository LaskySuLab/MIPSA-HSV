#!/bin/bash

viral="$1"
human_peptide="$2"

if [ -z "$viral" ] || [ -z "$human_peptide" ]; then
  echo "Usage: bash run_water_alignments_assoc.sh <viral> <human_peptide>"
  exit 1
fi

query_fasta="${viral}_${human_peptide}_assoc_seq.fasta"
subject_fasta="${human_peptide}_sequence.fasta"
outdir="water_alignments_assoc"

mkdir -p "$outdir"

# Counter for naming
i=0

# Read sequences from the multi-FASTA query file
awk '
/^>/{
    if (seq) {
        print header > f
        print seq >> f
        close(f)
    }

    match($0, /^>[^|]*\|[^|]*\|/)
    header = substr($0, RSTART, RLENGTH)

    seq=""
    f=sprintf("%s/query_%04d.fasta", "'$outdir'", ++i)
    next
}

{ seq=seq $0 }

END{
    if (seq) {
        print header > f
        print seq >> f
    }
}
' "$query_fasta"
# Loop through extracted query FASTA files
for q in "$outdir"/query_*.fasta; do
  base=$(basename "$q" .fasta)
  water -asequence "$q" -bsequence "$subject_fasta" \
        -gapopen 10 -gapextend 0.5 \
        -outfile "$outdir/${base}_vs_subject.txt"
done