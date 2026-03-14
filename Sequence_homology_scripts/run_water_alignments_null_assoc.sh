#!/bin/bash

viral="$1"
human_peptide="$2"

if [ -z "$viral" ] || [ -z "$human_peptide" ]; then
  echo "Usage: bash run_water_alignments_null_assoc.sh <viral> <human_peptide>"
  exit 1
fi

query_fasta="./${viral}_${human_peptide}_assoc_seq.fasta"
subject_fasta="./${human_peptide}_sequence_randomized_1000.fasta"
outdir="./water_alignments_randomized_null_assoc"


mkdir -p "$outdir"

# Function: print "header<TAB>sequence" for each FASTA record
fasta_to_tsv() {
  awk '
    BEGIN { header=""; seq="" }
    /^>/ {
      if (header != "") {
        gsub(/[ \t\r\n]/, "", seq)
        print header "\t" seq
      }
      header=$0
      seq=""
      next
    }
    { seq = seq $0 }
    END {
      if (header != "") {
        gsub(/[ \t\r\n]/, "", seq)
        print header "\t" seq
      }
    }
  ' "$1"
}

q=0
while IFS=$'\t' read -r header seq; do
  q=$((q+1))

  qfile=$(mktemp)
  printf "%s\n%s\n" "$header" "$seq" > "$qfile"

  s=0
  while IFS=$'\t' read -r sheader sseq; do
    s=$((s+1))

    sfile=$(mktemp)
    printf "%s\n%s\n" "$sheader" "$sseq" > "$sfile"

    water \
      -asequence "$qfile" \
      -bsequence "$sfile" \
      -gapopen 10 \
      -gapextend 0.5 \
      -outfile "$outdir/query_${q}_vs_subject_${s}.txt" \
      < /dev/null

    rm -f "$sfile"
  done < <(fasta_to_tsv "$subject_fasta")

  rm -f "$qfile"
done < <(fasta_to_tsv "$query_fasta")
 