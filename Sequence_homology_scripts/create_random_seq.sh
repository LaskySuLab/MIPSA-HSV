#!/bin/bash

human_peptide="$1"

if [ -z "$human_peptide" ]; then
  echo "Usage: bash create_random_seq.sh <human_peptide>"
  exit 1
fi

awk -v N=1000 -v seed=123 '
  BEGIN { srand(seed) }

  /^>/ { header=$0; next }

  {
    gsub(/\*/, "", $0)
    split($0, orig, "")

    for (n = 1; n <= N; n++) {
      for (i in orig) a[i] = orig[i]

      for (i = length($0); i > 1; i--) {
        j = int(rand() * i) + 1
        tmp = a[i]; a[i] = a[j]; a[j] = tmp
      }

      seq = ""
      for (i = 1; i <= length($0); i++) seq = seq a[i]

      print header "_shuffle_" n
      print seq "*"
    }
  }
' "./${human_peptide}_sequence.fasta" \
> "./${human_peptide}_sequence_randomized_1000.fasta"