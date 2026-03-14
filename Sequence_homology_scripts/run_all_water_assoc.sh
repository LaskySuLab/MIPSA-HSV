#!/bin/bash
#SBATCH --job-name=Sequence_homology
#SBATCH --output=./Sequence_homology%j.out
#SBATCH --error=./Sequence_homology%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4

viral="$1"
human_peptide="$2"

if [ -z "$viral" ] || [ -z "$human_peptide" ]; then
  echo "Usage: sbatch run_all_water_assoc.sh <viral> <human_peptide>"
  exit 1
fi

# conda activate MIPSA_homology
bash ./run_water_alignments_assoc.sh "$viral" "$human_peptide" &
bash ./create_random_seq.sh "$human_peptide" &
bash ./run_water_alignments_null_assoc.sh "$viral" "$human_peptide"

wait
