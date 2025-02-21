#!/usr/bin/env bash

# Operating on commit hash 0c9dfe6,
# see: https://github.com/naturalis/barcode_validator/commit/0c9dfe6252449409bd1524817d9a8e24c7fa5af0
# Clear previous output
rm data/pass/* data/fail/*

# Iterate over the files in the data dir to select fasta files
for f in $(ls data/*MGE-BGE_r1_1.3_1.5_s50_100.fasta*); do
 echo "$f"
 # Focal file is a TSV file (check if fasta has a tsv)
 if [[ $f == *.tsv ]]; then

   # Extract the basename without containing path or extensions
   STEM=$(basename $f .fasta.tsv)

   # Invoke the triage tool, forward STDERR to a log file
   python barcode_validator/triage-sweep-local.py \
     -i data/${STEM}.fasta \
     -v $f \
     -o data/pass/${STEM}.fasta \
     -f data/fail/${STEM}.tsv 2> data/fail/${STEM}.log

   # Append the passing sequences to the concatenation, replace tildas with spaces
   cat data/pass/${STEM}.fasta | sed -e 's/~/-/g' >> data/pass/concat.fasta
 fi
done

# Chunk the concatenated file into batches of 2000 sequences
perl data/fasta-splitter.pl \
  --part-size 2000 \
  --measure count \
  --line-length 0 \
  --out-dir data/pass \
  data/pass/concat.fasta