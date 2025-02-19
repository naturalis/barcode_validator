#!/usr/bin/env bash

# Operating on commit hash 0c9dfe6,
# see: https://github.com/naturalis/barcode_validator/commit/0c9dfe6252449409bd1524817d9a8e24c7fa5af0
# Remove the previous concatenation
rm data/pass/concat_0c9dfe6.fasta

# Iterate over the files in commit hash 0c9dfe6 and filter on those with data in the file path
for f in $(git show --name-only 0c9dfe6 | grep data); do

 # Focal file is a TSV file
 if [[ $f == *.fasta.tsv ]]; then

   # Extract the basename without containing path or extensions
   STEM=$(basename $f .fasta.tsv)

   # Invoke the triage tool, forward STDERR to a log file
   python barcode_validator/triage_sweep.py \
     -i data/${STEM}.fasta \
     -v $f \
     -o data/pass/${STEM}.fasta \
     -f data/fail/${STEM}.tsv 2> data/fail/${STEM}.log

   # Append the passing sequences to the concatenation, replace tildas with spaces
   cat data/pass/${STEM}.fasta | sed -e 's/~/-/g' >> data/pass/concat_0c9dfe6.fasta
 fi
done

# Chunk the concatenated file into batches of 2000 sequences
perl data/fasta-splitter.pl \
  --part-size 2000 \
  --measure count \
  --line-length 0 \
  --out-dir pass \
  data/pass/concat_0c9dfe6.fasta