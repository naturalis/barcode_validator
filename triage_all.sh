#!/usr/bin/env bash

# see: https://github.com/naturalis/barcode_validator/commit/0c9dfe6252449409bd1524817d9a8e24c7fa5af0
rm data/pass/concat_0c9dfe6.fasta
for f in $(git show --name-only 0c9dfe6 | grep data); do
 if [[ $f == *.fasta.tsv ]]; then
   STEM=$(basename $f .fasta.tsv)
   python barcode_validator/triage_sweep.py \
     -i data/${STEM}.fasta \
     -v $f \
     -o data/pass/${STEM}.fasta \
     -f data/fail/${STEM}.tsv 2> data/fail/${STEM}.log
   cat data/pass/${STEM}.fasta | sed -e 's/~/-/g' >> data/pass/concat_0c9dfe6.fasta
 fi
done
perl data/fasta-splitter.pl \
  --part-size 2000 \
  --measure count \
  --line-length 0 \
  --out-dir pass \
  data/pass/concat_0c9dfe6.fasta