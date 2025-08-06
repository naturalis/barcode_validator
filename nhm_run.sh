#!/usr/bin/env bash
INPUT_FILE=tests/data/dan_test.fasta
BOLD_EXCEL=tests/data/bold.xlsx

python barcode_validator \
  --input-file $INPUT_FILE \
  --mode both \
  --marker COI-5P \
  --exp-taxonomy-type bold \
  --exp-taxonomy $BOLD_EXCEL \
  --reflib-taxonomy-type bold \
  --reflib-taxonomy $BOLD_EXCEL \
  --output-format tsv \
  --taxon-validation method=bold \
  --taxon-validation rank=family \
  --taxon-validation extent=class \
  --triage-config group_id_separator=_ \
  --triage-config group_by_sample=true \
  --log-level ERROR > out.tsv