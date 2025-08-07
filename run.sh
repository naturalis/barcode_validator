#!/usr/bin/env bash

# Usage: ./run.sh <input_file> [<output_fasta> <output_tsv>]
INPUT_FILE=$1
OUTPUT_FASTA=${2:-out.fasta}
OUTPUT_TSV=${3:-out.tsv}
BOLD_EXCEL=tests/data/bold.xlsx

python barcode_validator \
  --input-file $INPUT_FILE \
  # --csv-file $CSV_FILE \ # if you have CSV output from MGE it will join its lines with the input by sequence ID
  # --yaml-file $YAML_FILE \ # if you have batch-level metadata in YAML format it will paste it across all output lines
  --mode both \
  --marker COI-5P \
  --exp-taxonomy-type bold \
  --exp-taxonomy $BOLD_EXCEL \
  # --reflib-taxonomy-type bold \   # this is for when you have a local BLAST DB that needs to be mapped back to taxonomy tree, i.e. not needed for BOLD web service actually
  # --reflib-taxonomy $BOLD_EXCEL \ # same, see one line above
  --output-fasta $OUTPUT_FASTA \
  --output-tsv $OUTPUT_TSV \
  --taxon-validation method=bold \
  --taxon-validation rank=family \
  # --taxon-validation extent=class \ # this is for when you have a local BLAST DB that is taxonomically aware so that you can constrain the search space. Not needed for BOLD
  --triage-config group_id_separator=_ \
  --triage-config group_by_sample=true \
  --log-level ERROR 2> out.log
