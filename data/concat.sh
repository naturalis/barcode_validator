#!/usr/bin/env bash

# This will become a FASTA file that has all mge_.*_nocontam.fasta in it.
# The challenge here is that there are both *_standard_* files and *_fastp_* files,
# and the definition lines don't yet have that information so there will be ID clashes
# unless we suffix that qualifier.
OUTFASTA=concatenated_untriaged.fasta

# This will become a TSV file that has all mge_.*_nocontam.tsv in it.
# The challenge here is that there are both *_standard_* files and *_fastp_* files,
# and the sequence_id columns don't yet have that information so there will be ID clashes
# unless we suffix that qualifier. Furthermore, the columns aren't always in the same
# order, and not all files have all the columns, at least insofar that some assembly
# runs apparently produced additional analytics (in csv) while others didn't. This
# needs to be aligned.
OUTTSV=concatenated_untriaged.tsv

# Clean up previous, no doubt failed, version
rm $OUTFASTA
rm $OUTTSV

# First concatenate the fastp FASTA files, updating the defline in the stream
cat mge_fastp_*_nocontam.fasta | sed '/^>/s/$/_fastp/' >> $OUTFASTA

# Now do the same for the standard FASTA files
cat mge_standard_*_nocontam.fasta | sed '/^>/s/$/_standard/' >> $OUTFASTA

# Now merge and update the TSV files
python concat_tsv.py -i ./ -o $OUTTSV