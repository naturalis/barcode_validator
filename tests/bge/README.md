These tests are to ensure proper functioning of the application for BGE. This means the following:
- input data is [FASTA](data/five_test_processids.fasta) where the first word is the process ID + '_' + assembly attempt ID
- ancillary input consists of [CSV](data/five_test_processids-stats.csv) files and YAML files
- input taxonomy is a BOLD spreadsheet with 'Lab Sheet' and 'Taxonomy' tabs in [XLSX](data/bold.xlsx) format
- structural validation is performed, initially for COI, and soon for other markers
- taxonomic validation is performed with the BOLD web service
- triage is performed with grouping by assembly attempt IDs

The invocation of the tool as applied to this use case is as follows:

```bash

python barcode_validator \
  --input-file <input_fasta> \
  --csv-file <input_csv> \
  --yaml-file <input_yaml> \
  --mode both \
  --marker COI-5P \
  --input-resolver format=bold \
  --input-resolver file=<bold_excel> \
  --output-fasta <output_fasta> \
  --output-tsv <output_tsv> \
  --taxon-validation method=bold \
  --taxon-validation rank=family \
  --taxon-validation min_identity=0.8 \
  --taxon-validation max_target_seqs=100 \
  --triage-config group_id_separator=_ \
  --triage-config group_by_sample=true \
  --log-level ERROR 2> <log_file>
```

Questions:

# Are the ancillary input files required?

The CSV and YAML files can be omitted but it is strongly recommended that they are there for consistency's sake.

# Is the BOLD spreadsheet required?

Yes, and it must be up to date.  Attempting to validate specimens that have not yet been registered will fail.
A theoretical fallback option is to Frankenstein a BOLD spreadsheet like file with your own 'Lab Sheet' and
'Taxonomy' tabs with the intended lineages. Thoughts and prayers for those that attempt this manually.

# What's the roadmap for other markers?

For protein-coding markers (e.g. matK, rbcL), additional HMMs need to be produced and added to the hmm_files folder,
with the naming convention matK.hmm, rbcL.hmm, and so on.

For non-coding markers, no HMMs are needed because we only use them to compute the amino acid translation to find
nonsense stop codons.

For both coding and non-coding markers, the validation criteria will need to be refined in criteria.py. There are
placeholder settings for sequence length, number of ambiguities, and so on. Whether these are appropriate is to
be determined. This will mean sticking some numbers in criteria.py in the right place so it is fairly painless.

# How about other taxon validation services?

Local BLAST is hard to get right and make it perform well enough. It works but since BOLD works better in basically
all aspects we don't plan to spend time on it.

Remote BLAST is an option that may be fairly easy to implement because the BLASTN wrapper from nbitk can be set to
remote, but this is untested and doesn't seem needed, for reasons stated above.

# How is triaging done?

Within each group of assembly attempts, we simply pick the longest valid sequence. Arguably there are more 
sophisticated criteria but how to rank them (e.g. is the next longest one with fewer ambiguities better) is
not all that straightforward.