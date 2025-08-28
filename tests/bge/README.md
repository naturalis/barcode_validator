# BGE Barcode Validator Tests

These tests are to ensure proper functioning of the application for BGE. This means the following:

- Input data is [FASTA](data/WK-3860_BSNHM190-con_seqs.fasta) where the first word is the process ID + '_' + assembly 
  attempt ID
- Ancillary input consists of a [CSV](data/WK-3860_BSNHM190-metrics.csv) file where the ID column is the process ID
- Input taxonomy is a BOLD spreadsheet with 'Lab Sheet' and 'Taxonomy' tabs in 
  [XLSX](data/bold_bge_container_plus-updated.xlsx) format
- Structural validation is performed and the best, valid sequence for each process ID is retained ('triage')
- Taxonomic validation is performed just on the triaged results, using the Galaxy BLAST web service

The invocation of the tool as applied to this use case is as follows:

## Set up the execution environment

In this step we set up the conda/mamba environment and the PYTHONPATH so that the module can be invoked. The exact
commands may differ depending on your setup. For example, if you use conda instead of mamba, or if you have a different
way of managing environments and PYTHONPATH. Here we assume that we are in the same folder as this README.md file.

```bash
mamba activate barcode_validator
export PYTHONPATH=$PYTHONPATH:$(pwd)/../../
```
## Specify input and output files for structural validation

Here we specify the input file from MGEE, as well is its analytics. We also specify the BOLD spreadsheet with the
intended taxonomic lineages. Finally, we specify the output files for structural validation, including the log file.

```bash
INPUT_FASTA=data/WK-3860_BSNHM190-con_seqs.fasta
INPUT_CSV=data/WK-3860_BSNHM190-metrics.csv
BOLD_EXCEL=data/bold_bge_container_plus-updated.xlsx
STRUCTVAL_FASTA=data/structval_out.fasta
STRUCTVAL_TSV=data/structval_out.tsv
STRUCTVAL_LOG=data/structval.log
```

## Execute structural validation

Here we invoke the structural validation step, which also performs triaging to retain the best valid sequence. We run
the tool as a module so that we don't have to worry about the PATH. However, this assumes that the PYTHONPATH is set
up correctly as shown above. In this run we specify the input files, the marker (COI-5P), the output files, and
triaging options. The log level is set to INFO to get a reasonable amount of detail in the log file. We resolve sequence
IDs to taxa using the BOLD spreadsheet. We need to do this here because we need to know the higher taxa to infer the
genetic code for translation.

```bash
python -m barcode_validator \
  --input-file $INPUT_FASTA \
  --csv-file $INPUT_CSV \
  --mode structural \
  --marker COI-5P \
  --input-resolver format=bold \
  --input-resolver file=$BOLD_EXCEL \
  --output-fasta $STRUCTVAL_FASTA \
  --output-tsv $STRUCTVAL_TSV \
  --triage-config group_id_separator=_ \
  --triage-config group_by_sample=true \
  --log-level INFO 2> $STRUCTVAL_LOG
```

## Set up Galaxy parameters

In the next run we will use the Galaxy BLAST web service for taxonomic validation. This requires setting up the Galaxy
API key and domain. The API key has been redacted here for security reasons. You will need to set this to your own key.

```bash
export GALAXY_API_KEY=******************************** # set your Galaxy API key
export GALAXY_DOMAIN=galaxy.naturalis.nl
```

## Specify output files for taxonomic validation

The input FASTA is now the output from the structural validation step. Note that we are continuing to use the 
$INPUT_CSV and $BOLD_EXCEL files from above. So here we just specify the output files for taxonomic validation, 
including the log file, and assume that the other locations (CSV, BOLD sheet) are still defined in our environment.

```bash
TAXVAL_FASTA=data/taxval_out.fasta
TAXVAL_TSV=data/taxval_out.tsv
TAXVAL_LOG=data/taxval.log
```
## Execute taxonomic validation

Here we invoke the taxonomic validation step, which uses the Galaxy BLAST web service. We run
the tool as a module so that we don't have to worry about the PATH. However, this assumes that the PYTHONPATH is set
up correctly as shown above. In this run we specify the input files, the marker (COI-5P), the output files, and
the taxonomic validation options. The log level is set to INFO to get a reasonable amount of detail in the log file. 
We resolve sequence IDs to taxa using the BOLD spreadsheet. We need to do this here because we check the observed 
taxonomic lineages against the intended ones from the BOLD sheet.

```bash
python -m barcode_validator \
  --input-file $STRUCTVAL_FASTA \
  --csv-file $INPUT_CSV \
  --mode taxonomic \
  --marker COI-5P \
  --input-resolver format=bold \
  --input-resolver file=$BOLD_EXCEL \
  --output-fasta $TAXVAL_FASTA \
  --output-tsv $TAXVAL_TSV \
  --taxon-validation method=galaxy \
  --taxon-validation rank=family \
  --taxon-validation min_identity=0.8 \
  --taxon-validation max_target_seqs=100 \
  --log-level INFO 2> $TAXVAL_LOG
```

This produces a [FASTA](data/taxval_out.fasta) file with the best, structurally and taxonomically valid sequence for 
each process ID, and a [TSV](data/taxval_out.tsv) file with the taxon validation results.

# Questions:

## Are the ancillary input files required?

The CSV file can be omitted but it is strongly recommended that it is there for consistency's sake.

## Is the BOLD spreadsheet required?

Yes, and it must be up to date.  Attempting to validate specimens that have not yet been registered will fail.
A theoretical fallback option is to Frankenstein a BOLD spreadsheet-like file with your own 'Lab Sheet' and
'Taxonomy' tabs with the intended lineages. Thoughts and prayers for those that attempt this manually.

## What's the roadmap for other markers?

For protein-coding markers (e.g. matK, rbcL), additional HMMs need to be produced and added to the hmm_files folder,
with the naming convention matK.hmm, rbcL.hmm, and so on.

For non-coding markers, no HMMs are needed because we only use them to compute the amino acid translation to find
nonsense stop codons.

For both coding and non-coding markers, the validation criteria will need to be refined in criteria.py. There are
placeholder settings for sequence length, number of ambiguities, and so on. Whether these are appropriate is to
be determined. This will mean sticking some numbers in criteria.py in the right place so it is fairly painless.

## How about other taxon validation services?

Local BLAST is hard to get right and make it perform well enough. It works but since BOLD works better in basically
all aspects we don't plan to spend time on it.

Remote BLAST is an option that may be fairly easy to implement because the BLASTN wrapper from nbitk can be set to
remote, but this is untested and doesn't seem needed, for reasons stated above.

## How is triaging done?

Within each group of assembly attempts, we simply pick the longest valid sequence. Arguably there are more 
sophisticated criteria but how to rank them (e.g. is the next longest one with fewer ambiguities better) is
not all that straightforward.