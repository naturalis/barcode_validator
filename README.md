# barcode-validator
## Standardized validation of barcodes

[BOLDigger-commandline](https://github.com/DominikBuchner/BOLDigger-commandline) is a Python program
that can check a FASTA file against the BOLD sequence ID service, producing a spreadsheet of near hits.
This can be used to validate barcodes produced, for example, by ARISE ONT barcoding, by BGE PacBio 
fresh barcoding, and by BGE genome skimming.

This repo contains a wrapper script with defaults to check the input file, understood to be COI. 
Functionally, the next thing it ought to do is:

1. parse the output excel sheet (e.g. with pandas)
2. parse an input sample sheet (which format / which columns?)
3. decide which taxonomic level to compare (genus / family)
4. produce a report
