Introduction
============

Sequences of the DNA barcode marker COI are produced in a variety of ways, including 
by genome skimming of old collection specimens, by nanopore sequencing on ONT platforms,
and by smrt sequencing on the PacBio Sequel platform. In all cases, the resulting 
sequence ought to be validated before publishing. This project aims to provide a 
universal solution for this, in order to meet the needs of the [BGE](https://biodiversitygenomics.eu/)
project and the [ARISE](https://www.arise-biodiversity.nl/) project. Broadly speaking,
what is provided here is a service that does the following:

1. Check to see if the expected taxon and that observed by checking the sequence against
   a reference database (GenBank or BOLD) match. The match can be performed at a higher
   taxonomic level, e.g. the family. Mismatches may indicate contaminations.
2. Check to see if there are no early stop codons in the amino acid translation. Early
   stop codons may indicate the erroneous amplification of a pseudogene copy, such as
   a NUMT sequence.
3. Persist the validation results somewhere. At present nothing is effectively done with
   this, but the hooks are there to persist to some data store, like a google sheet, a
   database, or a file. Now the output is simply printed to the console as TSV.

Quick start
-----------

```
python barcode_validator -h
```

For example:

```
python barcode_validator -f examples/mge.fa -c config/config.yml
```
