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

Usage
-----

## Web usage

The easiest way to use this functionality is by submitting a file to validate as per the
instructions given in the [data](data) section. This will trigger a bot to run the
validation and post the results as tabular output and summary reporting connected to 
the upload. This will only work if your upload is connected to the BGE project.

## Command line usage

The tooling can also be deployed locally. Dependency management is done via conda and
so installation is fairly painless. However, the local deployment will also require a
BLAST database that is indexed with taxon IDs from the NCBI taxonomy. This is the standard
`nr` database that can be downloaded using tools that come with the NCBI BLAST+ release
(which is part of the conda environment). Once that is downloaded and the 
[config file](config/config.yml) is updated as needed, the tool can be used as follows:

To get help:

```
python barcode_validator -h
```

To validate a FASTA file:

```
python barcode_validator -f examples/mge.fa -c config/config.yml
```

## Programmatic usage

The components of the toolset are programmed as object-oriented Python classes. The
classes have a reasonable separation of concerns such that other applications can
be composed out of them. To get started with this, peruse the Python files in the
[barcode_validator](barcode_validator) folder and the unit tests in the [tests](tests)
folder.
