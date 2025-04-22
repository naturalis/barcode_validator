#!/usr/bin/env python
from barcode_validator.cli import BarcodeValidatorCLI

"""
barcode_validator

A command line tool for DNA Barcode Validation.

This tool ingests FASTA/TSV sequence data, optionally merges CSV analytics and YAML configuration
into each result, performs structural and/or taxonomic validation, and outputs results as TSV.
Optionally, it can emit only valid sequences as FASTA.
"""

def main():
    cli = BarcodeValidatorCLI()
    cli.run()

if __name__ == "__main__":
    main()
