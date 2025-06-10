import os
import pathlib
import sys
from pathlib import Path

import pytest
from barcode_validator.cli import BarcodeValidatorCLI


@pytest.fixture
def data_dir():
    """Fixture to provide the test data directory."""
    # This assumes tests are run from the project root
    return Path("tests/data")

@pytest.fixture
def ncbi_taxdump(data_dir):
    """Fixture to locate the NCBI taxonomy dump"""
    ncbi_taxdump_archive = data_dir / "taxdump.tar.gz"
    return ncbi_taxdump_archive

@pytest.fixture
def structval_ncbi_tsv(data_dir):
    """Fixture to locate BCDM/TSV test file"""
    structval_ncbi_tsv = data_dir / "structval_ncbi.tsv"
    return structval_ncbi_tsv

@pytest.fixture
def test_cli_ncbi(structval_ncbi_tsv, ncbi_taxdump):
    """Fixture to run the CLI command."""

    # Save the original environment
    original_argv = sys.argv.copy()
    original_dir = os.getcwd()

    try:

        # Change to root of the repo
        script_dir = pathlib.Path(__file__).parent.parent
        os.chdir(script_dir)

        # Replace sys.argv with our test arguments
        sys.argv = [
            "barcode_validator",  # Program name
            "--input_file", str(structval_ncbi_tsv),
            "--mode", "structural",
            "--marker", "COI-5P",
            "--exp_taxonomy_type", "ncbi",
            "--exp_taxonomy", str(ncbi_taxdump),
            "--output_format", "tsv",
            "--log_level", "INFO"
        ]
        cli = BarcodeValidatorCLI()
        yield cli

    finally:
        # Restore the original environment
        sys.argv = original_argv
        os.chdir(original_dir)

def test_structural_validation_cli_bold(test_cli_ncbi):
    dars = test_cli_ncbi.run()
    assert dars is not None, "Validation results should not be None"
    assert len(dars.results) == 4, "Validation results should have 4 records"
    assert len(dars.results[0].stop_codons) == 1, "Second record should have 1 early stop"
    assert dars.results[1].seq_length == 76, "First record should have 1 stop codon"
    assert dars.results[2].full_ambiguities == 7, "Third record should have 7 ambiguities"