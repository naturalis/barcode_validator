import os
import pathlib
import subprocess
import pytest
from pathlib import Path
import pandas as pd
from Bio import SeqIO
import sys

from barcode_validator.cli import BarcodeValidatorCLI


@pytest.fixture
def data_dir():
    """Fixture to provide the test data directory."""
    # This assumes tests are run from the project root
    return Path("tests/data")


@pytest.fixture
def input_fasta_bold(data_dir):
    """Fixture to create a test FASTA file with sequences."""
    fasta_path = data_dir / "structval_early_stop.fasta"
    return str(fasta_path)


@pytest.fixture
def input_tsv_nsr(data_dir):
    """Fixture to create a test TSV file with sequences."""
    tsv_path = data_dir / "structval_nsr.tsv"
    return str(tsv_path)


@pytest.fixture
def bold_excel(data_dir):
    """Fixture to provide a test BOLD Excel file."""
    bold_path = data_dir / "bold.xlsx"
    return str(bold_path)


@pytest.fixture
def nsr_dwca(data_dir):
    """Fixture to provide a test NSR archive."""
    nsr_path = data_dir / "nsr-20250207.dwca.zip"
    return str(nsr_path)


@pytest.fixture
def test_cli_bold(input_fasta_bold, bold_excel):
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
            "--input-file", input_fasta_bold,
            "--mode", "structural",
            "--marker", "COI-5P",
            "--exp-taxonomy-type", "bold",
            "--exp-taxonomy", bold_excel,
            "--output-format", "tsv",
            "--log-level", "INFO"
        ]
        cli = BarcodeValidatorCLI()
        yield cli

    finally:
        # Restore the original environment
        sys.argv = original_argv
        os.chdir(original_dir)


@pytest.fixture
def test_cli_nsr(input_tsv_nsr, nsr_dwca):
    """Fixture to run the CLI command with NSR data."""

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
            "--input-file", input_tsv_nsr,
            "--mode", "structural",
            "--marker", "COI-5P",
            "--exp-taxonomy-type", "nsr",
            "--exp-taxonomy", nsr_dwca,
            "--output-format", "tsv",
            "--log-level", "INFO"
        ]
        cli = BarcodeValidatorCLI()
        yield cli

    finally:
        # Restore the original environment
        sys.argv = original_argv
        os.chdir(original_dir)

def test_structural_validation_cli_bold(test_cli_bold):
    dars = test_cli_bold.run()
    assert dars is not None, "Validation results should not be None"
    assert len(dars.results) == 4, "Validation results should have 4 records"
    assert len(dars.results[0].stop_codons) == 0, "Second record should have - early stop"
    assert dars.results[1].seq_length == 76, "First record is 76bp"
    assert dars.results[2].full_ambiguities == 7, "Third record should have 7 ambiguities"

def test_structural_validation_cli_nsr(test_cli_nsr):
    dars = test_cli_nsr.run()
    assert dars is not None, "Validation results should not be None"
    assert len(dars.results) == 4, "Validation results should have 4 records"
    assert len(dars.results[0].stop_codons) == 0, "Second record should have 0 early stop"
    assert dars.results[1].seq_length == 76, "First record s 76bp"
    assert dars.results[2].full_ambiguities == 7, "Third record should have 7 ambiguities"