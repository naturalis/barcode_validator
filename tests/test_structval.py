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
def config_file(data_dir):
    """Fixture to provide a test configuration file."""
    config_path = data_dir / "config_structural.yml"
    return config_path


@pytest.fixture
def input_fasta(data_dir):
    """Fixture to create a test FASTA file with sequences."""
    fasta_path = data_dir / "structval_early_stop.fasta"
    return fasta_path


@pytest.fixture
def bold_excel(data_dir):
    """Fixture to provide a test BOLD Excel file."""
    bold_path = data_dir / "bold.xlsx"
    return bold_path


@pytest.fixture
def test_cli(config_file, input_fasta, bold_excel):
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
            "--input_file", str(input_fasta),
            "--config", str(config_file),
            "--mode", "structural",
            "--exp_taxonomy_type", "bold",
            "--exp_taxonomy", str(bold_excel),
            "--output_format", "tsv",
            "--log_level", "INFO"
        ]
        cli = BarcodeValidatorCLI()
        yield cli

    finally:
        # Restore the original environment
        sys.argv = original_argv
        os.chdir(original_dir)

def test_structural_validation_cli(test_cli):
    dars = test_cli.run()
    assert dars is not None, "Validation results should not be None"
    assert len(dars.results) == 4, "Validation results should have 4 records"
    assert len(dars.results[0].stop_codons) == 1, "Second record should have 1 early stop"
    assert dars.results[1].seq_length == 76, "First record should have 1 stop codon"
    assert dars.results[2].full_ambiguities == 7, "Third record should have 7 ambiguities"