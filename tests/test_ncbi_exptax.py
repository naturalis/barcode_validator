import os
import pathlib
import sys
import tempfile
from pathlib import Path

import pytest
import requests

from barcode_validator.cli import BarcodeValidatorCLI


@pytest.fixture
def data_dir():
    """Fixture to provide the test data directory."""
    # This assumes tests are run from the project root
    return Path("tests/data")

@pytest.fixture
def ncbi_taxdump():
    """Fixture to locate the NCBI taxonomy dump"""
    tar_path = None

    # URL of the NCBI taxonomy dump
    url = "http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

    # Create a temporary file
    with tempfile.NamedTemporaryFile(suffix=".tar.gz", delete=False) as temp_file:
        # Download the file
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an exception for HTTP errors

        # Save the downloaded file
        for chunk in response.iter_content(chunk_size=8192):
            temp_file.write(chunk)

        # Get the path of the temporary file
        tar_path = temp_file.name

    return tar_path

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
            "--input-file", str(structval_ncbi_tsv),
            "--mode", "structural",
            "--marker", "COI-5P",
            "--exp-taxonomy-type", "ncbi",
            "--exp-taxonomy", str(ncbi_taxdump),
            "--output-format", "tsv",
            "--log-level", "INFO"
        ]
        cli = BarcodeValidatorCLI()
        yield cli

    finally:
        # Restore the original environment
        sys.argv = original_argv
        os.chdir(original_dir)

def test_structural_validation_cli_ncbi(test_cli_ncbi):
    dars = test_cli_ncbi.run()
    assert dars is not None, "Validation results should not be None"
    assert len(dars.results) == 4, "Validation results should have 4 records"
    assert len(dars.results[0].stop_codons) == 0, "Second record should have 0 early stop"
    assert dars.results[1].seq_length == 76, "First record is 76bp"
    assert dars.results[2].full_ambiguities == 7, "Third record should have 7 ambiguities"