import os
import pathlib
import tempfile

import pytest
from pathlib import Path
import sys

import requests

from barcode_validator.cli import BarcodeValidatorCLI

# do I still work?

@pytest.fixture(autouse=True)
def check_galaxy_key():
    """Verify Galaxy API key is available"""
    if not os.environ.get('GALAXY_API_KEY'):
        pytest.skip("GALAXY_API_KEY not set in environment")

@pytest.fixture(autouse=True)
def check_galaxy_domain():
    """Verify Galaxy API key is available"""
    if not os.environ.get('GALAXY_DOMAIN'):
        pytest.skip("GALAXY_DOMAIN not set in environment")

@pytest.fixture
def data_dir():
    """Fixture to provide the test data directory."""
    # This assumes tests are run from the project root
    return Path("tests/galaxy/data")

@pytest.fixture
def input_fasta(data_dir):
    """Fixture to create a test FASTA file with sequences."""
    fasta_path = data_dir / "NLVRT.fas"
    return str(fasta_path)

@pytest.fixture
def bold_excel(data_dir):
    """Fixture to provide a test BOLD Excel file."""
    bold_path = data_dir / "NLVRT.xlsx"
    return str(bold_path)

@pytest.fixture
def taxdump():
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

    # Yield the opened tar file
    yield tar_path

    # Clean up
    os.unlink(tar_path)



@pytest.fixture
def cli_prepare_coi(input_fasta, bold_excel, data_dir, taxdump):
    """Fixture to run the CLI command."""

    # Save the original environment
    original_argv = sys.argv.copy()
    original_dir = os.getcwd()

    try:

        # Change to root of the repo
        script_dir = pathlib.Path(__file__).parent.parent.parent
        os.chdir(script_dir)

        # Replace sys.argv with our test arguments
        sys.argv = [
            "barcode_validator",  # Program name
            "--input-file", input_fasta,
            "--mode", "both",
            "--marker", "COI-5P",
            "--input-resolver", "format=bold",
            "--input-resolver", f"file={bold_excel}",
            "--reflib-resolver", "format=ncbi",
            "--reflib-resolver", f"file={taxdump}",
            "--output-fasta", f"{data_dir}/arise_out.fasta",
            "--output-tsv", f"{data_dir}/arise_out.tsv",
            "--taxon-validation", "method=galaxy",
            "--taxon-validation", "rank=family",
            "--taxon-validation", "min_identity=0.8",
            "--taxon-validation", "max_target_seqs=100",
            "--triage-config", "group_by_sample=false",
            "--log-level", "ERROR"
        ]
        cli = BarcodeValidatorCLI()
        yield cli

    finally:
        # Restore the original environment
        sys.argv = original_argv
        os.chdir(original_dir)

def test_galaxy_coi(cli_prepare_coi):
    """Test struct/taxon validation with BOLD data against BOLD web service."""
    dars = cli_prepare_coi.run()
    assert dars is not None, "Validation results should not be None"

if __name__ == '__main__':
    pytest.main()
