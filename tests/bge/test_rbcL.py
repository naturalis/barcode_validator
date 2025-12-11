import os
import pathlib
import pytest
from pathlib import Path
import sys

from barcode_validator.cli import BarcodeValidatorCLI

@pytest.fixture(autouse=True)
def check_galaxy_key():
    """Verify Galaxy API key is available"""
    if not os.environ.get('GALAXY_API_KEY'):
        pytest.skip("GALAXY_API_KEY not set in environment")

@pytest.fixture(autouse=True)
def check_galaxy_domain():
    """Verify Galaxy domain is available"""
    if not os.environ.get('GALAXY_DOMAIN'):
        pytest.skip("GALAXY_DOMAIN not set in environment")

@pytest.fixture
def data_dir():
    """Fixture to provide the test data directory."""
    # This assumes tests are run from the project root
    return Path("tests/bge/data")

@pytest.fixture
def input_fasta_rbcl(data_dir):
    """Fixture to create a test FASTA file with rbcL sequences."""
    fasta_path = data_dir / "rbcl_barcodes_YG-4378.fasta"
    return str(fasta_path)

@pytest.fixture
def bold_excel(data_dir):
    """Fixture to provide a test BOLD Excel file."""
    bold_path = data_dir / "bold_bge_container_plus-updated.xlsx"
    return str(bold_path)

@pytest.fixture
def cli_structural_rbcl(input_fasta_rbcl, bold_excel, data_dir):
    """
    Structural validation for rbcL sequences.
    """

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
            "--input-file", input_fasta_rbcl,
            "--mode", "structural",
            "--marker", "rbcL",
            "--input-resolver", "format=bold",
            "--input-resolver", f"file={bold_excel}",
            "--output-fasta", f"{data_dir}/rbcl_structval_out.fasta",
            "--output-tsv", f"{data_dir}/rbcl_structval_out.tsv",
            "--log-level", "INFO"
        ]
        cli = BarcodeValidatorCLI()
        yield cli

    finally:
        # Restore the original environment
        sys.argv = original_argv
        os.chdir(original_dir)

@pytest.fixture
def cli_taxonomic_rbcl(bold_excel, data_dir):
    """Fixture to run the CLI command for taxonomic validation of rbcL sequences."""

    # Save the original environment
    original_argv = sys.argv.copy()
    original_dir = os.getcwd()

    try:

        # Change to root of the repo
        script_dir = pathlib.Path(__file__).parent.parent.parent
        os.chdir(script_dir)

        # Get input fasta from data_dir fixture and output
        # from previous structural validation step.
        input_fasta = str(data_dir / "rbcl_structval_out.fasta")

        # Replace sys.argv with our test arguments
        sys.argv = [
            "barcode_validator",  # Program name
            "--input-file", input_fasta,
            "--mode", "taxonomic",
            "--marker", "rbcL",
            "--input-resolver", "format=bold",
            "--input-resolver", f"file={bold_excel}",
            "--output-fasta", f"{data_dir}/rbcl_taxonval_out.fasta",
            "--output-tsv", f"{data_dir}/rbcl_taxonval_out.tsv",
            "--taxon-validation", "method=galaxy",
            "--taxon-validation", "rank=family",
            "--taxon-validation", "min_identity=0.8",
            "--taxon-validation", "max_target_seqs=1000",
            "--galaxy-blast", "database=Genbank rbcL (2023-11-15)",
            "--log-level", "INFO"
        ]
        cli = BarcodeValidatorCLI()
        yield cli

    finally:
        # Restore the original environment
        sys.argv = original_argv
        os.chdir(original_dir)


def test_rbcl(cli_structural_rbcl, cli_taxonomic_rbcl):
    """Test struct/taxon validation with rbcL data against Galaxy BLAST service."""
    structval_results = cli_structural_rbcl.run()
    assert len(structval_results.results) == 72, "Input FASTA had 72 rbcL sequences"

    taxonval_results = cli_taxonomic_rbcl.run()
    # The number of sequences after taxonomic validation will depend on the BLAST results
    # We just check that the test completes without errors
    assert taxonval_results is not None, "Taxonomic validation completed"


if __name__ == '__main__':
    pytest.main()
