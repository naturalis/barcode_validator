import os
import pathlib
import pytest
from pathlib import Path
import sys

from barcode_validator.cli import BarcodeValidatorCLI


@pytest.fixture
def data_dir():
    """Fixture to provide the test data directory."""
    # This assumes tests are run from the project root
    return Path("tests/bge/data")

@pytest.fixture
def input_fasta(data_dir):
    """Fixture to create a test FASTA file with sequences."""
    fasta_path = data_dir / "mge_fastp_r13s100_nocontam.fasta"
    return str(fasta_path)

@pytest.fixture
def input_csv(data_dir):
    csv_path = data_dir / "mge_fastp_r13s100_nocontam.csv"
    return str(csv_path)

@pytest.fixture
def input_yaml(data_dir):
    yaml_path = data_dir / "mge_fastp_r13s100_nocontam.yaml"
    return str(yaml_path)

@pytest.fixture
def bold_excel(data_dir):
    """Fixture to provide a test BOLD Excel file."""
    bold_path = data_dir / "bold.xlsx"
    return str(bold_path)

@pytest.fixture
def cli_prepare_coi(input_fasta, bold_excel, input_csv, input_yaml):
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
            "--csv-file", input_csv,
            "--yaml-file", input_yaml,
            "--mode", "both",
            "--marker", "COI-5P",
            "--input-resolver", "format=bold",
            "--input-resolver", f"file={bold_excel}",
            "--output-fasta", "bge_out.fasta",
            "--output-tsv", "bge_out.tsv",
            "--taxon-validation", "method=bold",
            "--taxon-validation", "rank=family",
            "--taxon-validation", "min_identity=0.8",
            "--taxon-validation", "max_target_seqs=100",
            "--triage-config", "group_id_separator=_",
            "--triage-config", "group_by_sample=true",
            "--log-level", "ERROR"
        ]
        cli = BarcodeValidatorCLI()
        yield cli

    finally:
        # Restore the original environment
        sys.argv = original_argv
        os.chdir(original_dir)

def test_bge_coi(cli_prepare_coi):
    """Test struct/taxon validation with BOLD data against BOLD web service."""
    dars = cli_prepare_coi.run()
    assert dars is not None, "Validation results should not be None"

if __name__ == '__main__':
    pytest.main()