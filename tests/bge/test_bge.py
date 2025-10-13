import os
import pathlib
import pytest
from pathlib import Path
import sys

from barcode_validator.cli import BarcodeValidatorCLI

# do I still work?

@pytest.fixture
def data_dir():
    """Fixture to provide the test data directory."""
    # This assumes tests are run from the project root
    return Path("tests/bge/data")

@pytest.fixture
def input_fasta(data_dir):
    """Fixture to create a test FASTA file with sequences."""
    fasta_path = data_dir / "WK-3860_BSNHM190-con_seqs.fasta"
    return str(fasta_path)

@pytest.fixture
def input_csv(data_dir):
    csv_path = data_dir / "WK-3860_BSNHM190-metrics.csv"
    return str(csv_path)

@pytest.fixture
def bold_excel(data_dir):
    """Fixture to provide a test BOLD Excel file."""
    bold_path = data_dir / "bold_bge_container_plus-updated.xlsx"
    return str(bold_path)

@pytest.fixture
def cli_structural_coi(input_fasta, bold_excel, input_csv, data_dir):
    """
    In this first pass, we only do structural validation. We triage
    the output based on the structural validation results at the
    level of the sample (group_by_sample=true). We do this by aggregating
    the results for each sample (based on group_id_separator=_) and picking
    the best sequence. Hence, the output FASTA will contain one sequence
    per sample, the best one based on structural validation.
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
            "--input-file", input_fasta,
            "--csv-file", input_csv,
            "--mode", "structural",
            "--marker", "COI-5P",
            "--input-resolver", "format=bold",
            "--input-resolver", f"file={bold_excel}",
            "--output-fasta", f"{data_dir}/structval_out.fasta",
            "--output-tsv", f"{data_dir}/structval_out.tsv",
            "--triage-config", "group_id_separator=_",
            "--triage-config", "group_by_sample=true",
            "--log-level", "INFO"
        ]
        cli = BarcodeValidatorCLI()
        yield cli

    finally:
        # Restore the original environment
        sys.argv = original_argv
        os.chdir(original_dir)

@pytest.fixture
def cli_taxonomic_coi(bold_excel, input_csv, data_dir):
    """Fixture to run the CLI command."""

    # Save the original environment
    original_argv = sys.argv.copy()
    original_dir = os.getcwd()

    try:

        # Change to root of the repo
        script_dir = pathlib.Path(__file__).parent.parent.parent
        os.chdir(script_dir)

        # Get input fasta from data_dir fixture and output
        # from previous structural validation step. The CSV
        # is the same as before.
        input_fasta = str(data_dir / "structval_out.fasta")

        # Replace sys.argv with our test arguments
        sys.argv = [
            "barcode_validator",  # Program name
            "--input-file", input_fasta,
            "--csv-file", input_csv,
            "--mode", "taxonomic",
            "--marker", "COI-5P",
            "--input-resolver", "format=bold",
            "--input-resolver", f"file={bold_excel}",
            "--output-fasta", f"{data_dir}/taxonval_out.fasta",
            "--output-tsv", f"{data_dir}/taxonval_out.tsv",
            "--taxon-validation", "method=galaxy",
            "--taxon-validation", "rank=family",
            "--taxon-validation", "min_identity=0.8",
            "--taxon-validation", "max_target_seqs=1000",
            "--log-level", "INFO"
        ]
        cli = BarcodeValidatorCLI()
        yield cli

    finally:
        # Restore the original environment
        sys.argv = original_argv
        os.chdir(original_dir)


def test_bge_coi(cli_structural_coi, cli_taxonomic_coi):
    """Test struct/taxon validation with BOLD data against BOLD web service."""
    structval_results = cli_structural_coi.run()
    assert len(structval_results.results) == 4421, "Input FASTA had 4421 sequences, including empty ones"

    taxonval_results = cli_taxonomic_coi.run()
    assert len(taxonval_results.results) == 145, "Structval output FASTA had 145 sequences"


if __name__ == '__main__':
    pytest.main()
