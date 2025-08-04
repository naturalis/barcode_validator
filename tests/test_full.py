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
    return Path("tests/data")

@pytest.fixture
def input_fasta_bold(data_dir):
    """Fixture to create a test FASTA file with sequences."""
    fasta_path = data_dir / "test_full.fasta"
    return str(fasta_path)

@pytest.fixture
def bold_excel(data_dir):
    """Fixture to provide a test BOLD Excel file."""
    bold_path = data_dir / "bold.xlsx"
    return str(bold_path)

@pytest.fixture
def input_tsv_nsr(data_dir):
    """Fixture to create a test TSV file with sequences."""
    tsv_path = data_dir / "structval_nsr.tsv"
    return str(tsv_path)

@pytest.fixture
def nsr_dwca(data_dir):
    """Fixture to provide a test NSR archive."""
    nsr_path = data_dir / "nsr-20250207.dwca.zip"
    return str(nsr_path)

@pytest.fixture
def test_cli_bold_bold(input_fasta_bold, bold_excel):
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
            "--mode", "both",
            "--marker", "COI-5P",
            "--exp-taxonomy-type", "bold",
            "--exp-taxonomy", bold_excel,
            "--reflib-taxonomy-type", "bold",
            "--reflib-taxonomy", bold_excel,
            "--output-format", "tsv",
            "--taxon-validation", "method=bold",
            "--taxon-validation", "rank=family",
            "--taxon-validation", "extent=class",
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


@pytest.fixture
def test_cli_nsr(input_tsv_nsr, nsr_dwca, bold_excel):
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
            "--mode", "both",
            "--marker", "COI-5P",
            "--exp-taxonomy-type", "nsr",
            "--exp-taxonomy", nsr_dwca,
            "--reflib-taxonomy-type", "bold",
            "--reflib-taxonomy", bold_excel,
            "--output-format", "tsv",
            "--taxon-validation", "method=bold",
            "--taxon-validation", "rank=family",
            "--taxon-validation", "extent=class",
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

def test_validation_bold_bold(test_cli_bold_bold):
    """Test struct/taxon validation with BOLD data against BOLD web service."""
    dars = test_cli_bold_bold.run()
    assert dars is not None, "Validation results should not be None"
    assert dars.results[0].full_length == 567, "Total length"
    assert dars.results[0].seq_length == 522, "Length within marker region"
    assert len(dars.results[0].stop_codons) == 0, "No early stop codons"
    assert dars.results[0].full_ambiguities == 0, "No ambiguities in sequence"
    assert dars.results[0].ambiguities == 0, "Taxon should be resolved"
    assert dars.results[0].obs_taxon[0].name == dars.results[0].exp_taxon.name, "Observed taxon matches expected taxon"
    assert dars.results[0].species.name == "Halesus tessellatus", "Asserted species name"

def test_validation_nsr_bold(test_cli_nsr):
    dars = test_cli_nsr.run()
    assert dars is not None, "Validation results should not be None"

if __name__ == '__main__':
    pytest.main()