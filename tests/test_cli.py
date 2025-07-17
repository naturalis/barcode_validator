import pytest
import os
from pathlib import Path
from unittest.mock import patch, MagicMock
from barcode_validator.cli import BarcodeValidatorCLI
from nbitk.config import Config
from barcode_validator.orchestrator import ValidationOrchestrator

# Test data path handling
TEST_DATA_DIR = Path(__file__).parent / "data"
CSC_SAMPLE = TEST_DATA_DIR / "csc_sample.tsv"
BOLD_SAMPLE = TEST_DATA_DIR / "BGE00146_MGE-BGE_r1_1.3_1.5_s50_100.fasta"
NSR_ARCHIVE = TEST_DATA_DIR / "nsr-20250207.dwca.zip"
BOLD_SHEET = TEST_DATA_DIR / "bold.xlsx"
OUTPUT_TSV = TEST_DATA_DIR / "output.tsv"
OUTPUT_FASTA = TEST_DATA_DIR / "valid.fasta"


@pytest.fixture
def base_args():
    """Provides base arguments for testing"""
    return ["--input_file", str(BOLD_SAMPLE)]


@pytest.fixture
def mock_orchestrator():
    """Provides mocked ValidationOrchestrator"""
    with patch('barcode_validator.cli.ValidationOrchestrator') as mock:
        instance = MagicMock()
        mock.return_value = instance
        yield instance


def test_basic_argument_parsing():
    """Test basic argument parsing with minimal required arguments"""
    with patch('sys.argv', ['script.py'] + ["--input-file", str(BOLD_SAMPLE)]):
        cli = BarcodeValidatorCLI()
        assert cli.args.input_file == BOLD_SAMPLE
        assert cli.args.mode == "both"  # Default value
        assert cli.args.exp_taxonomy_type == "nsr"  # Default value
        assert cli.args.output_format == "tsv"  # Default value
        assert cli.args.triage is False  # Default value


def test_all_arguments():
    """Test parsing of all possible arguments"""
    with patch('sys.argv', ['script.py'] + [
        "--input-file", str(BOLD_SAMPLE),
        "--csv-file", "analytics.csv",
        "--yaml-file", "metadata.yml",
        "--reflib-taxonomy", "ncbi.tar.gz",
        "--exp-taxonomy", str(BOLD_SHEET),
        "--exp-taxonomy-type", "bold",
        "--mode", "taxonomic",
        "--output-format", "fasta",
        "--triage",
        "--log-level", "DEBUG"
    ]):
        cli = BarcodeValidatorCLI()
        args = cli.args
        assert args.input_file == BOLD_SAMPLE
        assert args.csv_file == Path("analytics.csv")
        assert args.yaml_file == Path("metadata.yml")
        assert args.reflib_taxonomy == Path("ncbi.tar.gz")
        assert args.exp_taxonomy == BOLD_SHEET
        assert args.exp_taxonomy_type == "bold"
        assert args.mode == "taxonomic"
        assert args.output_format == "fasta"
        assert args.triage is True
        assert args.log_level == "DEBUG"


def test_validation_mode_options():
    """Test different validation mode options"""
    modes = ["structural", "taxonomic", "both"]
    for mode in modes:
        with patch('sys.argv', ['script.py'] + [
            "--input-file", str(BOLD_SAMPLE),
            "--mode", mode
        ]):
            cli = BarcodeValidatorCLI()
            assert cli.args.mode == mode


def test_taxonomy_type_validation():
    """Test validation of taxonomy type options"""
    with pytest.raises(SystemExit):
        with patch('sys.argv', ['script.py'] + [
            "--input-file", str(BOLD_SAMPLE),
            "--exp-taxonomy-type", "invalid"  # Invalid type
        ]):
            BarcodeValidatorCLI()


def test_required_arguments():
    """Test handling of missing required arguments"""
    required_args = ["--input-file"]

    for missing_arg in required_args:
        args = ["--input-file", str(BOLD_SAMPLE)]
        # Remove the argument and its value
        idx = args.index(missing_arg)
        test_args = args[:idx] + args[idx + 2:]

        with pytest.raises(SystemExit):
            with patch('sys.argv', ['script.py'] + test_args):
                BarcodeValidatorCLI()


def test_file_existence_checking():
    """Test handling of non-existent input files"""
    with pytest.raises(SystemExit):
        with patch('sys.argv', ['script.py'] + ["--input-file", str(BOLD_SAMPLE)]):
            cli = BarcodeValidatorCLI()
            cli.run()


def test_integration_with_orchestrator(mock_orchestrator):
    """Test integration with ValidationOrchestrator"""
    with patch('sys.argv', ['script.py'] + [
        "--input-file", str(BOLD_SAMPLE)
    ]):
        cli = BarcodeValidatorCLI()
        cli.run()

        # Verify orchestrator calls
        mock_orchestrator.validate_file.assert_called_once()
        mock_orchestrator.write_results.assert_called_once()

        # Verify call arguments
        validate_args = mock_orchestrator.validate_file.call_args[0]
        assert str(validate_args[0]) == str(BOLD_SAMPLE)
        assert validate_args[1] is None  # csv_file
        assert validate_args[2] is None  # yaml_file


def test_log_level_setting():
    """Test setting of different log levels"""
    log_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]

    for level in log_levels:
        with patch('sys.argv', ['script.py'] + [
            "--input-file", str(BOLD_SAMPLE),
            "--log-level", level
        ]):
            cli = BarcodeValidatorCLI()
            assert cli.config.get("log_level") == level

def test_error_handling(mock_orchestrator):
    """Test handling of various error conditions"""
    # Test orchestrator error
    mock_orchestrator.validate_file.side_effect = Exception("Validation failed")

    with pytest.raises(SystemExit) as exc_info:
        with patch('sys.argv', ['script.py'] + [
            "--input-file", str(BOLD_SAMPLE)
        ]):
            cli = BarcodeValidatorCLI()
            cli.run()

    assert exc_info.value.code == 1


if __name__ == '__main__':
    pytest.main(['-v'])