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
CONFIG_YML = TEST_DATA_DIR / "config.yml"
OUTPUT_TSV = TEST_DATA_DIR / "output.tsv"
OUTPUT_FASTA = TEST_DATA_DIR / "valid.fasta"


@pytest.fixture
def base_args():
    """Provides base arguments for testing"""
    return [
        "--config", str(CONFIG_YML),
        "--input_file", str(BOLD_SAMPLE),
        "--output_file", str(OUTPUT_TSV)
    ]


@pytest.fixture
def mock_orchestrator():
    """Provides mocked ValidationOrchestrator"""
    with patch('barcode_validator.cli.ValidationOrchestrator') as mock:
        instance = MagicMock()
        mock.return_value = instance
        yield instance


def test_basic_argument_parsing():
    """Test basic argument parsing with minimal required arguments"""
    with patch('sys.argv', ['script.py'] + [
        "--config", str(CONFIG_YML),
        "--input_file", str(BOLD_SAMPLE),
        "--output_file", str(OUTPUT_TSV)
    ]):
        cli = BarcodeValidatorCLI()
        assert cli.args.config == str(CONFIG_YML)
        assert cli.args.input_file == str(BOLD_SAMPLE)
        assert cli.args.output_file == str(OUTPUT_TSV)
        assert cli.args.mode == "both"  # Default value
        assert not cli.args.emit_valid_fasta  # Default False


def test_all_arguments():
    """Test parsing of all possible arguments"""
    with patch('sys.argv', ['script.py'] + [
        "--config", str(CONFIG_YML),
        "--input_file", str(BOLD_SAMPLE),
        "--csv_file", "analytics.csv",
        "--yaml_file", "metadata.yml",
        "--reflib_taxonomy", "ncbi.tar.gz",
        "--exp_taxonomy", str(BOLD_SHEET),
        "--exp_taxonomy_type", "bold",
        "--mode", "taxonomic",
        "--emit_valid_fasta",
        "--output_fasta", str(OUTPUT_FASTA),
        "--output_file", str(OUTPUT_TSV),
        "--log_level", "DEBUG"
    ]):
        cli = BarcodeValidatorCLI()
        args = cli.args
        assert args.csv_file == "analytics.csv"
        assert args.yaml_file == "metadata.yml"
        assert args.reflib_taxonomy == "ncbi.tar.gz"
        assert args.exp_taxonomy == str(BOLD_SHEET)
        assert args.exp_taxonomy_type == "bold"
        assert args.mode == "taxonomic"
        assert args.emit_valid_fasta
        assert args.output_fasta == str(OUTPUT_FASTA)
        assert args.log_level == "DEBUG"


def test_config_loading_and_override(tmp_path):
    """Test config loading and command line overrides"""
    # Create a temporary config file
    config_content = """
log_level: INFO
validate_taxonomy: true
validate_structure: true
"""
    config_file = tmp_path / "test_config.yml"
    config_file.write_text(config_content)

    with patch('sys.argv', ['script.py'] + [
        "--config", str(config_file),
        "--input_file", str(BOLD_SAMPLE),
        "--output_file", str(OUTPUT_TSV),
        "--log_level", "DEBUG",  # Override config
        "--reflib_taxonomy", "new_ncbi.tar.gz"  # Add new value
    ]):
        cli = BarcodeValidatorCLI()
        assert cli.config.get("log_level") == "DEBUG"  # Overridden
        assert cli.config.get("ncbi_taxonomy") == "new_ncbi.tar.gz"  # Added


def test_validation_mode_options():
    """Test different validation mode options"""
    modes = ["structural", "taxonomic", "both"]
    for mode in modes:
        with patch('sys.argv', ['script.py'] + [
            "--config", str(CONFIG_YML),
            "--input_file", str(BOLD_SAMPLE),
            "--output_file", str(OUTPUT_TSV),
            "--mode", mode
        ]):
            cli = BarcodeValidatorCLI()
            assert cli.args.mode == mode


def test_taxonomy_type_validation():
    """Test validation of taxonomy type options"""
    with pytest.raises(SystemExit):
        with patch('sys.argv', ['script.py'] + [
            "--config", str(CONFIG_YML),
            "--input_file", str(BOLD_SAMPLE),
            "--output_file", str(OUTPUT_TSV),
            "--exp_taxonomy_type", "invalid"  # Invalid type
        ]):
            BarcodeValidatorCLI()


def test_required_arguments():
    """Test handling of missing required arguments"""
    required_args = ["--config", "--input_file", "--output_file"]

    for missing_arg in required_args:
        args = [
            "--config", str(CONFIG_YML),
            "--input_file", str(BOLD_SAMPLE),
            "--output_file", str(OUTPUT_TSV)
        ]
        # Remove the argument and its value
        idx = args.index(missing_arg)
        test_args = args[:idx] + args[idx + 2:]

        with pytest.raises(SystemExit):
            with patch('sys.argv', ['script.py'] + test_args):
                BarcodeValidatorCLI()


def test_file_existence_checking():
    """Test handling of non-existent input files"""
    with pytest.raises(SystemExit):
        with patch('sys.argv', ['script.py'] + [
            "--config", "nonexistent.yml",
            "--input_file", str(BOLD_SAMPLE),
            "--output_file", str(OUTPUT_TSV)
        ]):
            cli = BarcodeValidatorCLI()
            cli.run()


def test_integration_with_orchestrator(mock_orchestrator):
    """Test integration with ValidationOrchestrator"""
    with patch('sys.argv', ['script.py'] + [
        "--config", str(CONFIG_YML),
        "--input_file", str(BOLD_SAMPLE),
        "--output_file", str(OUTPUT_TSV),
        "--emit_valid_fasta",
        "--output_fasta", str(OUTPUT_FASTA)
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
            "--config", str(CONFIG_YML),
            "--input_file", str(BOLD_SAMPLE),
            "--output_file", str(OUTPUT_TSV),
            "--log_level", level
        ]):
            cli = BarcodeValidatorCLI()
            assert cli.config.get("log_level") == level


def test_output_handling(mock_orchestrator, tmp_path):
    """Test handling of output files and directories"""
    output_tsv = tmp_path / "results.tsv"
    output_fasta = tmp_path / "valid.fasta"

    with patch('sys.argv', ['script.py'] + [
        "--config", str(CONFIG_YML),
        "--input_file", str(BOLD_SAMPLE),
        "--output_file", str(output_tsv),
        "--emit_valid_fasta",
        "--output_fasta", str(output_fasta)
    ]):
        cli = BarcodeValidatorCLI()
        cli.run()

        # Verify output parameters in orchestrator call
        write_args = mock_orchestrator.write_results.call_args[0]
        assert str(write_args[1]) == str(output_tsv)
        assert write_args[2] is True  # emit_valid_fasta
        assert str(write_args[3]) == str(output_fasta)


def test_error_handling(mock_orchestrator):
    """Test handling of various error conditions"""
    # Test orchestrator error
    mock_orchestrator.validate_file.side_effect = Exception("Validation failed")

    with pytest.raises(SystemExit) as exc_info:
        with patch('sys.argv', ['script.py'] + [
            "--config", str(CONFIG_YML),
            "--input_file", str(BOLD_SAMPLE),
            "--output_file", str(OUTPUT_TSV)
        ]):
            cli = BarcodeValidatorCLI()
            cli.run()

    assert exc_info.value.code == 1


if __name__ == '__main__':
    pytest.main(['-v'])