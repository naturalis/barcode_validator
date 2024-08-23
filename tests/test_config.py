import pytest
import os
import yaml
from pathlib import Path
from unittest.mock import patch, mock_open
from barcode_validator.config import Config


@pytest.fixture
def valid_config_data():
    return {
        "repo_owner": "naturalis",
        "repo_name": "barcode_validator",
        "hmm_file": "/path/to/hmm_file.hmm",
        "bold_sheet_file": "/path/to/bold_sheet.xlsx",
        "level": "family",
        "constrain": "class",
        "num_threads": 24,
        "evalue": 1e-5,
        "max_target_seqs": 10,
        "word_size": 28,
        "BLASTDB_LMDB_MAP_SIZE": 180000000000,
        "translation_table": 5,
        "log_level": "INFO"
    }


@pytest.fixture
def config_instance():
    return Config()


def test_file_existence(config_instance, valid_config_data, tmp_path):
    # Create temporary files
    hmm_file = tmp_path / "test.hmm"
    hmm_file.touch()
    bold_sheet = tmp_path / "test.xlsx"
    bold_sheet.touch()

    valid_config_data["hmm_file"] = str(hmm_file)
    valid_config_data["bold_sheet_file"] = str(bold_sheet)

    with patch("builtins.open", mock_open(read_data=yaml.dump(valid_config_data))):
        with patch("os.path.exists", lambda x: x in [str(hmm_file), str(bold_sheet), "dummy_path"]):
            config_instance.load_config("dummy_path")

    assert os.path.exists(config_instance.get("hmm_file"))
    assert os.path.exists(config_instance.get("bold_sheet_file"))


def test_github_repo_url(config_instance, valid_config_data):
    with patch("builtins.open", mock_open(read_data=yaml.dump(valid_config_data))):
        with patch("os.path.exists", return_value=True):
            config_instance.load_config("dummy_path")

    repo_url = f"https://github.com/{config_instance.get('repo_owner')}/{config_instance.get('repo_name')}"
    response = pytest.importorskip("requests").get(repo_url)
    assert response.status_code == 200


@pytest.mark.parametrize("key", ["level", "constrain"])
def test_valid_taxonomic_levels(config_instance, valid_config_data, key):
    with patch("builtins.open", mock_open(read_data=yaml.dump(valid_config_data))):
        with patch("os.path.exists", return_value=True):
            config_instance.load_config("dummy_path")

    valid_levels = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    assert config_instance.get(key).lower() in valid_levels


@pytest.mark.parametrize("key", ["num_threads", "evalue", "max_target_seqs", "word_size", "BLASTDB_LMDB_MAP_SIZE", "translation_table"])
def test_numeric_values(config_instance, valid_config_data, key):
    with patch("builtins.open", mock_open(read_data=yaml.dump(valid_config_data))):
        with patch("os.path.exists", return_value=True):
            config_instance.load_config("dummy_path")

    assert isinstance(config_instance.get(key), (int, float))


def test_valid_log_level(config_instance, valid_config_data):
    with patch("builtins.open", mock_open(read_data=yaml.dump(valid_config_data))):
        with patch("os.path.exists", return_value=True):
            config_instance.load_config("dummy_path")

    valid_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
    assert config_instance.get("log_level").upper() in valid_levels


def test_load_config_file_not_found(config_instance):
    with pytest.raises(FileNotFoundError):
        config_instance.load_config("non_existent_file.yml")


def test_load_config_invalid_yaml(config_instance):
    invalid_yaml = "invalid: yaml: content:"
    with patch("builtins.open", mock_open(read_data=invalid_yaml)):
        with patch("os.path.exists", return_value=True):
            with pytest.raises(yaml.YAMLError):
                config_instance.load_config("dummy_path")


def test_get_nonexistent_key(config_instance, valid_config_data):
    with patch("builtins.open", mock_open(read_data=yaml.dump(valid_config_data))):
        with patch("os.path.exists", return_value=True):
            config_instance.load_config("dummy_path")

    assert config_instance.get("non_existent_key") is None
    assert config_instance.get("non_existent_key", "default") == "default"


def test_contains_method(config_instance, valid_config_data):
    with patch("builtins.open", mock_open(read_data=yaml.dump(valid_config_data))):
        with patch("os.path.exists", return_value=True):
            config_instance.load_config("dummy_path")

    assert "repo_owner" in config_instance
    assert "non_existent_key" not in config_instance


def test_config_singleton():

    # Reset the singleton instance
    Config.reset()

    # Construct and check path to expected default config file
    current_dir = Path(__file__).parent
    default_config_path = current_dir.parent / 'config' / 'config.yml'
    if not default_config_path.exists():
        raise FileNotFoundError(f"Default config file not found at {default_config_path}")

    # Load the config file twice
    config1 = Config()
    config1.load_config(str(default_config_path))
    config2 = Config()

    # Additional checks
    assert config1 is config2, "Config objects are not the same instance"
    assert config1.config_path == str(default_config_path), "Config path is incorrect"
    assert config1.config_data is not None, "Config data is None"
    assert 'constrain' in config1.config_data, "'constrain' key is missing from config data"
