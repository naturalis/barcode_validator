import pytest
import io
from pathlib import Path
from unittest.mock import Mock
from nbitk.config import Config
from nbitk.Phylo.BOLDXLSXIO import Parser
from barcode_validator.barcode_validator import BarcodeValidator


@pytest.fixture
def mock_config():
    config = Mock(spec=Config)
    config.get.side_effect = lambda key: {
        'repo_owner': 'test_owner',
        'repo_name': 'test_repo',
        'repo_location': '/test/clone/path',
        'translation_table': 1,
        'ncbi_taxonomy': 'mock_ncbi.tar.gz',
        'bold_sheet_file': 'mock_bold.xlsx',
        'log_level': 'ERROR',
    }[key]
    return config

@pytest.fixture
def bold_xlsx_path():
    # Get the directory of the current test file
    current_dir = Path(__file__).parent
    # Construct the path to the bold.xlsx file
    return current_dir.parent / 'examples' / 'bold.xlsx'


@pytest.fixture
def parsed_tree(bold_xlsx_path):
    with open(bold_xlsx_path, 'rb') as file:
        parser = Parser(io.BytesIO(file.read()))
        return parser.parse()


def test_barcode_validator(parsed_tree, mock_config):
    validator = BarcodeValidator(mock_config)
    validator.bold_tree = parsed_tree

    # Test for the specific process ID
    node = validator.get_node_by_processid('BBIOP1028-24')

    assert node is not None, "Node not found for process ID 'BBIOP1028-24'"
    assert node.name == 'Bibio', f"Expected node name 'Bibio', but got '{node.name}'"
    assert node.taxonomic_rank == 'genus', f"Expected taxonomic rank 'genus', but got '{node.taxonomic_rank}'"


if __name__ == '__main__':
    pytest.main()