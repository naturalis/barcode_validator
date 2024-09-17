import pytest
import io
from pathlib import Path
from nbitk.Phylo.BOLDXLSXIO import Parser
from barcode_validator.core import BarcodeValidator


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


def test_barcode_validator(parsed_tree):
    validator = BarcodeValidator()
    validator.bold_tree = parsed_tree

    # Test for the specific process ID
    node = validator.get_node_by_processid('BBIOP1028-24')

    assert node is not None, "Node not found for process ID 'BBIOP1028-24'"
    assert node.name == 'Bibio', f"Expected node name 'Bibio', but got '{node.name}'"
    assert node.taxonomic_rank == 'genus', f"Expected taxonomic rank 'genus', but got '{node.taxonomic_rank}'"


if __name__ == '__main__':
    pytest.main([__file__])