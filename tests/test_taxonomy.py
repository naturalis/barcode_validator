import pytest
import requests
import tempfile
from pathlib import Path
from Bio import SeqRecord
from barcode_validator.config import Config
from barcode_validator.taxonomy import read_bold_taxonomy, read_ncbi_taxonomy, get_tip_by_processid, run_localblast


@pytest.fixture
def bold_tree():
    bold_file = Path(__file__).parent.parent / "examples" / "bold.xlsx"
    return read_bold_taxonomy(bold_file)


@pytest.fixture
def ncbi_tree():
    ncbi_url = "http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

    with tempfile.NamedTemporaryFile(suffix='.tar.gz', delete=False) as temp_file:
        response = requests.get(ncbi_url)
        response.raise_for_status()  # Raise an exception for bad status codes
        temp_file.write(response.content)
        temp_file_path = temp_file.name

    try:
        return read_ncbi_taxonomy(temp_file_path)
    finally:
        Path(temp_file_path).unlink()  # Delete the temporary file


def test_get_tip_by_processid(bold_tree):
    for tip in bold_tree.get_terminals():
        if tip.name == "Halesus tessellatus":
            assert tip.guids["BSNHM166-24"] == "BGE_00002_F11"
            assert tip.taxonomic_rank == "species"
            break
    tip = get_tip_by_processid("BSNHM166-24", bold_tree)
    assert tip.name == "Halesus tessellatus"


def test_none_sequence_blast():
    seq = SeqRecord.SeqRecord('')
    con = Config()
    con.initialized = True
    con.config_data = {"constrain": "family"}
    assert run_localblast(seq, None, None, con) is None


def test_bold_tree_ancestry(bold_tree):
    process_id = "BGENL1996-24"
    tip = get_tip_by_processid(process_id, bold_tree)

    class_ancestor = next(node for node in bold_tree.get_path(tip) if node.taxonomic_rank == "class")
    assert class_ancestor.name == "Insecta"


def test_regression_1(bold_tree):
    process_id = "BGENL191-23"
    tip = get_tip_by_processid(process_id, bold_tree)
    class_ancestor = next(node for node in bold_tree.root.get_path(tip) if node.taxonomic_rank == "class")
    assert class_ancestor.name == "Insecta"


def test_ncbi_tree_structure(ncbi_tree):
    insecta_node = next(node for node in ncbi_tree.get_nonterminals()
                        if node.name == "Insecta" and node.taxonomic_rank == "class")

    assert insecta_node is not None
    assert insecta_node.guids['taxon'] == "50557"

    terminal_descendants = list(insecta_node.get_terminals())
    assert len(terminal_descendants) > 0
