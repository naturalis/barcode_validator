import pytest
from nbitk.Taxon import Taxon

from barcode_validator.result import DNAAnalysisResult


@pytest.fixture
def dna_result():
    return DNAAnalysisResult("test_process_id")


@pytest.mark.parametrize(
    "seq_length,ambiguities,full_length,full_ambiguities,expected_barcode_rank,expected_full_rank,expected_messages", [
        (650, 0, 1501, 0, 1, 1,
         ["\U0001F947\U0001F31F BIN compliant, perfect", "\U0001F947 Excellent full length, no ambiguities"]),
        (500, 0, 1001, 0, 2, 2, ["\U0001F947 BIN compliant", "\U0001F948 Good full length, no ambiguities"]),
        (650, 3, 1501, 3, 3, 3, ["\U0001F948 BIN compliant", "\u2753 Marker may be chimeric",
                                 "\U0001F947 Excellent full length", "\u2757 Some ambiguities in full sequence"]),
        (500, 3, 1001, 3, 4, 4,
         ["\U0001F949 BIN compliant", "\u2753 Marker may be chimeric", "\U0001F948 Good full length",
          "\u2757 Some ambiguities in full sequence"]),
        (450, 0, 1000, 3, 5, 4, ["\u26A0 Useful marker sequence", "\u2757 Not BIN compliant",
                                 "\U0001F948 Good full length", "\u2757 Some ambiguities in full sequence"]),
        (350, 0, 900, 0, 6, 5, ["\u26A0 Useful marker sequence", "\u203C Not BIN compliant",
                                "\U0001F949 Acceptable full length"]),
        (250, 0, 800, 5, 7, 5, ["\u26A0 Useful marker sequence", "\u203C Not BIN compliant",
                                "\U0001F949 Acceptable full length", "\u2757 Some ambiguities in full sequence"]),
        (200, 10, 700, 20, 8, 6, ["\u26D4 Unacceptable marker sequence", "\u26D4 Unacceptable full sequence"]),
    ])
def test_ranks_and_messages(dna_result, seq_length, ambiguities, full_length, full_ambiguities,
                            expected_barcode_rank, expected_full_rank, expected_messages):
    dna_result.seq_length = seq_length
    dna_result.ambiguities = ambiguities
    dna_result.full_length = full_length
    dna_result.full_ambiguities = full_ambiguities

    barcode_rank, full_rank, messages = dna_result.calculate_ranks(verbosity=3)

    assert barcode_rank == expected_barcode_rank
    assert full_rank == expected_full_rank
    for message in expected_messages:
        assert message in messages


@pytest.mark.parametrize("verbosity,expected_message_count", [
    (1, 1),  # Only errors
    (2, 2),  # Errors and warnings
    (3, 4),  # Errors, warnings, and info
])
def test_verbosity_levels(dna_result, verbosity, expected_message_count):
    dna_result.seq_length = 650
    dna_result.ambiguities = 3
    dna_result.full_length = 1501
    dna_result.full_ambiguities = 3
    _, _, messages = dna_result.calculate_ranks(verbosity=verbosity)
    assert len(messages.split('\n')) == expected_message_count


def test_edge_case_marker(dna_result):
    # Test when seq_length or ambiguities are None
    dna_result.full_length = 1500
    dna_result.full_ambiguities = 0
    barcode_rank, full_rank, messages = dna_result.calculate_ranks(verbosity=1)
    assert barcode_rank == 8
    assert full_rank == 1
    assert "\u26D4 Unacceptable marker sequence" in messages


def test_columns(dna_result):
    columns = dna_result.result_fields()
    assert len(columns) == 13


def test_values(dna_result):
    dna_result.seq_length = 650
    dna_result.full_length = 1501
    dna_result.obs_taxon = [ Taxon(name="Family", taxonomic_rank='family') ]
    dna_result.exp_taxon = Taxon(name="Family", taxonomic_rank='family')
    dna_result.species = Taxon(name="Genus species", taxonomic_rank='species')
    dna_result.stop_codons = [ 678 ]
    dna_result.ambiguities = 3
    dna_result.full_ambiguities = 3
    dna_result.level = 'family'
    dna_result.ancillary = { 'foo': 'bar' }
    dna_result.error = 'error message'
    values = dna_result.get_values()
    assert len(values) == 14


def test_edge_case_full(dna_result):
    # Test when full_length or full_ambiguities are None
    dna_result.seq_length = 650
    dna_result.ambiguities = 0
    barcode_rank, full_rank, messages = dna_result.calculate_ranks(verbosity=1)
    assert barcode_rank == 1
    assert full_rank == 6
    assert "\u26D4 Unacceptable full sequence" in messages


if __name__ == '__main__':
    pytest.main()