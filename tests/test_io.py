import pytest
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from barcode_validator.orchestrator import ValidationOrchestrator
from barcode_validator.config.schema_config import SchemaConfig

# Test data path handling
TEST_DATA_DIR = Path(__file__).parent / "data"
CSC_SAMPLE = TEST_DATA_DIR / "csc_sample.tsv"
BOLD_SAMPLE = TEST_DATA_DIR / "BGE00146_MGE-BGE_r1_1.3_1.5_s50_100.fasta"


@pytest.fixture
def config():
    """Test fixture that provides a minimal Config object"""
    conf = SchemaConfig()
    conf.set('log_level', 'DEBUG')
    return conf


@pytest.fixture
def orchestrator(config):
    """Fixture providing orchestrator instance"""
    return ValidationOrchestrator(config)


def test_csc_file_parsing(orchestrator):
    """Test parsing of CSC format records"""
    records = list(orchestrator._parse_input(CSC_SAMPLE))

    # Test basic record count
    assert len(records) == 10, "Should parse 10 records from sample data"

    # Test record types
    assert all(isinstance(r, SeqRecord) for r in records), "All records should be SeqRecord objects"


def test_complete_record_parsing(orchestrator):
    """Test parsing of a complete record with all fields"""
    records = list(orchestrator._parse_input(CSC_SAMPLE))
    complete_record = next(r for r in records if r.id == "25281")

    # Test sequence content
    assert len(complete_record.seq) > 0, "Sequence should not be empty"
    assert str(complete_record.seq).startswith("AACAATATATCTAATCTTCGG"), "Sequence content mismatch"

    # Test annotations
    bcdm_fields = complete_record.annotations.get('bcdm_fields', {})
    assert bcdm_fields['marker_code'] == "COI-5P"
    assert bcdm_fields['identification'] == "Geophilus carpophagus"
    assert bcdm_fields['kingdom'] == "Animalia"
    assert bcdm_fields['identification_rank'].lower() == "species"


def test_null_value_handling(orchestrator):
    """Test handling of records with null values"""
    records = list(orchestrator._parse_input(CSC_SAMPLE))
    null_record = next(r for r in records if r.id == "507")

    bcdm_fields = null_record.annotations.get('bcdm_fields', {})
    assert bcdm_fields['kingdom'] == "null"
    assert bcdm_fields['identification_rank'] == "null"


def test_sequence_verification(orchestrator):
    """Test sequence content verification"""
    records = list(orchestrator._parse_input(CSC_SAMPLE))

    # Verify that all sequences are DNA
    for record in records:
        seq_str = str(record.seq).upper()
        valid_chars = set('ACGTN-')  # Basic DNA characters plus N and gaps
        assert all(c in valid_chars for c in seq_str), f"Invalid characters in sequence {record.id}"


def test_consistent_field_presence(orchestrator):
    """Test that all records have consistent field presence"""
    records = list(orchestrator._parse_input(CSC_SAMPLE))

    required_fields = {'marker_code', 'identification', 'kingdom', 'identification_rank', 'taxonomy_notes'}

    for record in records:
        bcdm_fields = record.annotations.get('bcdm_fields', {})
        assert all(field in bcdm_fields for field in required_fields), \
            f"Missing required fields in record {record.id}"


def test_marker_code_consistency(orchestrator):
    """Test consistency of marker codes"""
    records = list(orchestrator._parse_input(CSC_SAMPLE))

    # All records in sample should be COI-5P
    for record in records:
        assert record.annotations['bcdm_fields']['marker_code'] == "COI-5P", \
            f"Incorrect marker code in record {record.id}"


def test_invalid_file_handling(orchestrator):
    """Test handling of invalid file paths"""
    with pytest.raises(ValueError):
        list(orchestrator._parse_input(TEST_DATA_DIR / "nonexistent.xyz"))


def test_bold_fasta_parsing(orchestrator):
    """Test parsing of BOLD format FASTA records"""
    records = list(orchestrator._parse_input(BOLD_SAMPLE))

    # Should have two samples with six variants each
    assert len(records) == 10, "Should parse 10 records from sample data"
    assert all(isinstance(r, SeqRecord) for r in records), "All records should be SeqRecord objects"

    # Verify record identifiers follow BOLD pattern
    assert any(r.id.startswith("BHNHM001-24") for r in records), "Should find BHNHM001-24 records"
    assert any(r.id.startswith("BHNHM002-24") for r in records), "Should find BHNHM002-24 records"


def test_bold_sequence_content(orchestrator):
    """Test sequence content from BOLD FASTA records"""
    records = list(orchestrator._parse_input(BOLD_SAMPLE))
    record = next(r for r in records if r.id == "BHNHM001-24_r_1_s_50")

    # Test sequence basics
    seq_str = str(record.seq)
    assert len(seq_str) > 0, "Sequence should not be empty"
    assert seq_str.startswith("CGAAAATGAATATACTCAACAA"), "Sequence start mismatch"

    # Verify sequence characters
    valid_chars = set('ACGTN-')
    assert all(c.upper() in valid_chars for c in seq_str), "Invalid characters in sequence"

if __name__ == '__main__':
    pytest.main(['-v'])