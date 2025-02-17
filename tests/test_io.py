import pytest
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from barcode_validator.orchestrator import ValidationOrchestrator
from nbitk.config import Config

# Test data path handling
TEST_DATA_DIR = Path(__file__).parent / "data"
CSC_SAMPLE = TEST_DATA_DIR / "csc_sample.txt"


@pytest.fixture
def config():
    """Test fixture that provides a minimal Config object"""
    conf = Config()
    conf.config_data = {
        'log_level': 'DEBUG',
        'validate_taxonomy': False,  # No taxonomy needed
        'validate_structure': False  # No structure validation needed
    }
    conf.initialized = True
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
    assert bcdm_fields['rank'] == "Species"
    assert bcdm_fields['source'] == "ada"


def test_null_value_handling(orchestrator):
    """Test handling of records with null values"""
    records = list(orchestrator._parse_input(CSC_SAMPLE))
    null_record = next(r for r in records if r.id == "507")

    bcdm_fields = null_record.annotations.get('bcdm_fields', {})
    assert bcdm_fields['kingdom'] == "null"
    assert bcdm_fields['rank'] == "null"


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

    required_fields = {'marker_code', 'identification', 'kingdom', 'rank', 'source'}

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


def test_source_values(orchestrator):
    """Test the valid source values"""
    records = list(orchestrator._parse_input(CSC_SAMPLE))

    valid_sources = {'ada', 'sheet'}
    for record in records:
        source = record.annotations['bcdm_fields']['source']
        assert source in valid_sources, f"Invalid source value '{source}' in record {record.id}"


def test_invalid_file_handling(orchestrator):
    """Test handling of invalid file paths"""
    with pytest.raises(ValueError):
        list(orchestrator._parse_input(TEST_DATA_DIR / "nonexistent.xyz"))

def test_parsing_no_taxonomy_init(orchestrator):
    """Verify that parsing alone doesn't initialize taxonomy"""
    records = list(orchestrator._parse_input(CSC_SAMPLE))
    assert orchestrator.taxonomy_resolver.ncbi_tree is None
    assert orchestrator.taxonomy_resolver.backbone_tree is None