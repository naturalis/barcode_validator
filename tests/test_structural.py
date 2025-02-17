import pytest
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from barcode_validator.orchestrator import ValidationOrchestrator
from barcode_validator.validators.non_coding import NonCodingValidator
from barcode_validator.dna_analysis_result import DNAAnalysisResult
from nbitk.config import Config
from barcode_validator.taxonomy_resolver import Marker

# Test data path handling
TEST_DATA_DIR = Path(__file__).parent / "data"
BOLD_SAMPLE = TEST_DATA_DIR / "BGE00146_MGE-BGE_r1_1.3_1.5_s50_100.fasta"


@pytest.fixture
def config():
    """Test fixture that provides a Config object configured for ITS structural validation"""
    conf = Config()
    conf.config_data = {
        'log_level': 'DEBUG',
        'validate_taxonomy': False,  # No taxonomy needed
        'validate_structure': True,  # Enable structural validation
        'marker': 'ITS'  # Using ITS marker
    }
    conf.initialized = True
    return conf


@pytest.fixture
def orchestrator(config):
    """Fixture providing orchestrator instance"""
    return ValidationOrchestrator(config)


@pytest.fixture
def validator(config):
    """Fixture providing NonCodingValidator instance"""
    return NonCodingValidator(config)


@pytest.fixture
def fasta_records():
    """Fixture providing list of records from the FASTA file"""
    return list(SeqIO.parse(BOLD_SAMPLE, "fasta"))


def test_validator_initialization(validator):
    """Test basic validator initialization"""
    assert validator.marker == Marker.ITS
    assert validator.config is not None
    assert validator.logger is not None


def test_sequence_validation_basic(validator, fasta_records):
    """Test basic sequence validation functionality using real FASTA records"""
    record = fasta_records[0]  # Use first record from FASTA
    result = DNAAnalysisResult(record.id)

    validator.validate_sequence(record, result)

    # Check basic metrics
    assert result.seq_length > 0
    assert isinstance(result.ambiguities, int)
    assert isinstance(result.full_ambiguities, int)


def test_ambiguous_base_counting(validator, fasta_records):
    """Test counting of ambiguous bases in real sequences"""
    for record in fasta_records:
        result = DNAAnalysisResult(record.id)
        validator.validate_sequence(record, result)

        # Count ambiguous bases manually for verification
        seq_str = str(record.seq).upper()
        manual_ambig_count = sum(1 for base in seq_str if base not in 'ACGT-~')

        assert result.ambiguities == manual_ambig_count
        assert result.full_ambiguities == manual_ambig_count


def test_gap_handling(validator, fasta_records):
    """Test handling of gap characters in sequences"""
    for record in fasta_records:
        result = DNAAnalysisResult(record.id)
        validator.validate_sequence(record, result)

        seq_str = str(record.seq)
        gap_count = seq_str.count('-') + seq_str.count('~')

        # Length should exclude gaps
        assert result.seq_length == len(seq_str) - gap_count
        # Gaps should not count as ambiguities
        assert result.ambiguities >= 0


def test_real_fasta_validation(orchestrator):
    """Test validation with complete FASTA file"""
    results = orchestrator.validate_file(BOLD_SAMPLE)

    assert len(results.results) == 10  # We know there are 10 records

    # Check all results
    for result in results.results:
        assert result.seq_length > 0
        assert isinstance(result.ambiguities, int)
        assert isinstance(result.full_ambiguities, int)
        assert 'gc_content' in result.ancillary
        gc_content = float(result.ancillary['gc_content'])
        assert 0 <= gc_content <= 100


def test_gc_content_calculation(validator, fasta_records):
    """Test GC content calculation with real sequences"""
    for record in fasta_records:
        result = DNAAnalysisResult(record.id)
        validator.validate_sequence(record, result)

        # Calculate GC content manually for verification
        seq_str = str(record.seq).upper()
        gc_count = seq_str.count('G') + seq_str.count('C')
        total = len(seq_str) - seq_str.count('-') - seq_str.count('~')
        expected_gc = (gc_count / total * 100) if total > 0 else 0

        assert abs(float(result.ancillary['gc_content']) - expected_gc) < 0.01


def test_edge_cases(validator, caplog):
    """Test edge cases not present in the FASTA file"""
    # Empty sequence
    empty_record = SeqRecord(seq=Seq(""), id="empty_seq")
    empty_result = DNAAnalysisResult("empty_seq")
    validator.validate_sequence(empty_record, empty_result)
    assert empty_result.seq_length == 0

    # Invalid characters
    invalid_seq = "ACGT" * 10 + "RYNX"
    invalid_record = SeqRecord(seq=Seq(invalid_seq), id="invalid_seq")
    invalid_result = DNAAnalysisResult("invalid_seq")
    validator.validate_sequence(invalid_record, invalid_result)
    assert invalid_result.ambiguities == 3  # RYN are ambiguous, X is not