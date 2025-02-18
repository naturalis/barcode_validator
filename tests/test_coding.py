import pytest
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from barcode_validator.orchestrator import ValidationOrchestrator
from barcode_validator.validators.protein_coding import ProteinCodingValidator
from barcode_validator.dna_analysis_result import DNAAnalysisResult
from barcode_validator.taxonomy_resolver import TaxonomyResolver, Marker, TaxonomicBackbone
from nbitk.config import Config

# Test data path handling
TEST_DATA_DIR = Path(__file__).parent / "data"
BOLD_SAMPLE = TEST_DATA_DIR / "BGE00146_MGE-BGE_r1_1.3_1.5_s50_100.fasta"
BOLD_SHEET = TEST_DATA_DIR / "bold.xlsx"  # We'll need to create this
HMM_DIR = TEST_DATA_DIR / "hmm_profiles"  # We'll need to create this
NCBI_TAXDUMP = TEST_DATA_DIR / "taxdump.tar.gz"

@pytest.fixture
def config():
    """Test fixture that provides a Config object configured for COI-5P validation"""
    conf = Config()
    conf.config_data = {
        'log_level': 'DEBUG',
        'validate_taxonomy': False,  # No BLASTing just yet
        'validate_structure': True,  # Enable structural validation
        'marker': 'COI-5P',  # Using COI-5P marker
        'hmm_profile_dir': str(HMM_DIR),
        'bold_sheet_file': str(BOLD_SHEET),
        'taxonomic_backbone': 'bold',
        'ncbi_taxonomy': str(NCBI_TAXDUMP),
        'ncbi_taxonomy_url': 'https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'

    }
    conf.initialized = True
    return conf


@pytest.fixture
def taxonomy_resolver(config):
    """Fixture providing TaxonomyResolver instance"""
    resolver = TaxonomyResolver(config)
    resolver.setup_taxonomy(TaxonomicBackbone.BOLD)  # This will load BOLD taxonomy
    return resolver


@pytest.fixture
def validator(config, taxonomy_resolver):
    """Fixture providing ProteinCodingValidator instance"""
    return ProteinCodingValidator(config, str(HMM_DIR))


@pytest.fixture
def fasta_records():
    """Fixture providing list of records from the FASTA file"""
    return list(SeqIO.parse(BOLD_SAMPLE, "fasta"))


@pytest.fixture
def orchestrator(config):
    """Fixture providing orchestrator instance"""
    return ValidationOrchestrator(config)


def test_validator_initialization(validator):
    """Test basic validator initialization"""
    assert validator.marker == Marker.COI_5P
    assert validator.config is not None
    assert validator.logger is not None
    assert validator.hmm_profile_dir == HMM_DIR


def test_translation_table_determination(validator, taxonomy_resolver, fasta_records):
    """Test translation table determination for COI-5P sequences"""
    record = fasta_records[0]  # Use first record

    # Process ID format is BHNHM001-24_r_1_s_50, we want BHNHM001-24
    process_id = record.id.split('_')[0]

    # Resolve backbone taxonomy
    backbone_taxon = taxonomy_resolver.resolve_backbone(process_id, TaxonomicBackbone.BOLD)
    assert backbone_taxon is not None, "Failed to resolve backbone taxonomy"

    # Get validation level taxon (family)
    ncbi_valid = taxonomy_resolver.find_taxon_at_level(backbone_taxon, 'family')
    assert ncbi_valid is not None, "Failed to get NCBI taxon"

    # Get translation table
    trans_table = taxonomy_resolver.get_translation_table(Marker.COI_5P, ncbi_valid)
    assert isinstance(trans_table, int), "Translation table must be integer"
    assert trans_table in [2, 5], "Expected vertebrate (2) or invertebrate (5) mitochondrial code"


def test_hmm_alignment(validator, fasta_records):
    """Test HMM alignment of COI-5P sequences"""
    record = fasta_records[0]
    result = DNAAnalysisResult(record.id)

    # Test alignment process
    aligned_seq = validator._align_sequence(record)
    assert aligned_seq is not None, "HMM alignment failed"
    assert len(aligned_seq.seq) > 0, "Aligned sequence is empty"

    # Check alignment output format
    assert isinstance(aligned_seq, SeqRecord), "Alignment output should be SeqRecord"
    assert str(aligned_seq.seq).count('-') >= 0, "Alignment may contain gaps"


def test_reading_frame_detection(validator, fasta_records):
    """Test reading frame detection in aligned sequences"""
    record = fasta_records[0]
    result = DNAAnalysisResult(record.id)

    # First align sequence
    aligned_seq = validator._align_sequence(record)
    assert aligned_seq is not None, "HMM alignment failed"

    # Detect reading frame
    frame = validator._determine_reading_frame(aligned_seq, 5)
    assert frame in [0, 1, 2], "Reading frame must be 0, 1, or 2"

    # Verify frame produces fewest stops
    seq_nogaps = str(aligned_seq.seq).replace('-', '')
    for test_frame in range(3):
        if test_frame != frame:
            test_seq = Seq(seq_nogaps[test_frame:])
            test_stops = test_seq.translate(table=1).count('*')
            frame_seq = Seq(seq_nogaps[frame:])
            frame_stops = frame_seq.translate(table=1).count('*')
            assert frame_stops <= test_stops, "Selected frame does not minimize stop codons"


def test_stop_codon_detection(validator, fasta_records):
    """Test stop codon detection in COI-5P sequences"""
    record = fasta_records[0]
    result = DNAAnalysisResult(record.id)

    # First align sequence
    aligned_seq = validator._align_sequence(record)
    assert aligned_seq is not None, "HMM alignment failed"

    # Get reading frame
    frame = validator._determine_reading_frame(aligned_seq, 5)

    # Find stop codons,
    stops = validator._find_stop_codons(aligned_seq, 5, frame)  # Using invertebrate mitochondrial code
    assert isinstance(stops, list), "Stop codon positions should be list"
    assert all(isinstance(pos, int) for pos in stops), "Stop positions should be integers"

    # Verify stop positions
    seq_nogaps = str(aligned_seq.seq).replace('-', '')
    for pos in stops:
        codon_start = pos - (pos % 3)
        codon = seq_nogaps[codon_start:codon_start + 3]
        assert len(codon) == 3, f"Invalid codon position {pos}"
        assert Seq(codon).translate(table=5) == '*', f"Position {pos} is not a stop codon"


def test_complete_validation(validator, fasta_records):
    """Test complete protein-coding validation process"""
    for record in fasta_records:
        result = DNAAnalysisResult(record.id)
        result.add_ancillary('translation_table', '5')  # Using invertebrate mitochondrial code
        result.add_ancillary('marker_code', 'COI-5P')

        # Perform validation
        validator.validate_marker_specific(record, result)

        # Check results
        assert 'reading_frame' in result.ancillary, "Reading frame not determined"
        assert isinstance(result.stop_codons, list), "Stop codons not detected"
        frame = int(result.ancillary['reading_frame'])
        assert frame in [0, 1, 2], "Invalid reading frame"