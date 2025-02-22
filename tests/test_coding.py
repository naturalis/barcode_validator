import pytest
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from nbitk.Taxon import Taxon
from nbitk.Tools import Hmmalign

from barcode_validator.orchestrator import ValidationOrchestrator
from barcode_validator.validators.factory import StructureValidatorFactory
from barcode_validator.validators.protein_coding import ProteinCodingValidator
from barcode_validator.resolvers.factory import ResolverFactory
from barcode_validator.resolvers.taxonomy import TaxonResolver
from barcode_validator.dna_analysis_result import DNAAnalysisResult
from barcode_validator.constants import Marker, TaxonomicBackbone, TaxonomicRank
from nbitk.config import Config

# Test data path handling
TEST_DATA_DIR = Path(__file__).parent / "data"
BOLD_SAMPLE = TEST_DATA_DIR / "BGE00146_MGE-BGE_r1_1.3_1.5_s50_100.fasta"
BOLD_SHEET = TEST_DATA_DIR / "bold.xlsx"  # We'll need to create this
HMM_DIR = TEST_DATA_DIR / "hmm_profiles"
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
    }
    conf.initialized = True
    return conf


@pytest.fixture
def bold_resolver(config) -> TaxonResolver:
    """Fixture providing TaxonomyResolver instance"""
    resolver = ResolverFactory.create_resolver(config, TaxonomicBackbone.BOLD)
    resolver.load_tree(BOLD_SHEET)
    return resolver

@pytest.fixture
def ncbi_resolver(config) -> TaxonResolver:
    """Fixture providing TaxonomyResolver instance"""
    resolver = ResolverFactory.create_resolver(config, TaxonomicBackbone.NCBI)
    resolver.load_tree(NCBI_TAXDUMP)
    return resolver

@pytest.fixture
def validator(config) -> ProteinCodingValidator:
    """Fixture providing ProteinCodingValidator instance"""
    protein_coding_validator : ProteinCodingValidator = StructureValidatorFactory.create_validator(config, Marker.COI_5P)
    protein_coding_validator.set_hmm_profile_dir(HMM_DIR)
    return protein_coding_validator

@pytest.fixture
def hmmalign(config) -> Hmmalign:
    """Fixture providing a configured Hmmalign instance"""
    hmmalign = Hmmalign(config)
    hmm_file = f"{str(HMM_DIR)}/{config.get('marker', 'COI-5P')}.hmm"
    hmmalign.set_hmmfile(hmm_file)
    return hmmalign


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
    assert validator.logger is not None
    assert validator.hmm_profile_dir == HMM_DIR

@pytest.mark.skipif(not Path(NCBI_TAXDUMP).exists(), reason="NCBI taxonomy does not exist")
def test_translation_table_determination(validator, bold_resolver, ncbi_resolver, fasta_records):
    """Test translation table determination for COI-5P sequences"""
    record = fasta_records[0]  # Use first record

    # Process ID format is BHNHM001-24_r_1_s_50, we want BHNHM001-24
    process_id = bold_resolver.parse_id(record)

    # Resolve backbone taxonomy
    bold_matches = bold_resolver.find_nodes(process_id)
    assert bold_matches is not None, "Failed to resolve backbone taxonomy"
    assert len(bold_matches) == 1, "Must find one match"
    assert isinstance(bold_matches[0], Taxon), "Match is a Taxon object"

    # Get validation level taxon (family)
    validation_taxon = bold_resolver.find_ancestor_at_rank(bold_matches[0], TaxonomicRank.FAMILY)
    assert validation_taxon is not None, "Failed to get validation-level taxon (i.e. the family)"

    # Get NCBI equivalent
    ncbi_matches = ncbi_resolver.match_nodes(validation_taxon)
    assert ncbi_matches is not None, "Failed to find NCBI equivalent"
    assert len(ncbi_matches) == 1, "Must find one NCBI match"
    ncbi_valid = ncbi_matches[0]

    # Get translation table by way of NCBI taxonomy, resolver should be set here
    validator.set_taxonomy_resolver(ncbi_resolver)
    trans_table = validator.get_translation_table(Marker.COI_5P, ncbi_valid)
    assert isinstance(trans_table, int), "Translation table must be integer"
    assert trans_table in [2, 5], "Expected vertebrate (2) or invertebrate (5) mitochondrial code"


def test_hmm_alignment(validator, fasta_records, hmmalign: Hmmalign):
    """Test HMM alignment of COI-5P sequences"""
    record = fasta_records[0]
    result = DNAAnalysisResult(record.id)
    validator.set_hmmalign(hmmalign)
    validator.set_hmm_profile_dir(HMM_DIR)
    validator.set_marker(Marker.COI_5P)

    # Test alignment process
    aligned_seq = validator._align_sequence(record)
    assert aligned_seq is not None, "HMM alignment failed"
    assert len(aligned_seq.seq) > 0, "Aligned sequence is empty"

    # Check alignment output format
    assert isinstance(aligned_seq, SeqRecord), "Alignment output should be SeqRecord"
    assert str(aligned_seq.seq).count('-') >= 0, "Alignment may contain gaps"

@pytest.mark.skipif(not Path(NCBI_TAXDUMP).exists(), reason="NCBI taxonomy does not exist")
def test_reading_frame_detection(validator, fasta_records, hmmalign: Hmmalign, config: Config):
    """Test reading frame detection in aligned sequences"""
    record = fasta_records[0]
    result = DNAAnalysisResult(record.id)
    validator.set_hmmalign(hmmalign)
    validator.set_hmm_profile_dir(HMM_DIR)
    validator.set_marker(Marker.COI_5P)

    # First align sequence
    aligned_seq = validator._align_sequence(record)
    assert aligned_seq is not None, "HMM alignment failed"

    # Detect reading frame
    frame = validator._determine_reading_frame(aligned_seq, 5)
    assert frame == 1, "Most records are offset by 1 with this HMM"


def test_stop_codon_detection(validator, fasta_records, hmmalign: Hmmalign, config: Config):
    """Test stop codon detection in COI-5P sequences"""
    record = fasta_records[0]
    result = DNAAnalysisResult(record.id)
    validator.set_hmmalign(hmmalign)
    validator.set_hmm_profile_dir(HMM_DIR)
    validator.set_marker(Marker.COI_5P)

    # First align sequence
    aligned_seq = validator._align_sequence(record)
    assert aligned_seq is not None, "HMM alignment failed"

    # Get reading frame
    frame, protein = validator._determine_reading_frame(aligned_seq, 5)

    # Find stop codons,
    stops = validator._find_stop_codons(frame, protein)  # Using invertebrate mitochondrial code
    assert isinstance(stops, list), "Stop codon positions should be list"
    assert all(isinstance(pos, int) for pos in stops), "Stop positions should be integers"
    assert len(stops) == 0, "Valid sequence should not have stop codons"


def test_complete_validation(validator, fasta_records, bold_resolver, hmmalign: Hmmalign, config: Config):
    """Test complete protein-coding validation process"""
    validator.set_hmmalign(hmmalign)
    validator.set_hmm_profile_dir(HMM_DIR)
    validator.set_marker(Marker.COI_5P)
    validator.set_taxonomy_resolver(bold_resolver)

    for record in fasta_records:
        result = DNAAnalysisResult(record.id)
        result.add_ancillary('translation_table', '5')  # Using invertebrate mitochondrial code
        result.add_ancillary('marker_code', 'COI-5P')

        # Perform validation
        bold_resolver.enrich_result(record, result, TaxonomicRank.FAMILY)
        validator.validate_marker_specific(record, result)

        # Check results
        assert 'reading_frame' in result.ancillary, "Reading frame not determined"
        assert isinstance(result.stop_codons, list), "Stop codons not detected"
        frame = int(result.ancillary['reading_frame'])
        assert frame in [0, 1, 2], "Invalid reading frame"