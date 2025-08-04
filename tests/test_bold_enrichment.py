import pytest
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from nbitk.config import Config
from barcode_validator.resolvers.factory import ResolverFactory
from barcode_validator.dna_analysis_result import DNAAnalysisResult
from barcode_validator.resolvers.taxonomy import TaxonomicBackbone, TaxonomicRank
from barcode_validator.config.schema_config import SchemaConfig

# Test data path handling
TEST_DATA_DIR = Path(__file__).parent / "data"
BOLD_SAMPLE = TEST_DATA_DIR / "BGE00146_MGE-BGE_r1_1.3_1.5_s50_100.fasta"
BOLD_SHEET = TEST_DATA_DIR / "bold.xlsx"
NCBI_TAXDUMP = TEST_DATA_DIR / "taxdump.tar.gz"


@pytest.fixture
def config():
    """Provides base configuration for testing"""
    conf = SchemaConfig()
    conf.set('log_level', 'DEBUG')
    conf.set('taxon_validation.rank', 'family')
    conf.set('exp_taxonomy', BOLD_SHEET)
    conf.set('exp_taxonomy_type', 'bold')
    conf.set('reflib_taxonomy', NCBI_TAXDUMP)
    return conf

@pytest.fixture
def taxonomy_resolver(config):
    """Provides initialized TaxonomyResolver"""
    resolver = ResolverFactory.create_resolver(config, TaxonomicBackbone.BOLD)
    resolver.load_tree(BOLD_SHEET)
    return resolver

def test_bold_query_extraction(taxonomy_resolver):
    """Test extraction of process IDs from BOLD-style records"""
    # Test normal case
    record = SeqRecord(Seq(""), id="BHNHM001-24_r_1_s_50")
    assert taxonomy_resolver.parse_id(record) == "BHNHM001-24"

    # Test record without suffixes
    record = SeqRecord(Seq(""), id="BHNHM001-24")
    assert taxonomy_resolver.parse_id(record) == "BHNHM001-24"

def test_bold_enrichment(taxonomy_resolver):
    """Test complete BOLD taxonomy enrichment process"""
    # Create test record
    record = SeqRecord(Seq(""), id="BHNHM001-24_r_1_s_50")
    result = DNAAnalysisResult(record.id)

    # Perform enrichment
    taxonomy_resolver.enrich_result(record, result, TaxonomicRank.FAMILY)

    # Verify result population
    assert result.error is None, "Enrichment should succeed"
    assert result.species is not None, "Species should be populated"
    assert result.exp_taxon is not None, "Expected taxon should be populated"
    assert result.level == TaxonomicRank.FAMILY.value, "Expected to have family level"

def test_bold_enrichment_errors(taxonomy_resolver):
    """Test error handling in BOLD enrichment"""
    # Test invalid process ID format
    record = SeqRecord(Seq(""), id="invalid_id")
    result = DNAAnalysisResult(record.id)
    taxonomy_resolver.enrich_result(record, result, TaxonomicRank.FAMILY)
    assert result.error is not None
    assert "Could not find nodes" in result.error

    # Test missing process ID in taxonomy
    record = SeqRecord(Seq(""), id="XXXXX999-99_r_1")
    result = DNAAnalysisResult(record.id)
    taxonomy_resolver.enrich_result(record, result, TaxonomicRank.FAMILY)
    assert result.error is not None
    assert "Could not find nodes" in result.error

def test_taxonomy_initialization(taxonomy_resolver):
    """Test parser taxonomy backbone initialization"""
    assert taxonomy_resolver.get_type() == TaxonomicBackbone.BOLD

def test_different_validation_levels(taxonomy_resolver):
    """Test enrichment with different validation levels"""

    # Test enrichment
    record = SeqRecord(Seq(""), id="BHNHM001-24_r_1_s_50")
    result = DNAAnalysisResult(record.id)
    taxonomy_resolver.enrich_result(record, result, TaxonomicRank.ORDER)

    assert result.error is None, "Enrichment should succeed"
    assert result.level.lower() == "order", "Validation level should be order"
    assert result.exp_taxon.taxonomic_rank.lower() == "order"

if __name__ == "__main__":
    pytest.main()