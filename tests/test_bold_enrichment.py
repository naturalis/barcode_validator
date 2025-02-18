import pytest
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from nbitk.config import Config
from barcode_validator.taxonomy_resolver import TaxonomyResolver
from barcode_validator.taxonomic_enrichment import (
    BOLDTaxonomyParser,
    NSRTaxonomyParser,
    TaxonomicBackbone
)
from barcode_validator.dna_analysis_result import DNAAnalysisResult

# Test data path handling
TEST_DATA_DIR = Path(__file__).parent / "data"
BOLD_SAMPLE = TEST_DATA_DIR / "BGE00146_MGE-BGE_r1_1.3_1.5_s50_100.fasta"
BOLD_SHEET = TEST_DATA_DIR / "bold.xlsx"
NCBI_TAXDUMP = TEST_DATA_DIR / "taxdump.tar.gz"


@pytest.fixture
def config():
    """Provides base configuration for testing"""
    conf = Config()
    conf.config_data = {
        'log_level': 'DEBUG',
        'level': 'family',
        'bold_sheet_file': str(BOLD_SHEET),
        'taxonomic_backbone': 'bold',
        'ncbi_taxonomy': str(NCBI_TAXDUMP),
        'ncbi_taxonomy_url': 'https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'
    }
    conf.initialized = True
    return conf

@pytest.fixture
def taxonomy_resolver(config):
    """Provides initialized TaxonomyResolver"""
    resolver = TaxonomyResolver(config)
    resolver.setup_taxonomy(TaxonomicBackbone.BOLD)
    return resolver

@pytest.fixture
def bold_parser(config, taxonomy_resolver):
    """Provides BOLD taxonomy parser"""
    return BOLDTaxonomyParser(config, taxonomy_resolver)

def test_bold_query_extraction(bold_parser):
    """Test extraction of process IDs from BOLD-style records"""
    # Test normal case
    record = SeqRecord(Seq(""), id="BHNHM001-24_r_1_s_50")
    assert bold_parser.get_query(record) == "BHNHM001-24"

    # Test record without suffixes
    record = SeqRecord(Seq(""), id="BHNHM001-24")
    assert bold_parser.get_query(record) == "BHNHM001-24"

def test_bold_enrichment(bold_parser):
    """Test complete BOLD taxonomy enrichment process"""
    # Create test record
    record = SeqRecord(Seq(""), id="BHNHM001-24_r_1_s_50")
    result = DNAAnalysisResult(record.id)
    result.add_ancillary('marker_code', 'COI-5P')
    result.level = "family"

    # Perform enrichment
    bold_parser.enrich_result(record, result, TaxonomicBackbone.BOLD)

    # Verify result population
    assert result.error is None, "Enrichment should succeed"
    assert result.species is not None, "Species should be populated"
    assert result.exp_taxon is not None, "Expected taxon should be populated"
    assert int(result.ancillary.get('translation_table')) == 5, "Should have inferred insect mitochondrial"

def test_bold_enrichment_errors(bold_parser):
    """Test error handling in BOLD enrichment"""
    # Test invalid process ID format
    record = SeqRecord(Seq(""), id="invalid_id")
    result = DNAAnalysisResult(record.id)
    bold_parser.enrich_result(record, result, TaxonomicBackbone.BOLD)
    assert result.error is not None
    assert "No entry found" in result.error

    # Test missing process ID in taxonomy
    record = SeqRecord(Seq(""), id="XXXXX999-99_r_1")
    result = DNAAnalysisResult(record.id)
    bold_parser.enrich_result(record, result, TaxonomicBackbone.BOLD)
    assert result.error is not None
    assert "No entry found" in result.error

def test_taxonomy_initialization():
    """Test parser taxonomy backbone initialization"""
    config = Config()
    config.config_data = {'log_level': 'DEBUG'}
    config.initialized = True
    resolver = TaxonomyResolver(config)

    # Test BOLD parser
    bold = BOLDTaxonomyParser(config, resolver)
    assert bold.taxonomy == TaxonomicBackbone.BOLD

    # Test NSR parser
    nsr = NSRTaxonomyParser(config, resolver)
    assert nsr.taxonomy == TaxonomicBackbone.DWC

def test_different_validation_levels(config, taxonomy_resolver):
    """Test enrichment with different validation levels"""
    # Modify config for order-level validation
    config.config_data['validation_rank'] = 'order'
    parser = BOLDTaxonomyParser(config, taxonomy_resolver)

    # Test enrichment
    record = SeqRecord(Seq(""), id="BHNHM001-24_r_1_s_50")
    result = DNAAnalysisResult(record.id)
    result.add_ancillary('marker_code', 'COI-5P')
    result.level = "order"
    parser.enrich_result(record, result, TaxonomicBackbone.BOLD)

    assert result.error is None, "Enrichment should succeed"
    assert result.level == "order", "Validation level should be order"
    assert result.exp_taxon.taxonomic_rank.lower() == "order"