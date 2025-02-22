import pytest
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from nbitk.config import Config

from barcode_validator.constants import TaxonomicBackbone, TaxonomicRank
from barcode_validator.dna_analysis_result import DNAAnalysisResult
from barcode_validator.resolvers.factory import ResolverFactory
from barcode_validator.resolvers.taxonomy import TaxonResolver

# Test data path handling
TEST_DATA_DIR = Path(__file__).parent / "data"
NSR_SAMPLE = TEST_DATA_DIR / "csc_sample.tsv"
NSR_ARCHIVE = TEST_DATA_DIR / "nsr-20250207.dwca.zip"
NCBI_TAXDUMP = TEST_DATA_DIR / "taxdump.tar.gz"


@pytest.fixture
def config():
    """Provides base configuration for testing"""
    conf = Config()
    conf.config_data = {
        'log_level': 'DEBUG',
        'level': 'family',
        'dwc_archive': str(NSR_ARCHIVE),
        'taxonomic_backbone': 'dwc',
        'ncbi_taxonomy': str(NCBI_TAXDUMP),
        'ncbi_taxonomy_url': 'https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz',
        'entrez_email': 'bioinformatics@naturalis.nl'
    }
    conf.initialized = True
    return conf


@pytest.fixture
def taxonomy_resolver(config):
    """Provides initialized TaxonomyResolver"""
    resolver = ResolverFactory.create_resolver(config, TaxonomicBackbone.NSR)
    resolver.load_tree(NSR_ARCHIVE)
    return resolver

def test_nsr_query_extraction(taxonomy_resolver):
    """Test extraction of identifications from NSR-style records"""
    # Test with species name
    record = SeqRecord(
        Seq(""),
        id="test_id",
        annotations={'bcdm_fields': {'identification': 'Homo sapiens'}}
    )
    assert taxonomy_resolver.parse_id(record) == "Homo sapiens"


def test_nsr_enrichment(taxonomy_resolver):
    """Test complete NSR taxonomy enrichment process"""
    # Create test record with DwC-style fields
    record = SeqRecord(
        Seq(""),
        id="test_id",
        annotations={'bcdm_fields': {
            'identification': 'Geophilus carpophagus',
            'marker_code': 'COI-5P'
        }}
    )
    result = DNAAnalysisResult(record.id)
    result.add_ancillary('marker_code', 'COI-5P')
    result.level = "family"

    # Perform enrichment
    taxonomy_resolver.enrich_result(record, result, TaxonomicRank.FAMILY)

    # Verify result population
    assert result.error is None, "Enrichment should succeed"
    assert result.species is not None, "Species should be populated"
    assert result.exp_taxon is not None, "Expected taxon should be populated"
    assert result.level == TaxonomicRank.FAMILY.value, "Expected to have family level"


def test_nsr_enrichment_errors(taxonomy_resolver):
    """Test error handling in NSR enrichment"""
    # Test missing identification
    record = SeqRecord(
        Seq(""),
        id="test_id",
        annotations={'bcdm_fields': {}}
    )
    result = DNAAnalysisResult(record.id)
    taxonomy_resolver.enrich_result(record, result, TaxonomicRank.FAMILY)
    assert result.error is not None
    assert "Could not parse ID" in result.error

    # Test invalid identification format
    record = SeqRecord(
        Seq(""),
        id="test_id",
        annotations={'bcdm_fields': {'identification': ''}}
    )
    result = DNAAnalysisResult(record.id)
    taxonomy_resolver.enrich_result(record, result, TaxonomicRank.FAMILY)
    assert result.error is not None
    assert "Could not parse ID" in result.error


def test_different_validation_levels(config, taxonomy_resolver):
    """Test enrichment with different validation levels"""
    # Test enrichment with order-level validation
    record = SeqRecord(
        Seq(""),
        id="test_id",
        annotations={'bcdm_fields': {
            'identification': 'Geophilus carpophagus',
            'marker_code': 'COI-5P'
        }}
    )
    result = DNAAnalysisResult(record.id)
    taxonomy_resolver.enrich_result(record, result, TaxonomicRank.ORDER)

    assert result.error is None, "Enrichment should succeed"
    assert result.level == "order", "Validation level should be order"
    assert result.exp_taxon.taxonomic_rank.lower() == "order"

def test_nsr_taxonomy_traversal(taxonomy_resolver):
    """Test taxonomy traversal in NSR backbone"""
    record = SeqRecord(
        Seq(""),
        id="test_id",
        annotations={'bcdm_fields': {
            'identification': 'Geophilus carpophagus',
            'marker_code': 'COI-5P'
        }}
    )
    result = DNAAnalysisResult(record.id)
    taxonomy_resolver.enrich_result(record, result, TaxonomicRank.FAMILY)

    # Test that we can traverse up to family
    assert result.exp_taxon is not None
    assert result.exp_taxon.taxonomic_rank.lower() == "family"

    # Verify the complete path exists
    ranks = taxonomy_resolver.get_lineage_dict(result.species)
    assert "kingdom" in ranks
    assert "phylum" in ranks
    assert "class" in ranks
    assert "order" in ranks
    assert "family" in ranks
    assert "genus" in ranks
    assert "species" in ranks

if __name__ == '__main__':
    pytest.main()