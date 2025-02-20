import pytest
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from nbitk.config import Config
from barcode_validator.taxonomy_resolver import TaxonomyResolver
from barcode_validator.taxonomic_enrichment import (
    NSRTaxonomyParser,
    TaxonomicBackbone
)
from barcode_validator.dna_analysis_result import DNAAnalysisResult

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
    resolver = TaxonomyResolver(config)
    resolver.setup_taxonomy(TaxonomicBackbone.DWC)
    return resolver


@pytest.fixture
def nsr_parser(config, taxonomy_resolver):
    """Provides NSR taxonomy parser"""
    return NSRTaxonomyParser(config, taxonomy_resolver)


def test_nsr_query_extraction(nsr_parser):
    """Test extraction of identifications from NSR-style records"""
    # Test with species name
    record = SeqRecord(
        Seq(""),
        id="test_id",
        annotations={'bcdm_fields': {'identification': 'Homo sapiens'}}
    )
    assert nsr_parser.get_query(record) == "Homo sapiens"

    # Test with only genus
    record = SeqRecord(
        Seq(""),
        id="test_id",
        annotations={'bcdm_fields': {'identification': 'Homo'}}
    )
    assert nsr_parser.get_query(record) == "Homo"


def test_nsr_enrichment(nsr_parser):
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
    nsr_parser.enrich_result(record, result, TaxonomicBackbone.DWC)

    # Verify result population
    assert result.error is None, "Enrichment should succeed"
    assert result.species is not None, "Species should be populated"
    assert result.exp_taxon is not None, "Expected taxon should be populated"
    assert int(result.ancillary.get('translation_table')) == 5, "Should have inferred insect mitochondrial"


def test_nsr_enrichment_errors(nsr_parser):
    """Test error handling in NSR enrichment"""
    # Test missing identification
    record = SeqRecord(
        Seq(""),
        id="test_id",
        annotations={'bcdm_fields': {}}
    )
    result = DNAAnalysisResult(record.id)
    nsr_parser.enrich_result(record, result, TaxonomicBackbone.DWC)
    assert result.error is not None
    assert "No identification provided" in result.error

    # Test invalid identification format
    record = SeqRecord(
        Seq(""),
        id="test_id",
        annotations={'bcdm_fields': {'identification': ''}}
    )
    result = DNAAnalysisResult(record.id)
    nsr_parser.enrich_result(record, result, TaxonomicBackbone.DWC)
    assert result.error is not None
    assert "No entry found" in result.error


def test_different_validation_levels(config, taxonomy_resolver):
    """Test enrichment with different validation levels"""
    # Modify config for order-level validation
    config.config_data['validation_rank'] = 'order'
    parser = NSRTaxonomyParser(config, taxonomy_resolver)

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
    result.add_ancillary('marker_code', 'COI-5P')
    result.level = "order"
    parser.enrich_result(record, result, TaxonomicBackbone.DWC)

    assert result.error is None, "Enrichment should succeed"
    assert result.level == "order", "Validation level should be order"
    assert result.exp_taxon.taxonomic_rank.lower() == "order"


def test_backbone_mismatches(nsr_parser):
    """Test handling of backbone type mismatches"""
    record = SeqRecord(
        Seq(""),
        id="test_id",
        annotations={'bcdm_fields': {
            'identification': 'Geophilus carpophagus',
            'marker_code': 'COI-5P'
        }}
    )
    result = DNAAnalysisResult(record.id)

    # Test with wrong backbone type
    with pytest.raises(SystemExit):
        nsr_parser.enrich_result(record, result, TaxonomicBackbone.BOLD)


def test_nsr_taxonomy_traversal(nsr_parser):
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
    result.add_ancillary('marker_code', 'COI-5P')
    result.level = "family"

    nsr_parser.enrich_result(record, result, TaxonomicBackbone.DWC)

    # Test that we can traverse up to family
    assert result.exp_taxon is not None
    assert result.exp_taxon.taxonomic_rank.lower() == "family"

    # Verify the complete path exists
    path = nsr_parser.taxonomy_resolver.backbone_tree.root.get_path(result.species)
    ranks = [node.taxonomic_rank.lower() for node in path]
    assert "kingdom" in ranks
    assert "phylum" in ranks
    assert "class" in ranks
    assert "order" in ranks
    assert "family" in ranks
    assert "genus" in ranks
    assert "species" in ranks