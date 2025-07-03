import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from nbitk.config import Config
from barcode_validator.constants import TaxonomicRank, RefDB
from barcode_validator.idservices.bold import BOLD
from barcode_validator.idservices.factory import IDServiceFactory

@pytest.fixture()
def config():
    """Fixture to provide a Config object for tests."""
    config = Config()
    config.config_data = { "log_level": "DEBUG" }
    config.initialized = True
    return config

class TestBOLDInstantiation:
    """Test BOLD service instantiation and configuration."""

    def test_direct_instantiation(self, config):
        """Test direct instantiation of BOLD service."""
        bold = BOLD(config)

        assert bold is not None
        assert isinstance(bold, BOLD)
        assert bold.base_url == "https://id.boldsystems.org"
        assert bold.default_database == "all.tax-derep"
        assert bold.default_min_identity == 0.94
        assert bold.default_max_hits == 25

    def test_custom_base_url_instantiation(self, config):
        """Test instantiation with custom base URL."""
        custom_url = "https://custom.bold.org"
        bold = BOLD(config, base_url=custom_url)

        assert bold.base_url == custom_url

    def test_factory_instantiation(self, config):
        """Test instantiation through IDServiceFactory."""
        service = IDServiceFactory.create_idservice(config, RefDB.BOLD)

        assert service is not None
        assert isinstance(service, BOLD)

    def test_requires_resolver_static_method(self):
        """Test requires_resolver static method."""
        assert BOLD.requires_resolver() is True

    def test_requires_blastn_static_method(self):
        """Test requires_blastn static method."""
        assert BOLD.requires_blastn() is False


class TestBOLDConfiguration:
    """Test BOLD service configuration methods."""

    def test_set_search_parameters(self, config):
        """Test setting search parameters."""
        bold = BOLD(config)

        bold.set_search_parameters(
            database="public.tax-derep",
            min_identity=0.90,
            min_overlap=200,
            max_hits=50
        )

        assert bold.default_database == "public.tax-derep"
        assert bold.default_min_identity == 0.90
        assert bold.default_min_overlap == 200
        assert bold.default_max_hits == 50

    def test_set_partial_search_parameters(self, config):
        """Test setting only some search parameters."""
        bold = BOLD(config)
        original_database = bold.default_database

        bold.set_search_parameters(min_identity=0.85)

        assert bold.default_database == original_database  # unchanged
        assert bold.default_min_identity == 0.85  # changed

    def test_get_available_databases(self, config):
        """Test getting available databases."""
        bold = BOLD(config)
        databases = bold.get_available_databases()

        assert isinstance(databases, dict)
        assert "all.tax-derep" in databases
        assert "public.tax-derep" in databases
        assert "public.plants" in databases
        assert "public.fungi" in databases
        assert len(databases) >= 8


class TestBOLDValidation:
    """Test BOLD service validation methods."""

    def test_valid_fasta_single_sequence(self, config):
        """Test validation of valid single-sequence FASTA."""
        bold = BOLD(config)
        fasta = ">seq1\nATCGATCGATCG"

        assert bold._is_valid_fasta(fasta) is True

    def test_valid_fasta_multiple_sequences(self, config):
        """Test validation of valid multi-sequence FASTA."""
        bold = BOLD(config)
        fasta = ">seq1\nATCGATCG\n>seq2\nGGCCTTAA"

        assert bold._is_valid_fasta(fasta) is True

    def test_valid_fasta_with_gaps(self, config):
        """Test validation of FASTA with gaps and ambiguous bases."""
        bold = BOLD(config)
        fasta = ">seq1\nATCG-ATCGN"

        assert bold._is_valid_fasta(fasta) is True

    def test_invalid_fasta_no_header(self, config):
        """Test validation fails for FASTA without header."""
        bold = BOLD(config)
        fasta = "ATCGATCGATCG"

        assert bold._is_valid_fasta(fasta) is False

    def test_invalid_fasta_no_sequence(self, config):
        """Test validation fails for FASTA without sequence."""
        bold = BOLD(config)
        fasta = ">seq1"

        assert bold._is_valid_fasta(fasta) is False

    def test_invalid_fasta_empty(self, config):
        """Test validation fails for empty string."""
        bold = BOLD(config)
        fasta = ""

        assert bold._is_valid_fasta(fasta) is False

    def test_invalid_fasta_invalid_characters(self, config):
        """Test validation fails for invalid DNA characters."""
        bold = BOLD(config)
        fasta = ">seq1\nATCGXYZ"

        assert bold._is_valid_fasta(fasta) is False


class TestBOLDSequenceHandling:
    """Test BOLD service sequence handling methods."""

    def test_seqrecord_to_fasta_conversion(self, config):
        """Test conversion of SeqRecord to FASTA format."""
        bold = BOLD(config)

        seq = Seq("ATCGATCGATCG")
        record = SeqRecord(seq, id="test_seq", description="Test sequence")

        fasta = bold._seqrecord_to_fasta(record)

        assert fasta.startswith(">test_seq")
        assert "ATCGATCGATCG" in fasta

    def test_extract_taxonomy_from_match_basic(self, config):
        """Test extracting taxonomy from basic match structure."""
        bold = BOLD(config)

        match = {
            "family": "Noctuidae",
            "genus": "Spodoptera",
            "species": "frugiperda"
        }

        taxonomy = bold._extract_taxonomy_from_match(match)

        assert taxonomy is not None
        assert taxonomy["family"] == "Noctuidae"
        assert taxonomy["genus"] == "Spodoptera"
        assert taxonomy["species"] == "frugiperda"

    def test_extract_taxonomy_from_match_nested(self, config):
        """Test extracting taxonomy from nested match structure."""
        bold = BOLD(config)

        match = {
            "taxonomy": {
                "Family": "Noctuidae",
                "Genus": "Spodoptera"
            }
        }

        taxonomy = bold._extract_taxonomy_from_match(match)

        assert taxonomy is not None
        assert taxonomy["family"] == "Noctuidae"
        assert taxonomy["genus"] == "Spodoptera"

    def test_extract_taxonomy_from_match_empty(self, config):
        """Test extracting taxonomy from empty match."""
        bold = BOLD(config)

        match = {}
        taxonomy = bold._extract_taxonomy_from_match(match)

        assert taxonomy is None

    def test_get_taxon_at_level_family(self, config):
        """Test getting taxon at family level."""
        bold = BOLD(config)

        taxonomy = {
            "family": "Noctuidae",
            "genus": "Spodoptera",
            "species": "frugiperda"
        }

        taxon = bold._get_taxon_at_level(taxonomy, TaxonomicRank.FAMILY)

        assert taxon is not None
        assert taxon.name == "Noctuidae"
        assert taxon.taxonomic_rank == TaxonomicRank.FAMILY.value

    def test_get_taxon_at_level_genus(self, config):
        """Test getting taxon at genus level."""
        bold = BOLD(config)

        taxonomy = {
            "family": "Noctuidae",
            "genus": "Spodoptera",
            "species": "frugiperda"
        }

        taxon = bold._get_taxon_at_level(taxonomy, TaxonomicRank.GENUS)

        assert taxon is not None
        assert taxon.name == "Spodoptera"
        assert taxon.taxonomic_rank == TaxonomicRank.GENUS.value

    def test_get_taxon_at_level_missing(self, config):
        """Test getting taxon at level not present in taxonomy."""
        bold = BOLD(config)

        taxonomy = {
            "genus": "Spodoptera",
            "species": "frugiperda"
        }

        taxon = bold._get_taxon_at_level(taxonomy, TaxonomicRank.FAMILY)

        assert taxon is None


class TestBOLDDatabases:
    """Test BOLD database configurations."""

    def test_all_databases_have_descriptions(self, config):
        """Test that all databases have non-empty descriptions."""
        bold = BOLD(config)
        databases = bold.get_available_databases()

        for db_code, description in databases.items():
            assert db_code is not None
            assert len(db_code) > 0
            assert description is not None
            assert len(description) > 0

    def test_animal_databases_present(self, config):
        """Test that expected animal databases are present."""
        bold = BOLD(config)
        databases = bold.get_available_databases()

        expected_animal_dbs = [
            "public.tax-derep",
            "species",
            "all.tax-derep",
            "DS-CANREF22",
            "DS-IUCNPUB"
        ]

        for db in expected_animal_dbs:
            assert db in databases

    def test_plant_fungi_databases_present(self, config):
        """Test that plant and fungi databases are present."""
        bold = BOLD(config)
        databases = bold.get_available_databases()

        assert "public.plants" in databases
        assert "public.fungi" in databases
        assert "all.animal-alt" in databases