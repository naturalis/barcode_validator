import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from barcode_validator.config.schema_config import SchemaConfig
from barcode_validator.idservices.bold import BOLD

@pytest.fixture()
def config():
    """Fixture to provide a Config object for tests."""
    config = SchemaConfig()
    return config


class TestBOLDSearchs:
    """Test BOLD search methods."""

    def test_search_valid_fasta(self, config):
        """Test searching with valid FASTA input."""
        bold = BOLD(config)
        record = SeqRecord(
            Seq("AATATTATACTTTATTTTTGCTATATGATCAGGAATAATTGGTTCATCTATAAGATTATTAATTCGAATAGAATTAAGACATCCAGGTATATGAATTAATAATGATCAAATTTATAATTCTTTAGTAACAAGACATGCATTTTTAATAATTTTTTTTATAGTTATACCTTTTATAATTGGTGGATTTGGAAATTATCTAATTCCATTAATATTAGGATCCCCAGATATAGCTTTTCCTCGAATAAATAATATTAGATTTTGACTTCTACCTCCATCATTATTCATATTATTATTAAGAAATATATTTACACCTAATGTAGGTACAGGATGAACTGTATATCCTCCTTTATCTTCTTATTTATTTCATTCATCACCTTCAATTGATATTGCAATCTTTTCTTTACATATATCAGGAATCTCTTCAATTATTGGATCATTAAATTTTATCGTTACTATTTTATTAATAAAAAATTTTTCATTAAATTATGATCAAATTAATTTATTTTCATGATCAGTATGTATTACAGTAATTTTATTAATTCTATCTTTACCAGTATTAGCCGGCGCAATTACTATATTATTATTTGATCGAAATTTTAATACTTCATTTTTTGACCCAATAGGAGGAGGAGATCCAATCCTTTATCAACATTTATTT"),
            id="test_seq"
        )
        results = bold.identify_record(record)

        assert results is not None
        assert isinstance(results, set)
        assert len(results) > 0

if __name__ == '__main__':
    pytest.main()