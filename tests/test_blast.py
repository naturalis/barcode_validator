import os
from pathlib import Path
import pytest
from Bio.SeqIO import parse
from nbitk.config import Config
from nbitk.Tools import Blastn
from barcode_validator.idservices.ncbi import NCBI
from barcode_validator.resolvers.factory import ResolverFactory
from barcode_validator.resolvers.taxonomy import TaxonomicRank, TaxonResolver, TaxonomicBackbone

# Test data path handling
TEST_DATA_DIR = Path(__file__).parent / "data"
BOLD_SAMPLE = TEST_DATA_DIR / "BGE00146_MGE-BGE_r1_1.3_1.5_s50_100.fasta"
BOLD_SHEET = TEST_DATA_DIR / "bold.xlsx"  # We'll need to create this
HMM_DIR = TEST_DATA_DIR / "hmm_profiles"  # We'll need to create this
NCBI_TAXDUMP = TEST_DATA_DIR / "taxdump.tar.gz"
GENBANK = '/home/rutger.vos/data/ncbi/nt/nt'
BLASTDB = '/home/rutger.vos/data/ncbi/nt'
BLAST_REPORT = TEST_DATA_DIR / "tmp_o8h137e.fasta.tsv"

@pytest.fixture
def config():
    """Test fixture that provides a Config object configured for COI-5P validation"""
    conf = Config()
    conf.config_data = {
        'log_level': 'DEBUG',
        'validate_taxonomy': True,  # Test BLASTing just yet
        'validate_structure': False,  # Disable structural validation
        'marker': 'COI-5P',  # Using COI-5P marker
        'ncbi_file': str(NCBI_TAXDUMP),
        'db': GENBANK,
        'num_threads': 56,
        'evalue': 1e-5,
        'max_target_seqs': 10,
        'word_size': 28,
        'constraint_rank': 'class',
        'BLASTDB': BLASTDB,
        'BLASTDB_LMDB_MAP_SIZE': '1000G',
        'bold_file': str(BOLD_SHEET),
        'exp_taxonomy_type': 'bold',
        'reference_taxonomy': 'ncbi'
    }
    conf.initialized = True
    os.environ['BLASTDB'] = conf.get('BLASTDB')
    os.environ['BLASTDB_LMDB_MAP_SIZE'] = conf.get('BLASTDB_LMDB_MAP_SIZE')    
    return conf

@pytest.fixture
def blastn(config):
    """Fixture providing Blastn instance"""
    blastn = Blastn(config)
    blastn.set_db(config.get('db'))
    blastn.set_num_threads(config.get('num_threads'))
    blastn.set_evalue(config.get('evalue'))
    blastn.set_max_target_seqs(config.get('max_target_seqs'))
    blastn.set_word_size(config.get('word_size'))
    blastn.set_task('megablast')
    blastn.set_outfmt("6 qseqid sseqid pident length qstart qend sstart send evalue bitscore staxids")
    return blastn

@pytest.fixture
def taxonomy_resolver(config):
    """Provides initialized TaxonomyResolver"""
    resolver = ResolverFactory.create_resolver(config, TaxonomicBackbone.NCBI)
    resolver.load_tree(NCBI_TAXDUMP)
    return resolver

@pytest.fixture
def ncbi(config, blastn, taxonomy_resolver):
    """Fixture providing NCBI instance"""
    ncbi = NCBI(config)
    ncbi.set_blastn(blastn)
    ncbi.set_taxonomy_resolver(taxonomy_resolver)
    return ncbi

@pytest.fixture
def test_data():
    """Fixture providing test data"""
    records = []
    with open(BOLD_SAMPLE, "r") as f:
        for record in parse(f, "fasta"):
            records.append(record)
    return records

def test_blastn_initialization(blastn, config):
    """Test initialization of Blastn instance"""
    # This test is very stupid because here we're just verifying that the configurable values indeed 
    # are passed into the Blastn instance in the fixture.
    assert blastn is not None
    assert config.get('db') == blastn.get_parameter('db')
    assert config.get('num_threads') == int(blastn.get_parameter('num_threads'))
    assert config.get('evalue') == float(blastn.get_parameter('evalue'))
    assert config.get('max_target_seqs') == int(blastn.get_parameter('max_target_seqs'))
    assert config.get('word_size') == int(blastn.get_parameter('word_size'))
    assert blastn.get_parameter('task') == 'megablast'
    assert blastn.get_parameter('outfmt') == "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore staxids"

@pytest.mark.skipif(not Path(NCBI_TAXDUMP).exists(), reason="NCBI taxonomy does not exist")
def test_ncbi_initialization(ncbi, config):
    """Test initialization of NCBI instance"""
    assert ncbi is not None
    assert config.get('db') == ncbi.blastn.get_parameter('db')

# Skip test if GENBANK does not exist
@pytest.mark.skipif(not Path(BLASTDB).exists(), reason="GENBANK does not exist")
def test_run_localblast(ncbi, test_data):
    """Test identify_record method of NCBI instance"""
    assert ncbi is not None
    blast_report = ncbi.run_localblast(test_data[0], 51653)
    assert blast_report is not None
    assert Path(blast_report).exists()
    assert Path(f"{blast_report}.tsv").exists()

    # Clean up temporary files
    os.unlink(blast_report)
    os.unlink(f"{blast_report}.tsv")

@pytest.mark.skipif(not Path(NCBI_TAXDUMP).exists(), reason="NCBI taxonomy does not exist")
def test_parse_blast_result(ncbi):
    """Test parse_blast_result method of NCBI instance"""
    assert ncbi is not None
    taxids = ncbi.parse_blast_result(BLAST_REPORT.as_posix())
    assert taxids is not None
    assert len(taxids) > 0

@pytest.mark.skipif(not Path(NCBI_TAXDUMP).exists(), reason="NCBI taxonomy does not exist")
def test_collect_higher_taxa(ncbi):
    """Test collect_higher_taxa method of NCBI instance"""
    assert ncbi is not None
    taxids = ncbi.parse_blast_result(BLAST_REPORT.as_posix())
    taxa = ncbi.collect_higher_taxa(taxids, TaxonomicRank.FAMILY)
    assert taxa is not None
    assert len(taxa) > 0

# Skip test if GENBANK does not exist
@pytest.mark.skipif(not Path(BLASTDB).exists(), reason="GENBANK does not exist")
def test_identify_record(ncbi, test_data):
    """Test identify_record method of NCBI instance"""
    assert ncbi is not None
    record = test_data[0]
    taxa = ncbi.identify_record(record, TaxonomicRank.FAMILY)
    assert taxa is not None
    assert len(taxa) > 0

if __name__ == "__main__":
    pytest.main()