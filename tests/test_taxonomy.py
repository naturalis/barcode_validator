import pytest
from unittest.mock import Mock, patch, mock_open
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Phylo.BaseTree import Tree
from nbitk.config import Config
from nbitk.Tools import Blastn
from barcode_validator.blast_runner import BlastRunner
from nbitk.Taxon import Taxon

@pytest.fixture
def mock_config():
    config = Mock(spec=Config)
    config.get.side_effect = lambda key, default=None: {
        'blast_db': 'mock_db',
        'num_threads': 4,
        'evalue': 0.001,
        'max_target_seqs': 10,
        'word_size': 11,
        'BLASTDB_LMDB_MAP_SIZE': 1000000000,
        'BLASTDB': '/path/to/blastdb',
        'log_level': 'ERROR',
        'tool_name': 'blastn',
    }.get(key, default)
    return config

@pytest.fixture
def blast_runner(mock_config):
    with patch('barcode_validator.blast_runner.Blastn') as MockBlastn:
        MockBlastn.return_value = Mock(spec=Blastn)
        return BlastRunner(mock_config)

def test_init(blast_runner, mock_config):
    assert isinstance(blast_runner.blastn, Mock)
    blast_runner.blastn.set_db.assert_called_once_with('mock_db')
    blast_runner.blastn.set_num_threads.assert_called_once_with(4)
    blast_runner.blastn.set_evalue.assert_called_once_with(0.001)
    blast_runner.blastn.set_max_target_seqs.assert_called_once_with(10)
    blast_runner.blastn.set_word_size.assert_called_once_with(11)
    blast_runner.blastn.set_task.assert_called_once_with('megablast')
    blast_runner.blastn.set_outfmt.assert_called_once_with("6 qseqid sseqid pident length qstart qend sstart send evalue bitscore staxids")
#    assert blast_runner.BLASTDB_LMDB_MAP_SIZE == 1000000000
    assert blast_runner.BLASTDB == '/path/to/blastdb'

@patch('barcode_validator.blast_runner.tempfile.NamedTemporaryFile')
@patch('barcode_validator.blast_runner.SeqIO.write')
@patch('os.environ')
def test_run_localblast(mock_environ, mock_seqio_write, mock_temp_file, blast_runner):
    mock_temp_file.return_value.__enter__.return_value.name = 'temp_file.fasta'
    blast_runner.blastn.run.return_value = 0

    sequence = SeqRecord(Seq('ATCG'), id='test_seq')
    constraint = '1234'

    with patch.object(blast_runner, 'parse_blast_result', return_value=['Family1', 'Family2']):
        result = blast_runner.run_localblast(sequence, constraint)

    assert result == ['Family1', 'Family2']
    mock_seqio_write.assert_called_once()
    blast_runner.blastn.set_query.assert_called_once_with('temp_file.fasta')
    blast_runner.blastn.set_taxids.assert_called_once_with(['1234'])
    blast_runner.blastn.set_out.assert_called_once_with('temp_file.fasta.tsv')
    blast_runner.blastn.run.assert_called_once()
    #mock_environ.__setitem__.assert_any_call('BLASTDB_LMDB_MAP_SIZE', '1000G')
    mock_environ.__setitem__.assert_any_call('BLASTDB', '/path/to/blastdb')

def test_run_localblast_empty_sequence(blast_runner):
    sequence = SeqRecord(Seq(''), id='empty_seq')
    constraint = '1234'

    result = blast_runner.run_localblast(sequence, constraint)

    assert result is None

@patch('builtins.open', new_callable=mock_open, read_data='seq1\tseq2\t100\t100\t1\t100\t1\t100\t0.0\t200\t9606\n')
def test_parse_blast_result(mock_file, blast_runner):
    with patch.object(blast_runner, 'collect_higher_taxa', return_value=['Family1']):
        result = blast_runner.parse_blast_result('mock_blast_result.tsv', 'family')

    assert result == ['Family1']
    mock_file.assert_called_once_with('mock_blast_result.tsv', 'r')

def test_collect_higher_taxa(blast_runner):
    # Create a mock tree structure
    root = Taxon(name='root')
    family1 = Taxon(name='Family1', taxonomic_rank='family')
    family2 = Taxon(name='Family2', taxonomic_rank='family')
    species1 = Taxon(name='Species1', taxonomic_rank='species')
    species2 = Taxon(name='Species2', taxonomic_rank='species')

    root.clades = [family1, family2]
    family1.clades = [species1]
    family2.clades = [species2]

    species1.guids = {'taxon': '1'}
    species2.guids = {'taxon': '2'}

    blast_runner.ncbi_tree = Tree(root)

    taxids = {'1', '2'}
    result = blast_runner.collect_higher_taxa(taxids, 'family')

    assert len(result) == 2
    assert all(node.taxonomic_rank == 'family' for node in result)
    assert set(node.name for node in result) == {'Family1', 'Family2'}

if __name__ == '__main__':
    pytest.main()