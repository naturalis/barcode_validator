import pytest
from unittest.mock import Mock, patch, mock_open
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Phylo.BaseTree import Tree
from nbitk.config import Config
from barcode_validator.taxonomy import BlastRunner
from nbitk.Taxon import Taxon


@pytest.fixture
def mock_config():
    config = Mock(spec=Config)
    config.get.side_effect = lambda key: {
        'ncbi_tree': Mock(spec=Tree),
        'blast_db': 'mock_db',
        'num_threads': 4,
        'evalue': 0.001,
        'max_target_seqs': 10,
        'word_size': 11,
        'BLASTDB_LMDB_MAP_SIZE': 1000000000,
        'BLASTDB': '/path/to/blastdb',
        'log_level': 'ERROR',
    }[key]
    return config


@pytest.fixture
def blast_runner(mock_config):
    return BlastRunner(mock_config)


def test_init(blast_runner, mock_config):
    assert blast_runner.blast_db == 'mock_db'
    assert blast_runner.num_threads == 4
    assert blast_runner.evalue == 0.001
    assert blast_runner.max_target_seqs == 10
    assert blast_runner.word_size == 11
    assert blast_runner.BLASTDB_LMDB_MAP_SIZE == 1000000000
    assert blast_runner.BLASTDB == '/path/to/blastdb'


@patch('barcode_validator.taxonomy.tempfile.NamedTemporaryFile')
@patch('barcode_validator.taxonomy.SeqIO.write')
@patch('barcode_validator.taxonomy.subprocess.Popen')
def test_run_localblast(mock_popen, mock_seqio_write, mock_temp_file, blast_runner):
    mock_temp_file.return_value.__enter__.return_value.name = 'temp_file.fasta'
    mock_process = Mock()
    mock_process.stdout = ['BLAST output']
    mock_process.stderr = []
    mock_process.wait.return_value = 0
    mock_popen.return_value = mock_process

    sequence = SeqRecord(Seq('ATCG'), id='test_seq')
    constraint = '1234'

    with patch.object(blast_runner, 'parse_blast_result', return_value=['Family1', 'Family2']):
        result = blast_runner.run_localblast(sequence, constraint)

    assert result == ['Family1', 'Family2']
    mock_seqio_write.assert_called_once()
    mock_popen.assert_called_once()


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