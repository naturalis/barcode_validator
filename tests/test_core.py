import pytest
from unittest.mock import Mock, patch, mock_open
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Phylo.BaseTree import Tree, Clade
from nbitk.Taxon import Taxon
from nbitk.config import Config
from barcode_validator.dna_analysis_result import DNAAnalysisResult
from barcode_validator.barcode_validator import BarcodeValidator


@pytest.fixture
def mock_config():
    class ConfigMock(Mock):
        def get(self, key, default=None):
            values = {
                'level': 'family',
                'constrain': 'order',
                'hmm_file': 'mock_hmm.hmm',
                'translation_table': 1,
                'ncbi_taxonomy': 'mock_ncbi.tar.gz',
                'bold_sheet_file': 'mock_bold.xlsx',
                'log_level': 'ERROR',
                'tool_name': 'hmmalign',
            }
            return values.get(key, default)

    return ConfigMock(spec=Config)

@pytest.fixture
def barcode_validator(mock_config):
    return BarcodeValidator(mock_config)


@pytest.fixture
def mock_trees(barcode_validator):
    # Create mock NCBI and BOLD trees
    mock_ncbi_tree = Mock(spec=Tree)
    mock_bold_tree = Mock(spec=Tree)

    # Set up the BOLD tree structure
    mock_bold_root = Mock(spec=Taxon)
    mock_bold_tree.root = mock_bold_root
    mock_bold_tree.get_terminals = Mock()

    barcode_validator.ncbi_tree = mock_ncbi_tree
    barcode_validator.bold_tree = mock_bold_tree

    return barcode_validator


def test_initialize(barcode_validator):
    with patch('barcode_validator.barcode_validator.NCBIParser') as mock_ncbi_parser, \
            patch('barcode_validator.barcode_validator.BOLDParser') as mock_bold_parser, \
            patch('barcode_validator.barcode_validator.tarfile.open'), \
            patch('builtins.open', mock_open(read_data=b'mock_data')):
        mock_ncbi_parser.return_value.parse.return_value = Mock(spec=Tree)
        mock_bold_parser.return_value.parse.return_value = Mock(spec=Tree)

        barcode_validator.initialize()

        assert isinstance(barcode_validator.ncbi_tree, Mock)
        assert isinstance(barcode_validator.bold_tree, Mock)


@patch('barcode_validator.sequence_handler.SequenceHandler.parse_fasta')
def test_validate_fasta(mock_parse_fasta, barcode_validator, mock_config):
    mock_parse_fasta.return_value = [
        (SeqRecord(Seq('ATCG'), id='seq1'), {}),
        (SeqRecord(Seq('GCTA'), id='seq2'), {})
    ]

    with patch.object(barcode_validator, 'validate_record', return_value=Mock(spec=DNAAnalysisResult)) as mock_validate:
        results = barcode_validator.validate_fasta('mock.fasta', mock_config)

        assert len(results) == 2
        assert mock_validate.call_count == 2


def test_validate_record(mock_trees, mock_config):
    record = SeqRecord(Seq('ATCG'), id='seq1')
    record.annotations['bcdm_fields'] = { 'processid': 'process1' }

    with patch.object(mock_trees, 'validate_sequence_quality') as mock_validate_quality, \
            patch.object(mock_trees, 'validate_taxonomy') as mock_validate_taxonomy:
        result = DNAAnalysisResult('process1')
        mock_trees.validate_record(record, mock_config, result)

        assert isinstance(result, DNAAnalysisResult)
        mock_validate_quality.assert_called_once()
        mock_validate_taxonomy.assert_called_once()


@pytest.mark.xfail(reason="This test is just mock object performance art. Marked for deletion.")
def test_validate_taxonomy(mock_trees, mock_config):
    record = SeqRecord(Seq('ATCG'), id='seq1')
    record.annotations['bcdm_fields'] = { 'processid': 'process1' }
    result = DNAAnalysisResult('process1')

    mock_species = Mock(spec=Taxon)
    mock_species.name = "TestSpecies"
    mock_exp_taxon = Mock(spec=Taxon)
    mock_exp_taxon.name = "TestFamily"
    mock_exp_taxon.taxonomic_rank = "family"
    mock_obs_taxon = [Mock(spec=Taxon)]
    mock_obs_taxon[0].name = "ObservedFamily"

    mock_trees.bold_tree.root.get_path.return_value = [mock_exp_taxon, mock_species]

    with patch.object(mock_trees, 'get_node_by_processid', return_value=mock_species), \
            patch.object(mock_trees, 'build_constraint', return_value='1234'), \
            patch('barcode_validator.blast_runner.BlastRunner') as MockBlastRunner:
        MockBlastRunner.return_value.run_localblast.return_value = mock_obs_taxon

        mock_trees.validate_taxonomy(mock_config, record, result)

        # These tests do nothing other than getting the mocks back that were assigned.
        # Might as well remove these because nothing functional is really being done.
        assert result.species == mock_species
        assert result.exp_taxon == mock_exp_taxon
        assert result.obs_taxon == mock_obs_taxon
        assert result.error is None

    # Test error case when species is not found
    with patch.object(mock_trees, 'get_node_by_processid', return_value=None):
        mock_trees.validate_taxonomy(mock_config, record, result)
        assert result.error == "process1 not in BOLD"

    # Test error case when BLAST fails
    with patch.object(mock_trees, 'get_node_by_processid', return_value=mock_species), \
            patch.object(mock_trees, 'build_constraint', return_value='1234'), \
            patch('barcode_validator.blast_runner.BlastRunner') as MockBlastRunner:
        MockBlastRunner.return_value.run_localblast.return_value = None

        mock_trees.validate_taxonomy(mock_config, record, result)
        assert result.error == "Local BLAST failed for sequence 'ATCG'"


def test_build_constraint(mock_trees):
    mock_bold_tip = Mock(spec=Taxon)
    mock_bold_tip.name = 'TestSpecies'
    mock_bold_tip.taxonomic_rank = 'species'

    mock_bold_anc = Mock(spec=Taxon)
    mock_bold_anc.name = 'TestOrder'
    mock_bold_anc.taxonomic_rank = 'order'

    mock_ncbi_anc = Mock(spec=Taxon)
    mock_ncbi_anc.name = 'TestOrder'
    mock_ncbi_anc.taxonomic_rank = 'order'
    mock_ncbi_anc.guids = {'taxon': '1234'}

    mock_trees.bold_tree.root.get_path.return_value = [mock_bold_anc, mock_bold_tip]
    mock_trees.ncbi_tree.get_nonterminals.return_value = [mock_ncbi_anc]

    result = mock_trees.build_constraint(mock_bold_tip, 'order')

    assert result == '1234'


@patch('barcode_validator.sequence_handler.SequenceHandler.align_to_hmm')
@patch('barcode_validator.sequence_handler.SequenceHandler.translate_sequence')
@patch('barcode_validator.sequence_handler.SequenceHandler.get_stop_codons')
@patch('barcode_validator.sequence_handler.SequenceHandler.marker_seqlength')
@patch('barcode_validator.sequence_handler.SequenceHandler.num_ambiguous')
def test_validate_sequence_quality(mock_num_ambiguous, mock_marker_seqlength, mock_get_stop_codons,
                                   mock_translate_sequence, mock_align_to_hmm, mock_config):
    record = SeqRecord(Seq('ATCG'), id='seq1')
    result = DNAAnalysisResult('process1')

    mock_align_to_hmm.return_value = 'ATCG----'
    mock_translate_sequence.return_value = 'M'
    mock_get_stop_codons.return_value = []
    mock_marker_seqlength.return_value = 4
    mock_num_ambiguous.side_effect = [0, 0]

    bv = BarcodeValidator(mock_config)
    bv.validate_sequence_quality(mock_config, record, result)

    assert result.full_length == 4
    assert result.full_ambiguities == 0
    assert result.stop_codons == []
    assert result.seq_length == 4
    assert result.ambiguities == 0


if __name__ == '__main__':
    pytest.main()