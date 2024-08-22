import pytest
from unittest.mock import patch, MagicMock
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from barcode_validator.core import BarcodeValidator
from barcode_validator.result import DNAAnalysisResult
from nbitk.Taxon import Taxon


@pytest.fixture
def mock_config():
    with patch('barcode_validator.core.Config') as MockConfig:
        config = MockConfig.return_value
        config.get.side_effect = lambda x: {
            'ncbi_taxonomy': 'path/to/ncbi_taxonomy',
            'bold_sheet_file': 'path/to/bold_sheet',
            'level': 'family'
        }.get(x)
        yield config


@pytest.fixture
def barcode_validator(mock_config):
    return BarcodeValidator(mock_config)


@pytest.fixture
def mock_trees():
    with patch('barcode_validator.core.read_ncbi_taxonomy') as mock_ncbi, \
            patch('barcode_validator.core.read_bold_taxonomy') as mock_bold:
        mock_ncbi.return_value = MagicMock()
        mock_bold.return_value = MagicMock()
        yield mock_ncbi.return_value, mock_bold.return_value


def test_init(barcode_validator, mock_config):
    assert barcode_validator.config == mock_config
    assert barcode_validator.ncbi_tree is None
    assert barcode_validator.bold_tree is None


def test_initialize(barcode_validator, mock_trees):
    barcode_validator.initialize()
    assert barcode_validator.ncbi_tree == mock_trees[0]
    assert barcode_validator.bold_tree == mock_trees[1]


@patch('barcode_validator.core.parse_fasta')
def test_validate_fasta(mock_parse_fasta, barcode_validator):
    mock_parse_fasta.return_value = [
        ('process_id1', SeqRecord(Seq('ATCG'), id='seq1')),
        ('process_id2', SeqRecord(Seq('GCTA'), id='seq2'))
    ]

    with patch.object(barcode_validator, 'validate_record') as mock_validate_record:
        mock_validate_record.side_effect = [
            DNAAnalysisResult('process_id1'),
            DNAAnalysisResult('process_id2')
        ]

        results = barcode_validator.validate_fasta('dummy.fasta')

        assert len(results) == 2
        assert all(isinstance(result, DNAAnalysisResult) for result in results)
        assert mock_validate_record.call_count == 2


@patch('barcode_validator.core.get_tip_by_processid')
@patch('barcode_validator.core.run_localblast')
@patch('barcode_validator.core.align_to_hmm')
@patch('barcode_validator.core.translate_sequence')
@patch('barcode_validator.core.get_stop_codons')
@patch('barcode_validator.core.marker_seqlength')
@patch('barcode_validator.core.num_ambiguous')
def test_validate_record(mock_num_ambiguous, mock_marker_seqlength, mock_get_stop_codons,
                         mock_translate_sequence, mock_align_to_hmm, mock_run_localblast,
                         mock_get_tip_by_processid, barcode_validator, mock_trees):
    barcode_validator.ncbi_tree, barcode_validator.bold_tree = mock_trees

    mock_species = Taxon(name="Mock Species", taxonomic_rank="species")
    mock_get_tip_by_processid.return_value = mock_species
    mock_run_localblast.return_value = [Taxon(name="FamilyA", taxonomic_rank="family"),
                                        Taxon(name="FamilyB", taxonomic_rank="family")]
    mock_align_to_hmm.return_value = SeqRecord(Seq('ATCG' * 125), id='aligned_seq1')
    mock_translate_sequence.return_value = SeqRecord(Seq('M' * 166), id='translated_seq1')
    mock_get_stop_codons.return_value = []
    mock_marker_seqlength.return_value = 500
    mock_num_ambiguous.side_effect = [0, 2]  # For full sequence and aligned sequence

    # Mock the get_path method of the bold_tree
    mock_family = Taxon(name="MockFamily", taxonomic_rank="family")
    barcode_validator.bold_tree.root.get_path.return_value = [
        Taxon(name="MockKingdom", taxonomic_rank="kingdom"),
        Taxon(name="MockPhylum", taxonomic_rank="phylum"),
        Taxon(name="MockClass", taxonomic_rank="class"),
        Taxon(name="MockOrder", taxonomic_rank="order"),
        mock_family,
        Taxon(name="MockGenus", taxonomic_rank="genus"),
        mock_species
    ]

    record = SeqRecord(Seq('ATCG' * 100), id='seq1')
    result = barcode_validator.validate_record('process_id', record)

    assert isinstance(result, DNAAnalysisResult)
    assert result.process_id == 'process_id'
    assert result.full_length == 400
    assert result.full_ambiguities == 0
    assert result.species == mock_species
    assert result.exp_taxon == mock_family
    assert [taxon.name for taxon in result.obs_taxon] == ['FamilyA', 'FamilyB']
    assert result.stop_codons == []
    assert result.seq_length == 500
    assert result.ambiguities == 2
