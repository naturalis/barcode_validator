import pytest
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from unittest.mock import patch, Mock, MagicMock
from nbitk.config import Config
from barcode_validator.sequence_handler import SequenceHandler

@pytest.fixture
def mock_config():
    config = Mock(spec=Config)
    config.get.side_effect = lambda key, default=None: {
        'level': 'family',
        'constrain': 'order',
        'hmm_file': 'mock_hmm.hmm',
        'translation_table': 1,
        'ncbi_taxonomy': 'mock_ncbi.tar.gz',
        'bold_sheet_file': 'mock_bold.xlsx',
        'log_level': 'ERROR',
        'tool_name': 'hmmalign',
    }.get(key, default)
    return config

@pytest.fixture
def sequence_handler(mock_config):
    return SequenceHandler(mock_config)

@pytest.fixture
def sample_sequence():
    return SeqRecord(Seq("ATGC-N-TGCA"), id="test_seq", name="Test Sequence")

@pytest.fixture
def sample_aligned_sequence():
    return SeqRecord(Seq("ATG-C-N-TG-CA"), id="test_seq", name="Test Aligned Sequence")

def test_align_to_hmm(sequence_handler, sample_sequence):
    with patch('tempfile.NamedTemporaryFile') as mock_temp_file, \
         patch('Bio.SeqIO.parse') as mock_parse:
        mock_temp_file().__enter__().name = 'temp_file'
        mock_parse.return_value = iter([sample_sequence])
        sequence_handler.hmmalign.run = MagicMock(return_value=0)
        result = sequence_handler.align_to_hmm(sample_sequence)

        assert sequence_handler.hmmalign.run.called
        assert isinstance(result, SeqRecord)

def test_marker_seqlength(sequence_handler, sample_aligned_sequence):
    result = sequence_handler.marker_seqlength(sample_aligned_sequence)
    assert result == 9  # see line 28, all gaps removed, but Ns count

def test_num_ambiguous(sequence_handler, sample_aligned_sequence):
    result = sequence_handler.num_ambiguous(sample_aligned_sequence)
    assert result == 1  # One 'N' in the sequence

def test_unalign_sequence(sequence_handler):
    aligned_seq = SeqRecord(Seq("A-T-G-C-N"), id="test", name="Test")
    result = sequence_handler.unalign_sequence(aligned_seq)
    assert str(result.seq) == "ATGCN"

def test_translate_sequence_standard_code(sequence_handler):
    dna_seq = SeqRecord(Seq("ATGGAATAA"), id="test", name="Test")
    result = sequence_handler.translate_sequence(dna_seq, 1)  # Standard genetic code
    assert str(result.seq) == "WN"

def test_translate_sequence_invertebrate_mitochondrial_code(sequence_handler):
    dna_seq = SeqRecord(Seq("ATGGAATAA"), id="test", name="Test")
    result = sequence_handler.translate_sequence(dna_seq, 5)  # Invertebrate Mitochondrial Code
    assert str(result.seq) == "WN"

def test_parse_fasta(sequence_handler):
    current_dir = Path(__file__).parent
    example_fasta = current_dir.parent / 'examples' / 'mge.fa'

    result = list(sequence_handler.parse_fasta(example_fasta))
    assert len(result) == 1
    assert result[0][0].annotations['bcdm_fields']['processid'] == "BGENL191-23"

def test_parse_fasta_empty_file(sequence_handler):
    current_dir = Path(__file__).parent
    example_fasta = current_dir.parent / 'examples' / 'empty_seq.fa'

    result = list(sequence_handler.parse_fasta(example_fasta))
    assert len(result) == 1
    assert len(result[0][0].seq) == 0

def test_align_fasta_empty_file(sequence_handler):
    current_dir = Path(__file__).parent
    example_fasta = current_dir.parent / 'examples' / 'empty_seq.fa'
    result = list(sequence_handler.parse_fasta(example_fasta))
    seq = result[0][0]
    assert sequence_handler.align_to_hmm(seq) is None

def test_parse_fasta_with_json(sequence_handler):
    current_dir = Path(__file__).parent
    example_fasta = current_dir.parent / 'examples' / 'mge_meta.fa'

    result = list(sequence_handler.parse_fasta(example_fasta))

    assert len(result) == 1, "Expected one entry in the FASTA file"

    sequence, json_config = result[0]

    assert sequence.annotations['bcdm_fields']['processid'] == "BGENL191-23", "Process ID mismatch"
    assert isinstance(sequence, SeqRecord), "Sequence should be a SeqRecord object"
    assert len(sequence.seq) > 0, "Sequence should not be empty"
    assert json_config is not None, "JSON configuration should be present"
    assert json_config["level"] == "order", "Incorrect 'level' in JSON config"
    assert json_config["translation_table"] == 5, "Incorrect 'translation_table' in JSON config"
    assert "{" not in sequence.description, "JSON should be removed from sequence description"

def test_get_stop_codons(sequence_handler):
    aa_seq = SeqRecord(Seq("M*EF*GH"), id="test", name="Test")
    result = sequence_handler.get_stop_codons(aa_seq)
    assert result == [1, 4]

def test_unalign_sequence_with_tilde(sequence_handler):
    aligned_seq = SeqRecord(Seq("A~T~G~C~N"), id="test", name="Test")
    result = sequence_handler.unalign_sequence(aligned_seq)
    assert str(result.seq) == "ATGCN"

def test_unalign_sequence_with_string_input(sequence_handler):
    aligned_seq = "A-T-G-C-N"
    result = sequence_handler.unalign_sequence(aligned_seq)
    assert result == "ATGCN"

def test_unalign_sequence_with_invalid_input(sequence_handler):
    with pytest.raises(TypeError):
        sequence_handler.unalign_sequence(123)

def test_translate_sequence_with_ambiguous_bases(sequence_handler):
    dna_seq = SeqRecord(Seq("-ATGNNATAA"), id="test", name="Test")
    result = sequence_handler.translate_sequence(dna_seq, 5)
    assert str(result.seq) == "M*"  # The NN codon is skipped

def test_marker_seqlength_with_tilde(sequence_handler):
    seq = SeqRecord(Seq("ATG~C~N~TGCA"), id="test", name="Test")
    result = sequence_handler.marker_seqlength(seq)
    assert result == 9  # 12 total - 3 tildes

def test_num_ambiguous_with_multiple_ambiguous_bases(sequence_handler):
    seq = SeqRecord(Seq("ATGCNRYSWKMBDHV"), id="test", name="Test")
    result = sequence_handler.num_ambiguous(seq)
    assert result == 11  # All except A, T, G, C are ambiguous

def test_get_stop_codons_no_stops(sequence_handler):
    aa_seq = SeqRecord(Seq("MEFGH"), id="test", name="Test")
    result = sequence_handler.get_stop_codons(aa_seq)
    assert result == []

def test_marker_seqlength_empty_sequence(sequence_handler):
    seq = SeqRecord(Seq(""), id="test", name="Test")
    result = sequence_handler.marker_seqlength(seq)
    assert result == 0

def test_num_ambiguous_empty_sequence(sequence_handler):
    seq = SeqRecord(Seq(""), id="test", name="Test")
    result = sequence_handler.num_ambiguous(seq)
    assert result == 0

def test_get_stop_codons_empty_sequence(sequence_handler):
    aa_seq = SeqRecord(Seq(""), id="test", name="Test")
    result = sequence_handler.get_stop_codons(aa_seq)
    assert result == []

def test_translate_sequence_empty_sequence(sequence_handler):
    dna_seq = SeqRecord(Seq(""), id="test", name="Test")
    result = sequence_handler.translate_sequence(dna_seq, 5)
    assert str(result.seq) == ""

if __name__ == "__main__":
    pytest.main()