import pytest
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from unittest.mock import patch, Mock
from barcode_validator.alignment import SequenceHandler


@pytest.fixture
def mock_config():
    class MockConfig:
        def get(self, key):
            if key == 'hmm_file':
                return '/path/to/hmm_file'
            elif key == 'translation_table':
                return 5

    return MockConfig()


@pytest.fixture
def sample_sequence():
    return SeqRecord(Seq("ATGC-N-TGCA"), id="test_seq", name="Test Sequence")


@pytest.fixture
def sample_aligned_sequence():
    return SeqRecord(Seq("ATG-C-N-TG-CA"), id="test_seq", name="Test Aligned Sequence")


def test_align_to_hmm(mock_config, sample_sequence):
    with patch('subprocess.run') as mock_run, \
            patch('tempfile.NamedTemporaryFile') as mock_temp_file, \
            patch('Bio.SeqIO.parse') as mock_parse:
        mock_temp_file().__enter__().name = 'temp_file'
        mock_parse.return_value = iter([sample_sequence])

        result = SequenceHandler.align_to_hmm(sample_sequence, mock_config)

        assert mock_run.called
        assert isinstance(result, SeqRecord)


def test_marker_seqlength(sample_aligned_sequence):
    result = SequenceHandler.marker_seqlength(sample_aligned_sequence)
    assert result == 9  # see line 28, all gaps removed, but Ns count


def test_num_ambiguous(sample_aligned_sequence):
    result = SequenceHandler.num_ambiguous(sample_aligned_sequence)
    assert result == 1  # One 'N' in the sequence


def test_unalign_sequence():
    aligned_seq = SeqRecord(Seq("A-T-G-C-N"), id="test", name="Test")
    result = SequenceHandler.unalign_sequence(aligned_seq)
    assert str(result.seq) == "ATGCN"


def test_translate_sequence():
    dna_seq = SeqRecord(Seq("ATGGAATAA"), id="test", name="Test")
    with patch('Bio.Seq.Seq.translate') as mock_translate:
        mock_translate.return_value = Seq("ME*")
        result = SequenceHandler.translate_sequence(dna_seq, 5)
        assert str(result.seq) == "ME*"


def test_parse_fasta():
    # Get the directory of the current file (assumed to be in the tests folder)
    current_dir = Path(__file__).parent
    example_fasta = current_dir.parent / 'examples' / 'mge.fa'

    result = list(SequenceHandler.parse_fasta(example_fasta))
    assert len(result) == 1
    assert result[0][0] == "BGENL191-23"


def test_parse_fasta_empty_file():
    current_dir = Path(__file__).parent
    example_fasta = current_dir.parent / 'examples' / 'empty_seq.fa'

    result = list(SequenceHandler.parse_fasta(example_fasta))
    assert len(result) == 1
    assert len(result[0][1].seq) == 0


def test_align_fasta_empty_file():
    current_dir = Path(__file__).parent
    example_fasta = current_dir.parent / 'examples' / 'empty_seq.fa'
    result = list(SequenceHandler.parse_fasta(example_fasta))
    seq = result[0][1]
    assert SequenceHandler.align_to_hmm(seq, None) is None


def test_parse_fasta_with_json():
    # Get the directory of the current file (assumed to be in the tests folder)
    current_dir = Path(__file__).parent
    example_fasta = current_dir.parent / 'examples' / 'mge_meta.fa'

    result = list(SequenceHandler.parse_fasta(example_fasta))

    # Check the number of entries
    assert len(result) == 1, "Expected one entry in the FASTA file"

    # Unpack the single result
    process_id, sequence, json_config = result[0]

    # Check the process ID
    assert process_id == "BGENL191-23", "Process ID mismatch"

    # Check the sequence
    assert isinstance(sequence, SeqRecord), "Sequence should be a SeqRecord object"
    assert len(sequence.seq) > 0, "Sequence should not be empty"

    # Check the JSON configuration
    assert json_config is not None, "JSON configuration should be present"
    assert json_config["level"] == "order", "Incorrect 'level' in JSON config"
    assert json_config["translation_table"] == 5, "Incorrect 'translation_table' in JSON config"

    # Check that the JSON part is removed from the sequence description
    assert "{" not in sequence.description, "JSON should be removed from sequence description"


def test_get_stop_codons():
    aa_seq = SeqRecord(Seq("M*EF*GH"), id="test", name="Test")
    result = SequenceHandler.get_stop_codons(aa_seq)
    assert result == [1, 4]


def test_unalign_sequence_with_tilde():
    aligned_seq = SeqRecord(Seq("A~T~G~C~N"), id="test", name="Test")
    result = SequenceHandler.unalign_sequence(aligned_seq)
    assert str(result.seq) == "ATGCN"


def test_unalign_sequence_with_string_input():
    aligned_seq = "A-T-G-C-N"
    result = SequenceHandler.unalign_sequence(aligned_seq)
    assert result == "ATGCN"


def test_unalign_sequence_with_invalid_input():
    with pytest.raises(TypeError):
        SequenceHandler.unalign_sequence(123)


def test_translate_sequence_with_ambiguous_bases(mock_config):
    # When we get the sequence back from alignment by hmmalign, it contains an initial character
    # that shifts the phase of the sequence. This is why the sequence is shifted by one base and
    # why we insert a gap at the beginning of this test sequence.
    dna_seq = SeqRecord(Seq("-ATGNNATAA"), id="test", name="Test")
    result = SequenceHandler.translate_sequence(dna_seq, 5)
    assert str(result.seq) == "M*"  # The NN codon is skipped


def test_marker_seqlength_with_tilde():
    seq = SeqRecord(Seq("ATG~C~N~TGCA"), id="test", name="Test")
    result = SequenceHandler.marker_seqlength(seq)
    assert result == 9  # 12 total - 3 tildes


def test_num_ambiguous_with_multiple_ambiguous_bases():
    seq = SeqRecord(Seq("ATGCNRYSWKMBDHV"), id="test", name="Test")
    result = SequenceHandler.num_ambiguous(seq)
    assert result == 11  # All except A, T, G, C are ambiguous


def test_get_stop_codons_no_stops():
    aa_seq = SeqRecord(Seq("MEFGH"), id="test", name="Test")
    result = SequenceHandler.get_stop_codons(aa_seq)
    assert result == []


def test_marker_seqlength_empty_sequence():
    seq = SeqRecord(Seq(""), id="test", name="Test")
    result = SequenceHandler.marker_seqlength(seq)
    assert result == 0


def test_num_ambiguous_empty_sequence():
    seq = SeqRecord(Seq(""), id="test", name="Test")
    result = SequenceHandler.num_ambiguous(seq)
    assert result == 0


def test_get_stop_codons_empty_sequence():
    aa_seq = SeqRecord(Seq(""), id="test", name="Test")
    result = SequenceHandler.get_stop_codons(aa_seq)
    assert result == []


def test_translate_sequence_empty_sequence(mock_config):
    dna_seq = SeqRecord(Seq(""), id="test", name="Test")
    result = SequenceHandler.translate_sequence(dna_seq, 5)
    assert str(result.seq) == ""
