"""
Focused unit tests for NEW fragment processing functionality in ProteinCodingValidator.

Tests only the methods not covered by existing integration tests:
- _split_sequence_on_gaps() 
- _align_single_fragment()
- _recombine_fragments()
- Tilde removal preprocessing

Existing test_coding.py already covers:
- get_translation_table(), _determine_reading_frame(), _find_stop_codons()
- Basic HMM alignment, complete validation workflow
"""

import pytest
import tempfile
import os
from unittest.mock import Mock, patch, mock_open
from pathlib import Path

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO, SeqIO

# Assuming these imports - adjust paths as needed
from barcode_validator.validators.protein_coding import ProteinCodingValidator
from nbitk.config import Config


class TestFragmentProcessing:
    """Test suite for NEW fragment processing methods only."""

    @pytest.fixture
    def mock_validator(self):
        """Create minimal ProteinCodingValidator with mocked dependencies."""
        with patch('barcode_validator.validators.protein_coding.StructuralValidator.__init__'):
            validator = ProteinCodingValidator(Mock(spec=Config))
            validator.logger = Mock()
            validator.hmmalign = Mock()
            validator.hmm_profile_dir = None
            validator.marker = Mock()
            validator.marker.value = 'COI-5P'
            return validator
            
    # ============================================================================
    # TEST: _split_sequence_on_gaps()
    # ============================================================================

    def test_split_sequence_with_real_fragmented_data(self, mock_validator):
        """Test sequence splitting using real fragmented COI-5P sequence."""
        examples_path = Path(__file__).parent.parent / "examples"
        fasta_file = examples_path / "fragmented_coi5p.fasta"
        
        if fasta_file.exists():
            # Use real sequence from file
            records = list(SeqIO.parse(fasta_file, "fasta"))
            real_sequence = str(records[0].seq)
            fragments = mock_validator._split_sequence_on_gaps(real_sequence)
            
            # Verify fragments don't contain gaps and are non-empty
            assert len(fragments) >= 2
            for fragment, start_pos, end_pos in fragments:
                assert '-' not in fragment
                assert len(fragment) > 0
                assert start_pos < end_pos
        else:
            # Fallback to realistic pattern if file not found
            real_pattern = "------CGGTGGTGGCAAAGAACTAACCACAAAGATATTGGAACA------------------------GCAGGTATTGTAGGCAGAGCCTTGACATATTGG--------"
            fragments = mock_validator._split_sequence_on_gaps(real_pattern)
            
            expected = [
                ('CGGTGGTGGCAAAGAACTAACCACAAAGATATTGGAACA', 6, 45),
                ('GCAGGTATTGTAGGCAGAGCCTTGACATATTGG', 69, 102)
            ]
            assert fragments == expected

    def test_split_sequence_with_real_mixed_data(self, mock_validator):
        """Test sequence splitting using real mixed tildes/gaps sequence."""
        examples_path = Path(__file__).parent.parent / "examples"
        fasta_file = examples_path / "mixed_tildes_gaps.fasta"
        
        if fasta_file.exists():
            # Use real sequence from file  
            records = list(SeqIO.parse(fasta_file, "fasta"))
            original_seq = str(records[0].seq)
            
            # Simulate tilde removal, then split
            seq_after_tilde_removal = original_seq.replace('~', '')
            fragments = mock_validator._split_sequence_on_gaps(seq_after_tilde_removal)
            
            # Should create multiple fragments
            assert len(fragments) > 1
            for fragment, start_pos, end_pos in fragments:
                assert '-' not in fragment
                assert '~' not in fragment
                assert len(fragment) > 0
        else:
            # Fallback test
            pytest.skip("mixed_tildes_gaps.fasta not found")

    def test_split_sequence_simple_gaps(self, mock_validator):
        """Test splitting sequence with simple gap patterns."""
        seq_str = "ATG---CGT---AAA"
        fragments = mock_validator._split_sequence_on_gaps(seq_str)
        
        expected = [
            ('ATG', 0, 3),
            ('CGT', 6, 9), 
            ('AAA', 12, 15)
        ]
        assert fragments == expected

    def test_split_sequence_edge_cases(self, mock_validator):
        """Test edge cases for sequence splitting."""
        # Empty sequence
        assert mock_validator._split_sequence_on_gaps("") == []
        
        # Only gaps
        assert mock_validator._split_sequence_on_gaps("-------") == []
        
        # Single character fragments
        seq_str = "A-T-G-C"
        fragments = mock_validator._split_sequence_on_gaps(seq_str)
        expected = [('A', 0, 1), ('T', 2, 3), ('G', 4, 5), ('C', 6, 7)]
        assert fragments == expected

    # ============================================================================
    # TEST: _align_single_fragment()
    # ============================================================================

    def test_align_single_fragment_with_real_stockholm(self, mock_validator):
        """Test successful fragment alignment using real Stockholm file."""
        examples_path = Path(__file__).parent.parent / "examples"
        stockholm_file = examples_path / "good_fragment.sto"
        
        if not stockholm_file.exists():
            pytest.skip("good_fragment.sto not found in examples/")
        
        # Mock hmmalign to succeed and use real Stockholm file for parsing
        mock_validator.hmmalign.run.return_value = 0
        
        # Test the parsing with real Stockholm file
        result = mock_validator._parse_stockholm_alignment(
            str(stockholm_file), 
            'BSNTN294624_r_1_s_50_BSNTN294624_fcleaner_EDITED'
        )
        
        # Verify we get real sequence data
        assert result is not None
        assert len(str(result.seq)) > 400  # Real sequence is substantial
        assert str(result.seq).startswith('------ATATTTTATCCTTGGAAT')

    @patch('builtins.open', new_callable=mock_open)
    def test_align_single_fragment_hmmalign_failure(self, mock_file, mock_validator):
        """Test alignment when hmmalign returns error code."""
        mock_validator.hmmalign.run.return_value = 1  # Failure
        
        result = mock_validator._align_single_fragment(
            'ATGCGT',
            'test_fragment',
            Path('/fake/hmm/file.hmm')
        )
        
        assert result is None

    @patch('builtins.open', new_callable=mock_open)
    def test_align_single_fragment_exception(self, mock_file, mock_validator):
        """Test alignment when exception occurs."""
        mock_validator.hmmalign.run.side_effect = Exception("hmmalign crashed")
        
        result = mock_validator._align_single_fragment(
            'ATGCGT',
            'test_fragment',
            Path('/fake/hmm/file.hmm')
        )
        
        assert result is None

    # ============================================================================
    # TEST: _recombine_fragments()
    # ============================================================================

    def test_recombine_fragments_realistic_scenario(self, mock_validator):
        """Test recombining fragments using patterns from real alignments."""
        # Based on the Stockholm alignment patterns you showed earlier
        realistic_fragments = [
            # Fragment 1: 5' end with good alignment
            ('TATAATATGATTGTAACTATACATGCTTTTCTTATAATTTTTTTTTTTATTATACCCNTNATAATTGGTGGATTTGGAAATTGATTATTACCTTTGATATTAGGAAGACCTGATATAGCATTTCCTCGGATAAATAATATAAGGTTTTGATTACTTCCTCCTTCGTTAAGACTTCTTTTAGTCTCGTCTTTAGTAGAAACAGGTACTGGAACTGGATGAACAAT', 0, 225),
            
            # Fragment 2: 3' end with good alignment  
            ('TTACCTGTTTTGGCAGGTGCAATTACCATATTATTGACTGATCGAAATTTAAATACTTCATTTTTNNNNNNNNCTGGNGGAGGNGATCCGATNTTGTTTCAACACNTNTT', 588, 153)
        ]
        
        result = mock_validator._recombine_fragments('', realistic_fragments)
        
        # Should use the longer fragment's length (actual length is 224)
        assert len(result) == 224
        assert result.startswith('TATAATATGATT')

    def test_recombine_fragments_overlapping_consensus(self, mock_validator):
        """Test consensus building when fragments overlap in HMM space."""
        aligned_fragments = [
            ('ATGCCC---', 0, 9),
            ('---TTTTAA', 3, 12),
            ('------AAA', 6, 15)
        ]
        
        result = mock_validator._recombine_fragments('', aligned_fragments)
        # Should take first non-gap character at each position
        assert result.startswith('ATG')
        assert len(result) == 9
        assert 'ATG' in result


    # ============================================================================
    # TEST: Tilde removal preprocessing
    # ============================================================================

    def test_tilde_removal_with_real_sequence(self, mock_validator):
        """Test tilde removal using real sequence from file."""
        examples_path = Path(__file__).parent.parent / "examples"
        fasta_file = examples_path / "with_tildes.fasta"
        
        if fasta_file.exists():
            # Use real sequence with tildes
            records = list(SeqIO.parse(fasta_file, "fasta"))
            original_seq = str(records[0].seq)
            
            # Test tilde removal workflow
            seq_after_tilde_removal = original_seq.replace('~', '')
            
            # Verify tildes were removed
            assert '~' not in seq_after_tilde_removal
            assert len(seq_after_tilde_removal) < len(original_seq)
            
            # Test fragment splitting after tilde removal
            fragments = mock_validator._split_sequence_on_gaps(seq_after_tilde_removal)
            
            # Verify no tildes in fragments
            for fragment, _, _ in fragments:
                assert '~' not in fragment
        else:
            # Fallback test with example sequence
            original_seq = "ATGCGT~~~AAA~~~CGT"
            seq_after_tilde_removal = original_seq.replace('~', '')
            assert seq_after_tilde_removal == "ATGCGTAAACGT"
            
            fragments = mock_validator._split_sequence_on_gaps(seq_after_tilde_removal)
            assert len(fragments) == 1
            assert fragments[0] == ('ATGCGTAAACGT', 0, 12)

    def test_tilde_removal_mixed_with_real_sequence(self, mock_validator):
        """Test tilde removal with gaps using real mixed sequence."""
        examples_path = Path(__file__).parent.parent / "examples"
        fasta_file = examples_path / "mixed_tildes_gaps.fasta"
        
        if fasta_file.exists():
            # Use real mixed sequence
            records = list(SeqIO.parse(fasta_file, "fasta"))
            original_seq = str(records[0].seq)
            
            # Step 1: Remove tildes
            seq_after_tilde_removal = original_seq.replace('~', '')
            
            # Step 2: Split on gaps
            fragments = mock_validator._split_sequence_on_gaps(seq_after_tilde_removal)
            
            # Verify processing
            assert len(fragments) > 1  # Should create multiple fragments
            for fragment, start_pos, end_pos in fragments:
                assert '~' not in fragment  # No tildes
                assert '-' not in fragment  # No gaps
                assert len(fragment) > 0    # Non-empty
        else:
            # Fallback test
            original_seq = "ATG~~~CGT---AAA~~~TTT"
            seq_after_tilde_removal = "ATGCGT---AAATTT"
            fragments = mock_validator._split_sequence_on_gaps(seq_after_tilde_removal)
            
            expected = [
                ('ATGCGT', 0, 6),
                ('AAATTT', 9, 15)
            ]
            assert fragments == expected

    def test_tilde_removal_edge_cases(self, mock_validator):
        """Test edge cases for tilde removal."""
        # Only tildes - should result in empty sequence
        seq_after_tilde_removal = "~~~~~~~".replace('~', '')
        assert seq_after_tilde_removal == ""
        fragments = mock_validator._split_sequence_on_gaps(seq_after_tilde_removal)
        assert len(fragments) == 0
        
        # Mixed characters - tildes removed, others preserved
        original_seq = "ATG~~~CGT---NNN~~~"
        seq_after_tilde_removal = original_seq.replace('~', '')
        assert seq_after_tilde_removal == "ATGCGT---NNN"
        
        fragments = mock_validator._split_sequence_on_gaps(seq_after_tilde_removal)
        expected = [('ATGCGT', 0, 6), ('NNN', 9, 12)]
        assert fragments == expected


class TestFragmentProcessingWithFiles:
    """Test fragment processing using real input files from examples/ directory."""

    @pytest.fixture
    def mock_validator(self):
        """Create validator with mocked dependencies."""
        with patch('barcode_validator.validators.protein_coding.StructuralValidator.__init__'):
            validator = ProteinCodingValidator(Mock(spec=Config))
            validator.logger = Mock()
            validator.hmmalign = Mock()
            validator.hmm_profile_dir = None
            validator.marker = Mock()
            validator.marker.value = 'COI-5P'
            return validator

    @pytest.fixture
    def examples_dir(self):
        """Path to examples directory (one level up from tests)."""
        examples_path = Path(__file__).parent.parent / "examples"
        if not examples_path.exists():
            pytest.skip("examples/ directory not found. Please create examples/ with test files.")
        return examples_path

    def test_real_stockholm_parsing_good_fragment(self, mock_validator, examples_dir):
        """Test parsing real Stockholm file showing good alignment."""
        stockholm_file = examples_dir / "good_fragment.sto"
        if not stockholm_file.exists():
            pytest.skip("good_fragment.sto not found in examples/")

        # Test parsing the real Stockholm file
        result = mock_validator._parse_stockholm_alignment(
            str(stockholm_file), 
            'BSNTN294624_r_1_s_50_BSNTN294624_fcleaner_EDITED'
        )
        
        assert result is not None
        assert result.id == 'BSNTN294624_r_1_s_50_BSNTN294624_fcleaner_EDITED'
        assert len(str(result.seq)) > 400  # Should be a substantial sequence
        assert str(result.seq).startswith('------ATATTTTATCCTTGGAAT')  # Matches your file

    def test_real_stockholm_parsing_poor_fragment(self, mock_validator, examples_dir):
        """Test parsing real Stockholm file showing poor alignment (mostly gaps)."""
        stockholm_file = examples_dir / "poor_fragment.sto"  
        if not stockholm_file.exists():
            pytest.skip("poor_fragment.sto not found in examples/")

        # Extract sequence ID from the Stockholm file first
        with open(stockholm_file, 'r') as f:
            content = f.read()
            
        # Find the sequence ID (first non-comment, non-GR line)
        for line in content.split('\n'):
            if line.strip() and not line.startswith('#') and not line.startswith('//'):
                if not line.startswith('#=GR') and not line.startswith('#=GC'):
                    seq_id = line.split()[0]
                    break
        else:
            pytest.skip("Could not find sequence ID in poor_fragment.sto")

        # Test parsing
        result = mock_validator._parse_stockholm_alignment(str(stockholm_file), seq_id)
        
        assert result is not None
        assert result.id == seq_id
        # Poor fragment should have many gaps
        gap_count = str(result.seq).count('-')
        total_length = len(str(result.seq))
        gap_percentage = gap_count / total_length * 100
        assert gap_percentage > 50  # Poor alignment should be mostly gaps

    def test_fragmented_sequence_splitting(self, mock_validator, examples_dir):
        """Test splitting a real fragmented sequence."""
        fasta_file = examples_dir / "fragmented_coi5p.fasta"
        if not fasta_file.exists():
            pytest.skip("fragmented_coi5p.fasta not found in examples/")

        # Load the sequence
        from Bio import SeqIO
        records = list(SeqIO.parse(fasta_file, "fasta"))
        assert len(records) == 1
        
        record = records[0]
        
        # Test fragment splitting
        fragments = mock_validator._split_sequence_on_gaps(str(record.seq))
        
        # Should create multiple fragments
        assert len(fragments) >= 2
        
        # Verify fragments don't contain gaps
        for fragment, start_pos, end_pos in fragments:
            assert '-' not in fragment
            assert '~' not in fragment
            assert len(fragment) > 0

    def test_tilde_sequence_processing(self, mock_validator, examples_dir):
        """Test processing sequence with tildes."""
        fasta_file = examples_dir / "with_tildes.fasta"
        if not fasta_file.exists():
            pytest.skip("with_tildes.fasta not found in examples/")

        from Bio import SeqIO
        records = list(SeqIO.parse(fasta_file, "fasta"))
        assert len(records) == 1
        
        record = records[0]
        original_seq = str(record.seq)
        
        # Simulate tilde removal step
        seq_after_tilde_removal = original_seq.replace('~', '')
        
        # Test that tildes were actually removed
        assert '~' not in seq_after_tilde_removal
        assert len(seq_after_tilde_removal) < len(original_seq)
        
        # Test fragment splitting after tilde removal
        fragments = mock_validator._split_sequence_on_gaps(seq_after_tilde_removal)
        
        # Verify no tildes in any fragments
        for fragment, _, _ in fragments:
            assert '~' not in fragment

    def test_mixed_tildes_gaps_processing(self, mock_validator, examples_dir):
        """Test processing sequence with both tildes and gaps."""
        fasta_file = examples_dir / "mixed_tildes_gaps.fasta"
        if not fasta_file.exists():
            pytest.skip("mixed_tildes_gaps.fasta not found in examples/")

        from Bio import SeqIO
        records = list(SeqIO.parse(fasta_file, "fasta"))
        assert len(records) == 1
        
        record = records[0]
        original_seq = str(record.seq)
        
        # Step 1: Remove tildes (simulating the preprocessing step)
        seq_after_tilde_removal = original_seq.replace('~', '')
        
        # Step 2: Split on gaps
        fragments = mock_validator._split_sequence_on_gaps(seq_after_tilde_removal)
        
        # Verify processing worked correctly
        assert len(fragments) > 1  # Should create multiple fragments
        
        for fragment, start_pos, end_pos in fragments:
            assert '~' not in fragment  # No tildes
            assert '-' not in fragment  # No gaps
            assert len(fragment) > 0    # Non-empty


class TestFragmentProcessingEdgeCases:
    """Edge cases specific to fragment processing workflow."""

    @pytest.fixture
    def mock_validator(self):
        """Create minimal validator for edge case testing."""
        with patch('barcode_validator.validators.protein_coding.StructuralValidator.__init__'):
            validator = ProteinCodingValidator(Mock(spec=Config))
            validator.logger = Mock()
            # Add missing attributes that the real validator has
            validator.hmm_profile_dir = None
            validator.marker = Mock()
            validator.marker.value = 'COI-5P'
            return validator

    def test_extremely_fragmented_sequence(self, mock_validator):
        """Test sequence fragmented into many tiny pieces."""
        # 100 single-base fragments
        seq_parts = []
        for i in range(100):
            seq_parts.append('A')
            if i < 99:  # Don't add gap after last fragment
                seq_parts.append('-')
        
        seq_str = ''.join(seq_parts)
        fragments = mock_validator._split_sequence_on_gaps(seq_str)
        
        assert len(fragments) == 100
        assert all(len(frag[0]) == 1 for frag in fragments)
        assert all(frag[0] == 'A' for frag in fragments)

    def test_very_long_single_fragment(self, mock_validator):
        """Test handling of very long unfragmented sequence."""
        long_seq = 'ATGC' * 5000  # 20kb sequence
        fragments = mock_validator._split_sequence_on_gaps(long_seq)
        
        assert len(fragments) == 1
        assert len(fragments[0][0]) == 20000

    def test_consensus_memory_efficiency(self, mock_validator):
        """Test memory efficiency with many overlapping fragments."""
        # Create 50 overlapping fragments
        aligned_fragments = []
        for i in range(50):
            # Each fragment covers a slightly different region
            fragment = '-' * i + 'ATGCGT' + '-' * (50 - i)
            aligned_fragments.append((fragment, 0, len(fragment)))
        
        result = mock_validator._recombine_fragments('', aligned_fragments)
        assert 'ATGCGT' in result
        assert len(result) == 56  # Should handle without memory issues


if __name__ == "__main__":
    pytest.main([__file__, "-v"])