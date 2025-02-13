from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
from typing import Optional, List, Tuple
import tempfile
from pathlib import Path
from nbitk.config import Config
from nbitk.Tools import Hmmalign
from barcode_validator.structural_validator import StructuralValidator
from barcode_validator.taxonomy_resolver import Marker
from barcode_validator.dna_analysis_result import DNAAnalysisResult


class ProteinCodingValidator(StructuralValidator):
    """
    Validator for protein-coding markers like COI-5P.

    This class extends StructuralValidator to add protein-coding specific validation,
    including translation table determination, HMM alignment, and stop codon checking.

    Examples:
        >>> from nbitk.config import Config
        >>> config = Config()
        >>> config.load_config('/path/to/config.yaml')
        >>> validator = ProteinCodingValidator(config, '/path/to/hmm_profiles')
        >>> record = SeqRecord(...)
        >>> is_valid, details = validator.validate_marker_specific(record)

    :param config: Configuration object containing validation parameters
    :param hmm_profile_dir: Directory containing HMM profiles for markers
    """

    def __init__(self, config: Config, hmm_profile_dir: str):
        """Initialize the protein coding validator."""
        super().__init__(config)
        self.hmm_profile_dir = Path(hmm_profile_dir)
        self.marker = Marker(config.get('marker', 'COI-5P'))
        self.hmmalign = Hmmalign(config)

    def validate_marker_specific(self, record: SeqRecord) -> Tuple[bool, DNAAnalysisResult]:
        """
        Validate protein-coding specific features including reading frame and stop codons.

        :param record: The DNA sequence record to validate
        :return: Tuple of (validation_success, result_object)
        """
        result = DNAAnalysisResult(record.id)

        # Step 1: Align to marker-specific HMM
        aligned_seq = self._align_sequence(record)
        if aligned_seq is None:
            result.error = 'HMM alignment failed'
            return False, result

        # Step 2: Determine reading frame
        reading_frame = self._determine_reading_frame(aligned_seq)
        result.add_ancillary('reading_frame', str(reading_frame))

        # Step 3: Get translation table
        trans_table = self._determine_translation_table(record)
        result.add_ancillary('translation_table', str(trans_table))

        # Step 4: Translate and check for stop codons
        stop_codons = self._find_stop_codons(aligned_seq, trans_table, reading_frame)
        result.stop_codons = stop_codons

        # A valid protein-coding sequence should have no internal stop codons
        is_valid = len(stop_codons) == 0
        result.add_ancillary('is_valid', str(is_valid))

        return is_valid, result

    def _align_sequence(self, record: SeqRecord) -> Optional[SeqRecord]:
        """
        Align sequence to marker-specific HMM profile using nbitk's Hmmalign wrapper.

        :param record: The DNA sequence record to align
        :return: Aligned sequence record or None if alignment fails
        """
        hmm_file = self.hmm_profile_dir / f"{self.marker.value}.hmm"
        if not hmm_file.exists():
            self.logger.error(f"HMM profile not found: {hmm_file}")
            return None

        try:
            with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as temp_input, \
                    tempfile.NamedTemporaryFile(mode='w+', suffix='.sto', delete=False) as temp_output:

                # Write sequence to temporary FASTA
                temp_input.write(f">{record.id}\n{str(record.seq)}\n")
                temp_input.flush()

                # Configure hmmalign
                self.hmmalign.set_hmmfile(str(hmm_file))
                self.hmmalign.set_seqfile(temp_input.name)
                self.hmmalign.set_output(temp_output.name)
                self.hmmalign.set_outformat('Stockholm')
                self.hmmalign.set_dna()

                # Run hmmalign
                return_code = self.hmmalign.run()
                if return_code != 0:
                    self.logger.error("hmmalign failed")
                    return None

                # Parse Stockholm output
                aligned_seq = self._parse_stockholm_alignment(temp_output.name, record.id)
                return aligned_seq

        except Exception as e:
            self.logger.error(f"Error in HMM alignment: {str(e)}")
            return None

    def _parse_stockholm_alignment(self, stockholm_file: str, seq_id: str) -> Optional[SeqRecord]:
        """
        Parse Stockholm format alignment output from HMMER using BioPython's AlignIO.

        :param stockholm_file: Path to Stockholm format alignment file
        :param seq_id: ID of the sequence to extract
        :return: Aligned sequence record or None if parsing fails
        """
        try:
            # Parse the alignment using AlignIO
            alignment = AlignIO.read(stockholm_file, "stockholm")

            # Find the record with matching ID
            for record in alignment:
                if record.id == seq_id:
                    return record

            self.logger.error(f"Sequence {seq_id} not found in alignment")
            return None

        except Exception as e:
            self.logger.error(f"Error parsing Stockholm alignment: {str(e)}")
            return None

    def _determine_translation_table(self, record: SeqRecord) -> int:
        """
        Determine the appropriate translation table for the sequence based on marker and taxonomy.

        :param record: The DNA sequence record
        :return: Translation table number
        """
        # Use configured translation table if specified
        if 'translation_table' in self.config:
            return self.config.get('translation_table')

        # Get taxonomy from record annotations
        taxonomy = record.annotations.get('taxonomy', [])

        if self.marker in [Marker.MATK, Marker.RBCL]:
            # Handle chloroplast markers
            return 11  # Bacterial and plant plastid code
        elif self.marker == Marker.COI_5P:
            # Handle mitochondrial markers
            if 'Vertebrata' in taxonomy:
                return 2  # Vertebrate mitochondrial code
            elif 'Ascidiacea' in taxonomy:
                return 13  # Ascidian mitochondrial code
            return 5  # Invertebrate mitochondrial code

        # Default to standard code
        return 1

    def _determine_reading_frame(self, aligned_seq: SeqRecord) -> int:
        """
        Determine the reading frame from HMM-aligned sequence.

        :param aligned_seq: HMM-aligned sequence
        :return: Reading frame (0, 1, or 2)
        """
        # Remove gaps from aligned sequence
        seq_nogaps = str(aligned_seq.seq).replace('-', '')

        # Try each reading frame and score them
        best_frame = 0
        min_stops = float('inf')

        for frame in range(3):
            # Translate in this frame
            coding_seq = Seq(seq_nogaps[frame:])
            protein = coding_seq.translate(table=1)  # Use standard code for frame detection
            stops = protein.count('*')

            if stops < min_stops:
                min_stops = stops
                best_frame = frame

        return best_frame

    def _find_stop_codons(self, record: SeqRecord, trans_table: int, reading_frame: int) -> List[int]:
        """
        Find positions of premature stop codons in the sequence.

        :param record: The DNA sequence record
        :param trans_table: Translation table number
        :param reading_frame: Reading frame (0, 1, or 2)
        :return: List of stop codon positions
        """
        seq = str(record.seq).replace('-', '')  # Remove gaps
        coding_seq = Seq(seq[reading_frame:])

        stop_positions = []
        protein = coding_seq.translate(table=trans_table)

        # Find all stop codons except the last position
        for i, aa in enumerate(protein[:-1]):
            if aa == '*':
                # Convert amino acid position to nucleotide position
                nuc_pos = (i * 3) + reading_frame
                stop_positions.append(nuc_pos)

        return stop_positions

    def _translate_sequence(self, record: SeqRecord, trans_table: int) -> SeqRecord:
        """
        For COI-5P, remove first base to ensure proper phasing
        """
        # Take phasing logic from sequence_handler.py
        if self.marker == Marker.COI_5P:
            record = record[1:]  # Remove first base for proper phasing
        return record.translate(table=trans_table)