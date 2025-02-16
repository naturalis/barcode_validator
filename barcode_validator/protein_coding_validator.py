from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
from typing import Optional, List
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

    This class extends StructuralValidator to collect protein-coding specific data,
    including sequence alignment, reading frame determination, and stop codon positions.

    Examples:
        >>> from nbitk.config import Config
        >>> config = Config()
        >>> config.load_config('/path/to/config.yaml')
        >>> validator = ProteinCodingValidator(config, '/path/to/hmm_profiles')
        >>> result = DNAAnalysisResult("sequence_id")
        >>> validator.validate_marker_specific(record, result, trans_table=5)

    :param config: Configuration object containing validation parameters
    :param hmm_profile_dir: Directory containing HMM profiles for markers
    """

    def __init__(self, config: Config, hmm_profile_dir: str):
        """Initialize the protein coding validator."""
        super().__init__(config)
        self.hmm_profile_dir = Path(hmm_profile_dir)
        self.marker = Marker(config.get('marker', 'COI-5P'))
        self.hmmalign = Hmmalign(config)

    def validate_marker_specific(self, record: SeqRecord, result: DNAAnalysisResult) -> None:
        """
        Collect protein-coding specific measurements.
        Populates the result object with collected data.

        :param record: The DNA sequence record to validate
        :param result: Result object to populate with validation data
        """
        # Get translation table from result object
        trans_table = result.ancillary.get('translation_table')
        if trans_table is None:
            result.error = 'Translation table not available'
            return
        trans_table = int(trans_table)  # Convert from string since stored as ancillary data

        # Align sequence to HMM
        aligned_seq = self._align_sequence(record)
        if not aligned_seq:
            result.error = 'HMM alignment failed'
            return

        # Determine and store reading frame
        reading_frame = self._determine_reading_frame(aligned_seq)
        result.add_ancillary('reading_frame', str(reading_frame))
        self.logger.debug(f"Reading frame: {reading_frame}")

        # Find and store stop codon positions
        stop_positions = self._find_stop_codons(aligned_seq, trans_table, reading_frame)
        result.stop_codons = stop_positions
        self.logger.debug(f"Stop codon positions: {stop_positions}")

    def _align_sequence(self, record: SeqRecord) -> Optional[SeqRecord]:
        """
        Align sequence to marker-specific HMM profile.

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

                temp_input.write(f">{record.id}\n{str(record.seq)}\n")
                temp_input.flush()

                self.hmmalign.set_hmmfile(str(hmm_file))
                self.hmmalign.set_seqfile(temp_input.name)
                self.hmmalign.set_output(temp_output.name)
                self.hmmalign.set_outformat('Stockholm')
                self.hmmalign.set_dna()

                return_code = self.hmmalign.run()
                if return_code != 0:
                    self.logger.error("hmmalign failed")
                    return None

                return self._parse_stockholm_alignment(temp_output.name, record.id)

        except Exception as e:
            self.logger.error(f"Error in HMM alignment: {str(e)}")
            return None

    def _parse_stockholm_alignment(self, stockholm_file: str, seq_id: str) -> Optional[SeqRecord]:
        """
        Parse Stockholm format alignment output.

        :param stockholm_file: Path to Stockholm format alignment file
        :param seq_id: ID of the sequence to extract
        :return: Aligned sequence record or None if parsing fails
        """
        try:
            alignment = AlignIO.read(stockholm_file, "stockholm")
            for record in alignment:
                if record.id == seq_id:
                    return record
            self.logger.error(f"Sequence {seq_id} not found in alignment")
            return None
        except Exception as e:
            self.logger.error(f"Error parsing Stockholm alignment: {str(e)}")
            return None

    def _determine_reading_frame(self, aligned_seq: SeqRecord) -> int:
        """
        Determine the reading frame from HMM-aligned sequence.

        :param aligned_seq: HMM-aligned sequence
        :return: Reading frame (0, 1, or 2)
        """
        seq_nogaps = str(aligned_seq.seq).replace('-', '')
        best_frame = 0
        min_stops = float('inf')

        for frame in range(3):
            coding_seq = Seq(seq_nogaps[frame:])
            protein = coding_seq.translate(table=1)
            stops = protein.count('*')

            if stops < min_stops:
                min_stops = stops
                best_frame = frame

        self.logger.debug(f"Best reading frame: {best_frame}")
        return best_frame

    def _find_stop_codons(self, record: SeqRecord, trans_table: int, reading_frame: int) -> List[int]:
        """
        Find positions of stop codons in the sequence.

        :param record: The DNA sequence record
        :param trans_table: Translation table number
        :param reading_frame: Reading frame (0, 1, or 2)
        :return: List of stop codon positions
        """
        seq = str(record.seq).replace('-', '')
        coding_seq = Seq(seq[reading_frame:])
        protein = coding_seq.translate(table=trans_table)

        stop_positions = []
        for i, aa in enumerate(protein[:-1]):
            if aa == '*':
                nuc_pos = (i * 3) + reading_frame
                stop_positions.append(nuc_pos)

        self.logger.debug(f"Stop codon positions: {stop_positions}")
        return stop_positions