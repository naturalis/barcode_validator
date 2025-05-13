import importlib.resources
from pathlib import Path

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
from typing import Optional, List, Tuple
import tempfile

from nbitk.Taxon import Taxon
from nbitk.config import Config
from .structural import StructuralValidator
from barcode_validator.dna_analysis_result import DNAAnalysisResult
from barcode_validator.constants import Marker


class ProteinCodingValidator(StructuralValidator):
    """
    Validator for protein-coding markers like COI-5P.

    This class extends StructuralValidator to collect protein-coding specific data,
    including sequence alignment, reading frame determination, and stop codon positions.
    These calculations form part of the structural validation process. While other
    calculations can be performed directly on the DNA sequence record, this class
    delegates most of the work to the external tool `hmmalign`, which is wrapped by the
    Naturalis bioinformatics toolkit (NBITK).

    In turn, this class exposes a single method for validation (`validate_marker_specific`).
    This method is invoked by the overall orchestrator as part of its composed validation
    process.

    Examples:
        >>> from nbitk.config import Config
        >>> config = Config()
        >>> config.load_config('/path/to/config.yaml')
        >>> validator = ProteinCodingValidator(config, '/path/to/hmm_profiles')
        >>> result = DNAAnalysisResult("foo")
        >>> validator.validate_marker_specific(SeqRecord(Seq('acgatgctacgag'),id="foo"), result, trans_table=5)

    :param config: Configuration object containing validation parameters
    """

    def __init__(self, config: Config):
        """Initialize the protein coding validator."""
        super().__init__(config)

    def validate_marker_specific(self, record: SeqRecord, result: DNAAnalysisResult) -> None:
        """
        Collect protein-coding specific measurements.
        Populates the result object with collected data.

        :param record: The DNA sequence record to validate
        :param result: Result object to populate with validation data
        """

        # Make sure that result.exp_taxon is assigned through tr.enrich_result() because this may not have been done
        # if we aren't doing taxonomic validation
        if result.exp_taxon is None:
            self.taxonomy_resolver.enrich_result(record, result)
            if result.error:
                return
        trans_table = int(self.get_translation_table(self.marker, result.exp_taxon))

        # Align sequence to HMM
        aligned_seq = self._align_sequence(record)
        if not aligned_seq:
            result.error = 'HMM alignment failed'
            return
        result.add_ancillary('nuc', str(aligned_seq.seq))

        # In the base class, both the marker length and the full length are calculated
        # in the same way and result in the same value. Here we recalculate the
        # marker length after alignment to the HMM (and trimming). Same for ambiguity
        # calculation, btw.
        result.ambiguities = self._calc_ambiguities(aligned_seq.seq)
        result.seq_length = self._calc_length(aligned_seq.seq)

        # Determine and store reading frame
        reading_frame, best_prot = self._determine_reading_frame(aligned_seq, trans_table)
        result.add_ancillary('reading_frame', str(reading_frame))
        self.logger.debug(f"Reading frame: {reading_frame}")

        # Find and store stop codon positions
        stop_positions = self._find_stop_codons(reading_frame, best_prot)
        result.stop_codons = stop_positions
        self.logger.debug(f"Stop codon positions: {stop_positions}")

    def _align_sequence(self, record: SeqRecord) -> Optional[SeqRecord]:
        """
        Align sequence to marker-specific HMM profile.

        :param record: The DNA sequence record to align
        :return: Aligned sequence record or None if alignment fails
        """
        if self.hmm_profile_dir:
            hmm_file = Path(self.hmm_profile_dir) / Path(f"{self.marker.value}.hmm")
        else:
            with importlib.resources.path("barcode_validator.hmm_files", f"{self.marker.value}.hmm") as path:
                hmm_file = path
        if not hmm_file.exists():
            self.logger.error(f"HMM profile not found: {hmm_file}")
            return None

        try:
            with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as temp_input, \
                    tempfile.NamedTemporaryFile(mode='w+', suffix='.sto', delete=False) as temp_output:

                # Remove any gaps from the sequence
                seq = record.seq
                seq = seq.replace('-', '')

                temp_input.write(f">{record.id}\n{str(seq)}\n")
                temp_input.flush()

                self.hmmalign.set_hmmfile(str(hmm_file))
                self.hmmalign.set_seqfile(temp_input.name)
                self.hmmalign.set_output(temp_output.name)
                self.hmmalign.set_outformat('Stockholm')
                self.hmmalign.set_dna()
                self.hmmalign.set_trim()

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

    @staticmethod
    def _get_complete_codons(seq, offset):
        complete_codons = ""
        # Loop over the sequence in steps of 3, starting from the given offset.
        for i in range(offset, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            # Only add the codon if it is complete and contains no gap characters.
            if '-' not in codon:
                complete_codons += codon
        return complete_codons

    def _determine_reading_frame(self, aligned_seq: SeqRecord, trans_table: int) -> Tuple[int, SeqRecord]:
        """
        Determine the reading frame from HMM-aligned sequence.

        :param aligned_seq: HMM-aligned sequence
        :return: Reading frame (0, 1, or 2)
        """
        raw_seq = aligned_seq.seq
        best_frame = 0
        best_prot = ""
        min_stops = float('inf')

        # Brute force the 3 possible phases by counting
        # stop codons for each offset
        for frame in range(3):
            coding_seq = Seq("".join(self._get_complete_codons(raw_seq, frame)))
            protein = coding_seq.translate(table=trans_table)
            stops = protein.count('*')
            self.logger.debug(f"Frame {frame}: {stops} stop codons")

            if stops < min_stops:
                min_stops = stops
                best_frame = frame
                best_prot = protein

        self.logger.debug(f"Best reading frame: {best_frame}")
        return best_frame, best_prot

    def _find_stop_codons(self, reading_frame, protein) -> List[int]:
        """
        Find positions of stop codons in the sequence.

        :param reading_frame: The reading frame of the protein sequence (0, 1, or 2)
        :param protein: The protein sequence record
        :return: List of stop codon positions
        """
        stop_positions = []
        for i, aa in enumerate(protein[:-1]):
            if aa == '*':
                nuc_pos = (i * 3) + reading_frame
                stop_positions.append(nuc_pos)

        self.logger.debug(f"Stop codon positions: {stop_positions}")
        return stop_positions

    def get_translation_table(self, marker: Marker, taxon: Taxon) -> int:
        """
        Determine the appropriate translation table based on marker and taxonomy.

        :param marker: The genetic Marker (enum) being analyzed
        :param taxon: The NCBI taxon object representing at the highest a family in the taxonomy.
        :return: Translation table index (int) for use with Biopython
        """
        taxonomy_dict = self.taxonomy_resolver.get_lineage_dict(taxon)

        if marker in [Marker.MATK, Marker.RBCL]:
            if taxonomy_dict.get('family') == 'Balanophoraceae':
                return 32
            return 11

        elif marker == Marker.COI_5P:
            phylum = taxonomy_dict.get('phylum')
            tax_class = taxonomy_dict.get('class')
            family = taxonomy_dict.get('family')

            if phylum == 'Chordata':
                if tax_class == 'Ascidiacea':
                    return 13
                elif tax_class in ['Actinopteri', 'Amphibia', 'Mammalia', 'Aves', 'Reptilia']:
                    return 2

            elif phylum == 'Hemichordata':
                if family == 'Cephalodiscidae':
                    return 33
                elif family == 'Rhabdopleuridae':
                    return 24

            elif phylum in ['Echinodermata', 'Platyhelminthes']:
                return 9

            # Default invertebrate mitochondrial code for other invertebrates
            if phylum != 'Chordata':
                return 5

        # Fall back to standard code if no specific rules match
        return 1

    @staticmethod
    def requires_resolver() -> bool:
        return True

    @staticmethod
    def requires_marker() -> bool:
        return True

    @staticmethod
    def requires_hmmalign() -> bool:
        return True