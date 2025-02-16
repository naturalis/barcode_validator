from Bio.SeqRecord import SeqRecord
from typing import Dict, Tuple
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from barcode_validator.dna_analysis_result import DNAAnalysisResult


class StructuralValidator:
    """
    Base class for structural validation of DNA sequences.

    This class defines the interface and common functionality for validating
    the structural properties of DNA sequences, such as length and ambiguous
    bases. Specific marker types (protein-coding vs non-coding) implement
    their own validators inheriting from this class.

    Examples:
        >>> from nbitk.config import Config
        >>> config = Config()
        >>> config.load_config('/path/to/config.yaml')
        >>> validator = StructuralValidator(config)
        >>> result = DNAAnalysisResult("sequence_id")
        >>> validator.validate_sequence(record, result)

    :param config: Configuration object containing validation parameters
    """

    def __init__(self, config: Config):
        """Initialize the structural validator."""
        self.config = config
        self.logger = get_formatted_logger(self.__class__.__name__, config)

    def validate_sequence(self, record: SeqRecord, result: DNAAnalysisResult) -> None:
        """
        Validate a DNA sequence record structurally. Populates the provided
        result object with validation outcomes.

        :param record: The DNA sequence record to validate
        :param result: DNAAnalysisResult object to populate with validation results
        """
        # Validate sequence length
        self.validate_length(record, result)

        # Validate ambiguous bases
        self.validate_ambiguities(record, result)

        # Perform marker-specific validation
        self.validate_marker_specific(record, result)

    def validate_length(self, record: SeqRecord, result: DNAAnalysisResult) -> None:
        """
        Validate the sequence length against minimum requirements. Stores the
        sequence length in the result object.

        :param record: The DNA sequence record to validate
        :param result: Result object to populate with length data
        """
        seq_str = str(record.seq).replace('-', '').replace('~', '')
        length = len(seq_str)

        # Store lengths in result object
        result.nuc_basecount = length
        result.nuc_full_basecount = length  # For now these are the same
        self.logger.debug(f"Length: {length}")

    def validate_ambiguities(self, record: SeqRecord, result: DNAAnalysisResult) -> None:
        """
        Count and store ambiguous bases in the result object.

        :param record: The DNA sequence record to validate
        :param result: Result object to populate with ambiguity data
        """
        # Count ambiguous bases (anything not ACGT)
        ambig_count = len([base for base in record.seq if base not in 'acgtACGT-~'])

        # Store counts in result object
        result.ambig_basecount = ambig_count
        result.ambig_full_basecount = ambig_count
        self.logger.debug(f"Ambiguities: {ambig_count}")

    def validate_marker_specific(self, record: SeqRecord, result: DNAAnalysisResult) -> None:
        """
        Perform marker-specific validation steps.

        This method should be overridden by subclasses to handle validation
        specific to their marker type (e.g., codon validation for protein-coding).
        Base implementation does nothing.

        :param record: The DNA sequence record to validate
        :param result: Result object to populate with marker-specific data
        """
        pass