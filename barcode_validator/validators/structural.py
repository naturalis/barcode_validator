from Bio.SeqRecord import SeqRecord
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

    Subclasses should implement the `validate_marker_specific` method to perform
    marker-specific steps for their marker type. As a whole, this class hierarchy
    participates in the overall validation process through the appropriate instance
    from the hierarchy that has the required functionality for the marker type.
    This instance is invoked by the overall orchestrator as part of its composed
    validation steps.

    Examples:
        >>> from nbitk.config import Config
        >>> from Bio.SeqRecord import SeqRecord
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
        result.seq_length = length
        result.full_length = length  # For now these are the same
        self.logger.debug(f"Length: {length}")

    def validate_ambiguities(self, record: SeqRecord, result: DNAAnalysisResult) -> None:
        """
        Count and store ambiguous bases in the result object.

        :param record: The DNA sequence record to validate
        :param result: Result object to populate with ambiguity data
        """
        # Valid IUPAC codes (both cases)
        valid_bases = set('ACGTUacgtu-~')
        ambig_codes = set('RYKMSWBDHVNrykmswbdhvn')
        valid_chars = valid_bases | ambig_codes

        # Check for invalid characters
        seq_str = str(record.seq)
        invalid_chars = set(seq_str) - valid_chars
        if invalid_chars:
            self.logger.warning(f"Sequence {record.id} contains invalid characters: {', '.join(sorted(invalid_chars))}")

        # Count only valid IUPAC ambiguity codes
        ambig_count = sum(1 for base in seq_str if base in ambig_codes)
        result.ambiguities = ambig_count
        result.full_ambiguities = ambig_count
        self.logger.debug(f"Ambiguities: {ambig_count}")

    def validate_marker_specific(self, record: SeqRecord, result: DNAAnalysisResult) -> None:
        """
        Perform marker-specific validation steps.

        This method should be overridden by subclasses to handle validation
        specific to their marker type (e.g., codon validation for protein-coding
        or gc-content calculation for non-coding).
        Base implementation does nothing.

        :param record: The DNA sequence record to validate
        :param result: Result object to populate with marker-specific data
        """
        pass