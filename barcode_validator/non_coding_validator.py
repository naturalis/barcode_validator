from Bio.SeqRecord import SeqRecord
from nbitk.config import Config
from barcode_validator.structural_validator import StructuralValidator
from barcode_validator.taxonomy_resolver import Marker
from barcode_validator.dna_analysis_result import DNAAnalysisResult


class NonCodingValidator(StructuralValidator):
    """
    Validator for non-coding markers like ITS.

    This class extends StructuralValidator to collect non-coding specific
    measurements such as sequence composition data.

    Examples:
        >>> from nbitk.config import Config
        >>> config = Config()
        >>> config.load_config('/path/to/config.yaml')
        >>> validator = NonCodingValidator(config)
        >>> result = DNAAnalysisResult("sequence_id")
        >>> validator.validate_marker_specific(record, result)

    :param config: Configuration object containing validation parameters
    """

    def __init__(self, config: Config):
        """Initialize the non-coding validator."""
        super().__init__(config)
        self.marker = Marker(config.get('marker', 'ITS'))

    def validate_marker_specific(self, record: SeqRecord, result: DNAAnalysisResult) -> None:
        """
        Collect non-coding specific measurements.
        Populates the result object with collected data.

        :param record: The DNA sequence record to validate
        :param result: Result object to populate with validation data
        """
        # Calculate and store composition metrics
        gc_content = self.calculate_gc_content(record)
        result.add_ancillary('gc_content', str(gc_content))
        self.logger.debug(f"GC content: {gc_content}%")

    def calculate_gc_content(self, record: SeqRecord) -> float:
        """
        Calculate GC content of the sequence.

        :param record: The DNA sequence record
        :return: GC content as percentage
        """
        seq = str(record.seq).upper()
        gc_count = seq.count('G') + seq.count('C')
        total = len(seq)
        gc_perc = (gc_count / total * 100) if total > 0 else 0
        self.logger.debug(f"GC percentage: {gc_perc}%")
        return gc_perc


