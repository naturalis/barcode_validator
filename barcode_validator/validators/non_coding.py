from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from nbitk.config import Config
from .structural import StructuralValidator
from barcode_validator.dna_analysis_result import DNAAnalysisResult


class NonCodingValidator(StructuralValidator):
    """
    Validator for non-coding markers like ITS.

    This class extends StructuralValidator to collect non-coding specific
    measurements such as sequence composition data. At present, this class
    is not needed in production just yet because so far we have only COI-5P
    sequences, but this is likely to change in the future. The only concrete
    method that is therefore currently here specifically for non-coding genes
    is GC content calculation. Further implementation of this class is
    currently on hold until ITS data starts to arrive.

    Examples:
        >>> from nbitk.config import Config
        >>> config = Config()
        >>> config.load_config('/path/to/config.yaml')
        >>> validator = NonCodingValidator(config)
        >>> result = DNAAnalysisResult("foo")
        >>> validator.validate_marker_specific(SeqRecord(Seq('acgtacgatgcatct'),id="foo"), result)

    :param config: Configuration object containing validation parameters
    """

    def __init__(self, config: Config):
        """Initialize the non-coding validator."""
        super().__init__(config)

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

    # TODO determine validation criteria for non-coding genes
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

