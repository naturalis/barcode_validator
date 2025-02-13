from Bio.SeqRecord import SeqRecord
from typing import Dict, Tuple
from nbitk.config import Config
from barcode_validator.structural_validator import StructuralValidator
from barcode_validator.taxonomy_resolver import Marker


class NonCodingValidator(StructuralValidator):
    """
    Validator for non-coding markers like ITS.

    This class extends StructuralValidator to add non-coding specific validation
    such as sequence composition checks.

    Examples:
        >>> from nbitk.config import Config
        >>> config = Config()
        >>> config.load_config('/path/to/config.yaml')
        >>> validator = NonCodingValidator(config)
        >>> record = SeqRecord(...)
        >>> is_valid, details = validator.validate_marker_specific(record)

    :param config: Configuration object containing validation parameters
    """

    def __init__(self, config: Config):
        """Initialize the non-coding validator."""
        super().__init__(config)
        self.marker = Marker(config.get('marker', 'ITS'))

    def validate_marker_specific(self, record: SeqRecord) -> Tuple[bool, Dict]:
        """
        Validate non-coding specific features including sequence composition.

        :param record: The DNA sequence record to validate
        :return: Tuple of (validation_success, validation_details)
        """
        validation_results = {}

        # Verify sequence composition
        gc_content = self.calculate_gc_content(record)
        validation_results['gc_content'] = gc_content

        # For now, just check basic composition
        is_valid = True

        validation_results['is_valid'] = is_valid

        return is_valid, validation_results

    def calculate_gc_content(self, record: SeqRecord) -> float:
        """
        Calculate GC content of the sequence.

        :param record: The DNA sequence record
        :return: GC content as percentage
        """
        seq = str(record.seq).upper()
        gc_count = seq.count('G') + seq.count('C')
        total = len(seq)
        return (gc_count / total * 100) if total > 0 else 0
