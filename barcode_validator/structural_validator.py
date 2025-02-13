from Bio.SeqRecord import SeqRecord
from typing import Dict, Tuple
from nbitk.config import Config
from nbitk.logger import get_formatted_logger


class StructuralValidator:
    """
    Base class for structural validation of DNA sequences.
    
    This class defines the interface and common functionality for validating
    the structural properties of DNA sequences, such as length and ambiguous
    bases. Specific marker types (protein-coding vs non-coding) should
    implement their own subclasses.

    Examples:
        >>> from nbitk.config import Config
        >>> config = Config()
        >>> config.load_config('/path/to/config.yaml')
        >>> validator = StructuralValidator(config)
        >>> record = SeqRecord(...)
        >>> is_valid, details = validator.validate_sequence(record)

    :param config: Configuration object containing validation parameters
    """

    def __init__(self, config: Config):
        """Initialize the structural validator."""
        self.config = config
        self.logger = get_formatted_logger(self.__class__.__name__, config)
        
    def validate_sequence(self, record: SeqRecord) -> Tuple[bool, Dict]:
        """
        Validate a DNA sequence record structurally.
        
        :param record: The DNA sequence record to validate            
        :return: Tuple of (validation_success, validation_details)
        """
        results = {}
        
        # Validate sequence length
        length_valid, length_details = self.validate_length(record)
        results.update(length_details)
        
        # Validate ambiguous bases
        ambig_valid, ambig_details = self.validate_ambiguities(record)
        results.update(ambig_details)
        
        # Perform marker-specific validation
        marker_valid, marker_details = self.validate_marker_specific(record)
        results.update(marker_details)
        
        # Overall validation success requires all checks to pass
        success = all([length_valid, ambig_valid, marker_valid])
        
        return success, results

    def validate_length(self, record: SeqRecord) -> Tuple[bool, Dict]:
        """
        Validate the sequence length against minimum requirements.
        
        :param record: The DNA sequence record to validate            
        :return: Tuple of (length_valid, length_details)
        """
        min_length = self.config.get('min_sequence_length', 500)
        seq_length = len(record.seq)
        
        is_valid = seq_length >= min_length
        details = {
            'sequence_length': seq_length,
            'min_length_required': min_length
        }
        
        return is_valid, details

    def validate_ambiguities(self, record: SeqRecord) -> Tuple[bool, Dict]:
        """
        Count and validate ambiguous bases in the sequence.
        
        :param record: The DNA sequence record to validate            
        :return: Tuple of (ambiguities_valid, ambiguity_details)
        """
        max_ambiguities = self.config.get('max_ambiguities', 6)
        ambig_count = sum(1 for base in record.seq if base not in 'ATCGatcg')
        
        is_valid = ambig_count <= max_ambiguities
        details = {
            'ambiguous_bases': ambig_count,
            'max_ambiguities_allowed': max_ambiguities
        }
        
        return is_valid, details

    def validate_marker_specific(self, record: SeqRecord) -> Tuple[bool, Dict]:
        """
        Perform marker-specific validation steps.
        
        This method should be overridden by subclasses to handle validation
        specific to their marker type (e.g., codon validation for protein-coding).
        Base implementation always returns valid.
        
        :param record: The DNA sequence record to validate            
        :return: Tuple of (marker_valid, marker_details)
        """
        return True, {}

    def validate_length(self, record: SeqRecord) -> Tuple[bool, Dict]:
        """
        Using length calculation from sequence_handler.py
        """
        seq_str = str(record.seq).replace('-', '').replace('~', '')
        length = len(seq_str)
        return length >= self.config.get('min_length', 500), {'length': length}

    def validate_ambiguities(self, record: SeqRecord) -> Tuple[bool, Dict]:
        """
        Using ambiguity counting from sequence_handler.py
        """
        ambig_count = len([base for base in record.seq if base not in 'acgtACGT-~'])
        return ambig_count <= self.config.get('max_ambiguities', 6), {'ambiguities': ambig_count}