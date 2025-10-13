from typing import Dict, Any, List, Optional

from nbitk.config import Config

from barcode_validator.constants import Marker


class MarkerCriteria:
    """
    Base class for marker-specific validation criteria.

    This class defines the interface and default values for all marker validation criteria.
    Subclasses can override these defaults to provide marker-specific validation requirements.

    The criteria defined here match those used in the DNAAnalysisResult class for adjudicating
    whether a sequence is valid or not.

    Examples:
        >>> criteria = MarkerCriteria()
        >>> criteria.min_length
        500
        >>> criteria.max_ambiguities
        0
        >>> criteria.max_stop_codons
        0
    """

    def __init__(self, config: Config = None):
        """
        Initialize marker criteria with optional configuration data.

        :param config: Optional Config object containing configuration data to override defaults
        """
        self.min_length: int = 500  # Default minimum sequence length
        self.max_ambiguities: int = 0  # Default maximum allowed ambiguous bases
        self.max_stop_codons: int = 0  # Default maximum allowed stop codons
        self.marker_type: Marker = None  # Will be set by subclasses

        # Override defaults with provided configuration
        if config:
            self.update_from_config(config)

    def update_from_config(self, config: Config) -> None:
        """
        Update criteria from configuration data.

        :param config: Config object containing configuration data
        """
        if 'criteria' in config.config_data:
            criteria = config.config_data['criteria']

            if 'min_length' in criteria:
                self.min_length = int(criteria.get('min_length'))
            if 'max_ambiguities' in criteria:
                self.max_ambiguities = int(criteria.get('max_ambiguities'))
            if 'max_stop_codons' in criteria:
                self.max_stop_codons = int(criteria.get('max_stop_codons'))

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert criteria to a dictionary for configuration or serialization.

        :return: Dictionary representation of the criteria
        """
        return {
            'min_length': self.min_length,
            'max_ambiguities': self.max_ambiguities,
            'max_stop_codons': self.max_stop_codons,
            'marker_type': self.marker_type.value if self.marker_type else None
        }

    def __str__(self) -> str:
        """
        String representation of the criteria.

        :return: String representation
        """
        marker = self.marker_type.value if self.marker_type else "Unknown"
        return (f"{marker} Criteria: min_length={self.min_length}, "
                f"max_ambiguities={self.max_ambiguities}, "
                f"max_stop_codons={self.max_stop_codons}")


class COI5PCriteria(MarkerCriteria):
    """
    Validation criteria specific to COI-5P marker.

    This class defines the validation criteria for the COI-5P marker, which is a protein-coding
    mitochondrial marker commonly used for DNA barcoding in animals. The criteria are based on
    BOLD's requirements for COI-5P sequences.

    Examples:
        >>> criteria = COI5PCriteria()
        >>> criteria.min_length
        500
        >>> criteria.max_ambiguities
        0
        >>> criteria.max_stop_codons
        0
    """

    def __init__(self, config: Config = None):
        """
        Initialize COI-5P criteria with BOLD-specific defaults.

        :param config: Optional dictionary containing configuration data to override defaults
        """
        super().__init__(config)
        self.marker_type = Marker.COI_5P
        # BOLD defaults for COI-5P:
        self.min_length = 500  # 500 base pairs
        if config is not None:
            self.max_ambiguities = config.get("ambiguities", 0)  # No more than 0 ambiguous bases
        self.max_stop_codons = 0  # No stop codons allowed

        # Override with any provided configuration
        if config:
            self.update_from_config(config)


class MATKCriteria(MarkerCriteria):
    """
    Validation criteria specific to matK marker.

    This class defines the validation criteria for the matK marker, which is a protein-coding
    chloroplast marker commonly used for DNA barcoding in plants.

    Examples:
        >>> criteria = MATKCriteria()
        >>> criteria.min_length
        700
        >>> criteria.max_ambiguities
        20
        >>> criteria.max_stop_codons
        0
    """

    def __init__(self, config: Config = None):
        """
        Initialize matK criteria with plant-specific defaults.

        :param config: Optional Config object containing configuration data to override defaults
        """
        super().__init__(config)
        self.marker_type = Marker.MATK
        # Default criteria for matK (these are example values)
        self.min_length = 700  # 700 base pairs
        self.max_ambiguities = 20  # Plant sequences often have more ambiguities
        self.max_stop_codons = 0  # No stop codons allowed

        # Override with any provided configuration
        if config:
            self.update_from_config(config)


class RBCLCriteria(MarkerCriteria):
    """
    Validation criteria specific to rbcL marker.

    This class defines the validation criteria for the rbcL marker, which is a protein-coding
    chloroplast marker commonly used for DNA barcoding in plants.

    Examples:
        >>> criteria = RBCLCriteria()
        >>> criteria.min_length
        550
        >>> criteria.max_ambiguities
        15
        >>> criteria.max_stop_codons
        0
    """

    def __init__(self, config: Config = None):
        """
        Initialize rbcL criteria with plant-specific defaults.

        :param config: Optional Config object containing configuration data to override defaults
        """
        super().__init__(config)
        self.marker_type = Marker.RBCL
        # Default criteria for rbcL (these are example values)
        self.min_length = 200  # min length for BOLD
        self.max_ambiguities = 0  # Skims must have no ambiguities
        self.max_stop_codons = 0  # No stop codons allowed

        # Override with any provided configuration
        if config:
            self.update_from_config(config)


class ITSCriteria(MarkerCriteria):
    """
    Validation criteria specific to ITS marker.

    This class defines the validation criteria for the ITS marker, which is a non-coding
    nuclear marker commonly used for DNA barcoding in fungi.

    Examples:
        >>> criteria = ITSCriteria()
        >>> criteria.min_length
        500
        >>> criteria.max_ambiguities
        30
    """

    def __init__(self, config: Config = None):
        """
        Initialize ITS criteria with fungi-specific defaults.

        :param config: Optional Config object containing configuration data to override defaults
        """
        super().__init__(config)
        self.marker_type = Marker.ITS
        # Default criteria for ITS (these are example values)
        self.min_length = 500  # 500 base pairs
        self.max_ambiguities = 30  # ITS can have many ambiguities
        # Stop codons are irrelevant for non-coding regions
        self.max_stop_codons = None

        # Override with any provided configuration
        if config:
            self.update_from_config(config)


class ITS2Criteria(MarkerCriteria):
    """
    Validation criteria specific to ITS2 marker.

    This class defines the validation criteria for the ITS2 marker, which is a non-coding
    nuclear marker commonly used for DNA barcoding in fungi.

    Examples:
        >>> criteria = ITS2Criteria()
        >>> criteria.min_length
        350
        >>> criteria.max_ambiguities
        20
    """

    def __init__(self, config: Config = None):
        """
        Initialize ITS2 criteria with fungi-specific defaults.

        :param config: Optional Config object containing configuration data to override defaults
        """
        super().__init__(config)
        self.marker_type = Marker.ITS2
        # Default criteria for ITS2 (these are example values)
        self.min_length = 350  # 350 base pairs
        self.max_ambiguities = 20  # ITS2 can have many ambiguities
        # Stop codons are irrelevant for non-coding regions
        self.max_stop_codons = None

        # Override with any provided configuration
        if config:
            self.update_from_config(config)


class MarkerCriteriaFactory:
    """
    Factory class to create marker-specific criteria objects.

    This factory class provides methods to create instances of the appropriate MarkerCriteria
    subclass based on the marker type.

    Examples:
        >>> factory = MarkerCriteriaFactory()
        >>> criteria = factory.get_criteria(Marker.COI_5P)
        >>> isinstance(criteria, COI5PCriteria)
        True
    """

    @staticmethod
    def get_criteria(marker_type: Marker, config: Config = None) -> MarkerCriteria:
        """
        Get the appropriate criteria object for the given marker type.

        :param marker_type: The marker type enum
        :param config: Optional Config object to override defaults
        :return: An instance of the appropriate MarkerCriteria subclass
        :raises ValueError: If the marker type is unknown
        """
        if marker_type == Marker.COI_5P:
            return COI5PCriteria(config)
        elif marker_type == Marker.MATK:
            return MATKCriteria(config)
        elif marker_type == Marker.RBCL:
            return RBCLCriteria(config)
        elif marker_type == Marker.ITS:
            return ITSCriteria(config)
        elif marker_type == Marker.ITS2:
            return ITS2Criteria(config)
        else:
            raise ValueError(f"Unknown marker type: {marker_type}")

    @staticmethod
    def from_string(marker_string: str, config: Config = None) -> MarkerCriteria:
        """
        Get the appropriate criteria object from a marker string.

        :param marker_string: The marker name as a string
        :param config: Optional configuration data to override defaults
        :return: An instance of the appropriate MarkerCriteria subclass
        :raises ValueError: If the marker string is unknown
        """
        try:
            marker_type = Marker(marker_string)
            return MarkerCriteriaFactory.get_criteria(marker_type, config)
        except ValueError:
            raise ValueError(f"Unknown marker string: {marker_string}")

    @staticmethod
    def get_available_markers() -> List[str]:
        """
        Get a list of available marker types.

        :return: List of marker type names
        """
        return [marker.value for marker in Marker]
