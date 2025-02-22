from pathlib import Path
from typing import Optional

from Bio.SeqRecord import SeqRecord
from nbitk.Tools import Hmmalign
from nbitk.config import Config
from nbitk.logger import get_formatted_logger

from barcode_validator.dna_analysis_result import DNAAnalysisResult
from barcode_validator.idservices.idservice import IDService
from barcode_validator.resolvers.taxonomy import TaxonResolver
from barcode_validator.constants import Marker


class AbstractValidator:
    """
    Base class for all validators. Each validator can optionally declare whether it needs a
    taxonomy resolver and/or a marker by implementing requires_resolver() and requires_marker().
    """
    def __init__(self, config: Config):
        self.logger = get_formatted_logger(self.__class__.__name__, config)
        self.marker: Marker = Optional[Marker]
        self.taxonomy_resolver: TaxonResolver = Optional[TaxonResolver]
        self.idservice: IDService = Optional[IDService]
        self.hmmalign: Hmmalign = Optional[Hmmalign]
        self.hmm_profile_dir: Path = Optional[Path]

    def set_hmmalign(self, hmmalign: Hmmalign) -> None:
        """Set the hmmalign wrapper for this validator."""
        self.hmmalign = hmmalign

    def set_hmm_profile_dir(self, hmm_profile_dir: Path) -> None:
        """Set the directory containing HMM profiles for this validator."""
        self.hmm_profile_dir = hmm_profile_dir

    def set_marker(self, marker: Marker) -> None:
        """Set the marker type for this validator."""
        self.marker = marker

    def set_taxonomy_resolver(self, resolver: TaxonResolver) -> None:
        """Set the taxonomy resolver for this validator."""
        self.taxonomy_resolver = resolver

    def set_idservice(self, idservice: IDService) -> None:
        """Set the idservice for this validator."""
        self.idservice = idservice

    @staticmethod
    def requires_marker() -> bool:
        """Override to declare if validator needs a marker."""
        return False

    @staticmethod
    def requires_resolver() -> bool:
        """Override to declare if validator needs a taxonomy resolver."""
        return False

    @staticmethod
    def requires_idservice() -> bool:
        """Override to declare if validator needs an idservice."""
        return False

    @staticmethod
    def requires_hmmalign() -> bool:
        """Override to declare if validator needs a hmmalign wrapper."""
        return False

    def validate(self, record: SeqRecord, result: DNAAnalysisResult) -> None:
        """Validate a sequence record. Must be implemented by subclasses."""
        raise NotImplementedError