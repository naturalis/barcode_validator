from nbitk.config import Config

from barcode_validator.validators.non_coding import NonCodingValidator
from barcode_validator.validators.protein_coding import ProteinCodingValidator
from barcode_validator.validators.structural import StructuralValidator
from barcode_validator.constants import Marker


class StructureValidatorFactory:

    @staticmethod
    def create_validator(config: Config, marker: Marker) -> StructuralValidator:
        """
        Factory method to create a structural validator for a given marker.
        :param config: Configuration object containing validation parameters
        :param marker: Marker enum
        :return: A subclass of StructuralValidator
        """
        if marker in [Marker.COI_5P, Marker.MATK, Marker.RBCL]:
            instance = ProteinCodingValidator(config)
        else:
            instance = NonCodingValidator(config)
        if instance.requires_marker():
            instance.set_marker(marker)
        return instance