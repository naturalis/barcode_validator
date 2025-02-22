from typing import Optional
from nbitk.config import Config
from barcode_validator.resolvers.bold import BoldResolver
from barcode_validator.resolvers.ncbi import NCBIResolver
from barcode_validator.resolvers.nsr import NSRResolver
from barcode_validator.resolvers.taxonomy import TaxonResolver
from barcode_validator.constants import TaxonomicBackbone


class ResolverFactory:

    @staticmethod
    def create_resolver(config: Config, backbone_type: TaxonomicBackbone) -> Optional[TaxonResolver]:
        """
        Factory method to create a resolver based on the backbone taxonomy type.
        :param config: Config object.
        :param backbone_type: Enum value of TaxonomicBackbone.
        :return: TaxonResolver object or None if the backbone taxonomy type is not supported.
        """
        if backbone_type.value == TaxonomicBackbone.NSR.value:
            tr = NSRResolver(config)
        elif backbone_type.value == TaxonomicBackbone.BOLD.value:
            tr = BoldResolver(config)
        elif backbone_type.value == TaxonomicBackbone.NCBI.value:
            tr = NCBIResolver(config)
        else:
            return None
        return tr