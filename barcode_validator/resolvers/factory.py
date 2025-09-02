from typing import Optional
from nbitk.config import Config
from barcode_validator.resolvers.bold import BoldResolver
from barcode_validator.resolvers.ncbi import NCBIResolver
from barcode_validator.resolvers.nsr import NSRResolver
from barcode_validator.resolvers.taxonomy import TaxonResolver
from barcode_validator.constants import TaxonomicBackbone
from barcode_validator.resolvers.bolddistilled import BoldDistilledResolver


class ResolverFactory:

    @staticmethod
    def create_resolver(config: Config, backbone_type: TaxonomicBackbone) -> Optional[TaxonResolver]:
        """
        Factory method to create a resolver based on the backbone taxonomy type.
        :param config: Config object.
        :param backbone_type: Enum value of TaxonomicBackbone.
        :return: TaxonResolver object or None if the backbone taxonomy type is not supported.
        """
        if backbone_type == TaxonomicBackbone.DWC:
            tr = NSRResolver(config)
        elif backbone_type == TaxonomicBackbone.BOLD:
            tr = BoldResolver(config)
        elif backbone_type == TaxonomicBackbone.NCBI:
            tr = NCBIResolver(config)
        elif backbone_type == TaxonomicBackbone.BOLDDISTILLED:
            tr = BoldDistilledResolver(config)
        else:
            return None
        return tr