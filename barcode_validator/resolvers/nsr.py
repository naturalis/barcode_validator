from Bio.SeqRecord import SeqRecord
from nbitk.Taxon import Taxon
from nbitk.config import Config
from .taxonomy import TaxonResolver
from barcode_validator.constants import TaxonomicBackbone
from nbitk.Phylo.DwCATaxonomyIO import Parser

class NSRResolver(TaxonResolver):
    def __init__(self, config: Config):
        super().__init__(config)
        self.parser = Parser(None)

    def get_type(self) -> TaxonomicBackbone:
        return TaxonomicBackbone.NSR

    def check_id(self, id_string: str, node: Taxon) -> bool:
        """
        Does a simple string comparison.
        :param id_string: A name string.
        :param node: Focal node.
        :return: True if the name matches literally
        """
        return node.name == id_string

    def parse_id(self, record: SeqRecord) -> str:
        """
        Parses the identification field from the bcdm annotations
        provided in the CSC/BCDM format.
        :param record: A Bio.SeqRecord object.
        :return: The identification field, i.e. a literal taxon name
        """
        return record.annotations['bcdm_fields']['identification']