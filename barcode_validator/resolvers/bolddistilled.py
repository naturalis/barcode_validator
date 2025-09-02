from Bio.SeqRecord import SeqRecord
from nbitk.Taxon import Taxon
from nbitk.config import Config
from nbitk.Phylo.BOLDDistilled import Parser
from barcode_validator.resolvers.taxonomy import TaxonResolver
from barcode_validator.constants import TaxonomicBackbone


class BoldDistilledResolver(TaxonResolver):
    def __init__(self, config: Config):
        super().__init__(config)
        self.parser = Parser(None)

    def get_type(self) -> TaxonomicBackbone:
        return TaxonomicBackbone.BOLDDISTILLED

    def check_id(self, id_string: str, node: Taxon) -> bool:
        """
        Returns True if the focal node has the given ID as a key
        in the guids dictionary.
        :param id_string: A BOLD process ID.
        :param node: The focal node.
        :return: Whether the ID is found in the guids dictionary.
        """
        return id_string == node.name

    def parse_id(self, record: SeqRecord) -> str:
        """
        Fetches the BOLD process ID from the FASTA definition line
        :param record: A Bio.SeqRecord object.
        :return: A BOLD process ID.
        """
        return record.id.split("|")[1]