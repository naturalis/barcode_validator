import tarfile
from pathlib import Path
from typing import Optional
from Bio.SeqRecord import SeqRecord
from nbitk.Taxon import Taxon
from nbitk.config import Config
from nbitk.Phylo.NCBITaxdmp import Parser
from .taxonomy import TaxonResolver
from barcode_validator.constants import TaxonomicBackbone


class NCBIResolver(TaxonResolver):
    def __init__(self, config: Config):
        super().__init__(config)
        self.parser = Parser(None)

    def get_type(self) -> TaxonomicBackbone:
        return TaxonomicBackbone.NCBI

    def check_id(self, id_string: str, node: Taxon) -> bool:
        """
        Returns True if the focal node has the given ID under the
        taxon key in the guids dictionary.
        :param id_string: NCBI taxon ID to check.
        :param node: Focal node.
        :return:
        """
        return str(node.guids['taxon']) == id_string

    def parse_id(self, record: SeqRecord) -> Optional[str]:
        """
        Assumes the sequence was parsed from a GenBank file.
        :param record: A Bio.SeqRecord object.
        :return: Taxon ID
        """
        # Assuming you parsed a GenBank format file (not FASTA)
        for feature in record.features:
            if feature.type == "source":
                for xref in feature.qualifiers.get("db_xref", []):
                    if xref.startswith("taxon:"):
                        taxon_id = xref.split(":")[1]  # e.g., "taxon:9606" -> "9606"
                        return taxon_id
        return None

    def load_tree(self, file: Path) -> None:
        """
        Load a tree from a file.
        :param file: Location of the tree file.
        :return:
        """
        # Open the tar file
        tar = tarfile.open(file, "r:gz")
        self.parser.file = tar
        self.tree = self.parser.parse()