import tarfile
from pathlib import Path
from typing import Optional
from Bio.SeqRecord import SeqRecord
from nbitk.Taxon import Taxon
from nbitk.config import Config
from nbitk.Phylo.NCBITaxdmp import Parser
from tests.test_taxon import taxon

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

        # There seem to be two options here. Either we have a sequence record that is parsed from a GenBank file or a
        # BCDM-style file that has a taxon_id field, where in either case we have a taxon_id that is a string of digits,
        # or we have a string of letters. If it's a string of letters, it's just a taxon name.
        # These two scenarios mirror what happens in parse_id().
        if id_string.isdigit():
            return int(node.guids['taxon']) == int(id_string)
        else:
            return node.name == id_string

    def parse_id(self, record: SeqRecord) -> Optional[str]:
        """
        Assumes the sequence was parsed from a GenBank file.
        :param record: A Bio.SeqRecord object.
        :return: Taxon ID
        """
        # Assuming you parsed a GenBank format file (not FASTA)
        taxon_id = None
        for feature in record.features:
            if feature.type == "source":
                for xref in feature.qualifiers.get("db_xref", []):
                    if xref.startswith("taxon:"):
                        taxon_id = xref.split(":")[1]  # e.g., "taxon:9606" -> "9606"

        # Assuming you parsed a BCDM format file with a taxon_id field for NCBI
        if taxon_id is None:
            self.logger.warning(f"No taxon ID found in record {record.id}.")
            taxon_id = record.annotations.get('bcdm_fields', {}).get('taxon_id', None)

        # If still not found, try to get it from the identification field, being the name string fallback
        if taxon_id is None:
            self.logger.warning(f"No taxon ID found in record {record.id}.")
            taxon_id = record.annotations.get('bcdm_fields', {}).get('identification', None)

        return taxon_id

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