import logging
import re
from enum import Enum
from Bio import Entrez
from nbitk.Taxon import Taxon
from Bio.Phylo.BaseTree import Tree
from typing import Optional, Tuple

"""
SYNOPSIS
    >>> taxonomy = {
    >>>    'phylum': 'Chordata',
    >>>    'class': 'Mammalia',
    >>>    'family': 'Hominidae'
    >>> }
    >>> translation_table = get_translation_table(Marker.COI_5P, taxonomy)
"""

class Marker(Enum):
    COI_5P = "COI-5P"
    MATK = "matK"
    RBCL = "rbcL"

class TaxonomyResolver:
    """
    A class for resolving taxonomic names against a backbone taxonomy (NSR/BOLD/DwC)
    and mapping to NCBI taxonomy for validation.

    For each identification, first resolve against the configured backbone taxonomy,
    then map to NCBI for validation purposes. This ensures maximum coverage of local
    names while enabling validation against GenBank's reference data.

    Example usage:
        >>> resolver = TaxonomyResolver(email, logger, ncbi_tree, nsr_tree)
        >>> backbone_taxon = resolver.resolve_backbone("Homo sapiens")
        >>> ncbi_taxon = resolver.get_ncbi_taxon(backbone_taxon)
    """

    def __init__(self, email: str, logger: logging.Logger, ncbi_tree: Tree, backbone_tree: Tree):
        """
        Initialize the taxonomy resolver.
        :param email: Email address for NCBI Entrez queries
        :param logger: Logger instance
        :param ncbi_tree: NCBI taxonomy tree
        :param backbone_tree: Configured backbone taxonomy (NSR/BOLD/DwC)
        """
        Entrez.email = email
        self.logger = logger
        self.ncbi_tree = ncbi_tree
        self.backbone_tree = backbone_tree
        self.tax_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    def resolve_backbone(self, identification: str, rank: str = 'null') -> Optional[Taxon]:
        """
        Resolve a verbatim identification against the backbone taxonomy.

        :param identification: Verbatim taxonomic name
        :param rank: Provided rank of the identification
        :return: Resolved taxon in backbone taxonomy or None if not found
        """
        if not identification:
            return None

        # If no rank but looks like species, infer that
        rank = rank.lower()
        if rank == 'null':
            pattern = r'^[A-Z][a-z]+ [a-z]+$'
            if re.match(pattern, identification):
                self.logger.warning(f"Assuming {identification} is a species")

        # Search in backbone tree
        for node in self.backbone_tree.find_clades():
            if node.name.lower() == identification.lower():
                return node

        self.logger.warning(f"Taxon not found in backbone: {identification}")
        return None

    def get_ncbi_taxon(self, backbone_taxon: Taxon) -> Optional[Taxon]:
        """
        Map a backbone taxon to NCBI taxonomy.

        :param backbone_taxon: Taxon from backbone taxonomy
        :return: Corresponding NCBI taxon or None if not found
        """
        try:
            search_term = f"{backbone_taxon.name}[Scientific Name]"
            self.logger.info(f"Searching for '{search_term}' at Entrez")

            handle = Entrez.esearch(db="taxonomy", term=search_term)
            record = Entrez.read(handle)
            handle.close()

            if not record['IdList']:
                self.logger.warning(f"Taxon not found in NCBI: {backbone_taxon.name}")
                return None

            taxid = record['IdList'][0]
            for node in self.ncbi_tree.find_clades():
                if node.guids.get('taxon') == taxid:
                    return node

            return None

        except Exception as e:
            self.logger.error(f"Error resolving NCBI taxonomy: {str(e)}")
            return None

    def get_validation_taxon(self, backbone_taxon: Taxon, level: str) -> Tuple[Optional[Taxon], Optional[Taxon]]:
        """
        Get taxa at validation level in both backbone and NCBI taxonomies.

        :param backbone_taxon: Initially resolved taxon in backbone taxonomy
        :param level: Desired validation level
        :return: Tuple of (backbone_validation_taxon, ncbi_validation_taxon)
        """
        if not backbone_taxon:
            return None, None

        # First get the validation level taxon in the backbone
        backbone_validation = None
        for node in self.backbone_tree.root.get_path(backbone_taxon):
            if str(node.taxonomic_rank).lower() == level:
                backbone_validation = node
                break

        if not backbone_validation:
            return None, None

        # Then map to NCBI
        ncbi_validation = self.get_ncbi_taxon(backbone_validation)
        return backbone_validation, ncbi_validation

    def get_constraint_taxon(self, taxon: Taxon, constraint_level: str = 'class') -> str:
        """
        Get the NCBI taxon ID for BLAST constraint.

        :param taxon: The NCBI taxon to get constraint for
        :param constraint_level: Level at which to constrain (default: class)
        :return: NCBI taxon ID as string, defaults to Eukaryota (2759)
        """
        if not taxon:
            return '2759'

        for node in self.ncbi_tree.root.get_path(taxon):
            if str(node.taxonomic_rank).lower() == constraint_level:
                if 'taxon' in node.guids:
                    return node.guids['taxon']

        self.logger.warning("Using Eukaryota (taxon:2759) as BLAST constraint")
        return '2759'

    def get_translation_table(self, marker: Marker, taxon: Taxon) -> int:
        """
        Determine the appropriate translation table based on marker and taxonomy.

        :param marker: The genetic Marker (enum) being analyzed
        :param taxon: The NCBI taxon object representing the tip of the classification
        :return: Translation table index (int) for use with Biopython
        """
        taxonomy_dict = {}
        for node in self.ncbi_tree.root.get_path(taxon):
            if node.taxonomic_rank in self.tax_ranks:
                taxonomy_dict[node.taxonomic_rank] = node.name

        if marker in [Marker.MATK, Marker.RBCL]:
            if taxonomy_dict.get('family') == 'Balanophoraceae':
                return 32
            return 11

        elif marker == Marker.COI_5P:
            phylum = taxonomy_dict.get('phylum')
            tax_class = taxonomy_dict.get('class')
            family = taxonomy_dict.get('family')

            if phylum == 'Chordata':
                if tax_class == 'Ascidiacea':
                    return 13
                elif tax_class in ['Actinopteri', 'Amphibia', 'Mammalia', 'Aves', 'Reptilia']:
                    return 2

            elif phylum == 'Hemichordata':
                if family == 'Cephalodiscidae':
                    return 33
                elif family == 'Rhabdopleuridae':
                    return 24

            elif phylum in ['Echinodermata', 'Platyhelminthes']:
                return 9

            # Default invertebrate mitochondrial code for other invertebrates
            if phylum != 'Chordata':
                return 5

            # Fall back to standard code if no specific rules match
            return 1