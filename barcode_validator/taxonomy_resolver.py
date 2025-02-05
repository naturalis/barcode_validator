import logging
from enum import Enum
from Bio import Entrez
from nbitk.Taxon import Taxon
from Bio.Phylo.BaseTree import Tree
from typing import Optional

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
    A class for resolving taxonomic names to full taxonomic classifications using NCBI taxonomy database.
    The class uses the Entrez module from Biopython to query the NCBI taxonomy database to obtain a taxon ID.
    The web service is aware of synonyms and alternate spellings, so it is somewhat robus. From that initial
    taxon ID, a nbitk.Taxon object is retrieved from the resident NCBI taxonomy tree.
    SYNOPSIS
    # Example usage:
    >>> resolver = TaxonomyResolver("your.email@example.com")
    >>> taxon = resolver.get_taxon("Homo sapiens")
    """

    def __init__(self, email: str, logger: logging.Logger, ncbi_tree: Tree):
        """
        Initialize the taxonomy resolver.
        :param email: Email address for NCBI Entrez queries (required by NCBI)
        """
        Entrez.email = email
        self.logger = logger
        self.ncbi_tree = ncbi_tree
        self.tax_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    def get_taxon(self, taxon_name: str) -> Optional[Taxon]:
        """
        Resolve a taxon name to its full taxonomy using NCBI taxonomy database.
        :param taxon_name: The scientific taxonomic name to resolve
        :return: A nbitk.Taxon object or None if not found
        """
        try:
            # First, search for the taxon ID
            search_term = f"{taxon_name}[Scientific Name]"
            self.logger.info(f"Searching for '{search_term}' at Entrez")

            # Search in taxonomy database
            handle = Entrez.esearch(db="taxonomy", term=search_term)
            record = Entrez.read(handle)
            handle.close()
            self.logger.debug(f"Search results: {record}")
            if not record['IdList']:
                self.logger.warning(f"Taxon not found: {taxon_name}")
                return None

            # Get the first (most relevant) taxonomy ID and the associated taxon from the resident tree
            taxid = record['IdList'][0]
            self.logger.info(f"Found taxid: {taxid}")
            for node in self.ncbi_tree.find_clades():
                if node.guids['taxon'] == taxid:
                    self.logger.debug(f"Found taxon: {node}")
                    return node
            self.logger.warning(f"Taxon not found: {taxon_name}")
            return None

        except Exception as e:
            print(f"Error resolving taxonomy: {str(e)}")
            return None

    def get_translation_table(self, marker: Marker, taxon: Taxon) -> int:
        """
        Determine the appropriate translation table based on marker and taxonomy.

        :param marker: The genetic Marker (enum) being analyzed
        :param taxon: The Taxon object representing the tip of the classification
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

