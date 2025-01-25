import logging
from enum import Enum
from Bio import Entrez
from typing import Optional, Dict, List
import time

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


def get_translation_table(marker: Marker, taxonomy_dict: dict) -> int:
    """
    Determine the appropriate translation table based on marker and taxonomy.

    Args:
        marker: The genetic marker being analyzed
        taxonomy_dict: Dictionary containing taxonomic classification with keys:
            'phylum', 'class', 'family'

    Returns:
        int: Translation table index for use with Biopython
    """

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


class TaxonomyResolver:

    """
    SYNOPSIS
    # Example usage:
    >>> resolver = TaxonomyResolver("your.email@example.com")
    >>>
    >>> # Example 1: Get full taxonomy for a species
    >>> taxonomy = resolver.get_taxonomy_dict("Homo sapiens", kingdom="Metazoa")
    >>> print("Full taxonomy:", taxonomy)
    >>>
    >>> # Example 2: Get only specific ranks needed for genetic code determination
    >>> required_ranks = ['phylum', 'class', 'family']
    >>> specific_taxonomy = resolver.get_lineage_at_ranks("Homo sapiens", required_ranks)
    >>>print("Specific ranks:", specific_taxonomy)
    >>>
    >>> # Example 3: Resolve a genus name
    >>> genus_taxonomy = resolver.get_taxonomy_dict("Drosophila", kingdom="Metazoa")
    >>> print("Genus taxonomy:", genus_taxonomy)
    """

    def __init__(self, email: str, logger: logging.Logger):
        """
        Initialize the taxonomy resolver.
        :param email: Email address for NCBI Entrez queries (required by NCBI)

        """
        Entrez.email = email
        self.logger = logger
        self.tax_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    def get_taxonomy_dict(self, taxon_name: str, kingdom: Optional[str] = None) -> Optional[Dict[str, str]]:
        """
        Resolve a taxon name to its full taxonomy using NCBI taxonomy database.

        Args:
            taxon_name: The taxonomic name to resolve
            kingdom: Optional kingdom name to filter results

        Returns:
            Dictionary containing taxonomic classification or None if not found
        """
        try:
            # First, search for the taxon ID
            search_term = f"{taxon_name}[Scientific Name]"
            if kingdom:
                search_term = f"{search_term} AND {kingdom}[Kingdom]"
            self.logger.info(f"Searching for '{search_term}' at Entrez")

            # Search in taxonomy database
            handle = Entrez.esearch(db="taxonomy", term=search_term)
            record = Entrez.read(handle)
            handle.close()
            self.logger.debug(f"Search results: {record}")

            if not record['IdList']:
                return None

            # Get the first (most relevant) taxonomy ID
            taxid = record['IdList'][0]

            # Fetch the full taxonomy record
            # TODO: Use the memory-resident version of the taxonomy database
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            if not records:
                return None

            record = records[0]

            # Extract lineage information
            lineage_ex = record.get('LineageEx', [])
            taxonomy_dict = {}

            # Build dictionary of rank -> name
            for entry in lineage_ex:
                rank = entry.get('Rank', '').lower()
                if rank in self.tax_ranks:
                    taxonomy_dict[rank] = entry.get('ScientificName', '')

            # Add the query taxon's rank and name
            query_rank = record.get('Rank', '').lower()
            if query_rank in self.tax_ranks:
                taxonomy_dict[query_rank] = record.get('ScientificName', '')

            return taxonomy_dict

        except Exception as e:
            print(f"Error resolving taxonomy: {str(e)}")
            return None

    def get_lineage_at_ranks(self, taxon_name: str, required_ranks: List[str],
                             kingdom: Optional[str] = None) -> Optional[Dict[str, str]]:
        """
        Get taxonomy information for specific required ranks.

        Args:
            taxon_name: The taxonomic name to resolve
            required_ranks: List of taxonomic ranks required
            kingdom: Optional kingdom name to filter results

        Returns:
            Dictionary containing only the required ranks or None if not found
        """
        full_taxonomy = self.get_taxonomy_dict(taxon_name, kingdom)
        if not full_taxonomy:
            return None

        return {rank: full_taxonomy.get(rank, '') for rank in required_ranks}


