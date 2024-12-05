from enum import Enum

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


