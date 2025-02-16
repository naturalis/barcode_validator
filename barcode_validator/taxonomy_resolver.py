import sys
import tarfile
from enum import Enum
from pathlib import Path

import requests
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from nbitk.Taxon import Taxon
from nbitk.Phylo.NCBITaxdmp import Parser as NCBIParser
from nbitk.Phylo.DwCATaxonomyIO import Parser as DwCParser
from nbitk.Phylo.BOLDXLSXIO import Parser as BOLDParser
from typing import Optional, Tuple

from nbitk.config import Config
from nbitk.logger import get_formatted_logger

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
    Central resolver for taxonomic names against different backbone taxonomies.

    The resolver contributes to the overall validation process in the steps where
    the provided identifications of DNA barcodes are mapped onto a configured taxonomy
    (currently either the NSR taxonomy or the BOLD taxonomy).

    The overall procedure then advances from those initial names (often species) up to
    a configured higher level (e.g., family). The taxon at this level is then reconciled
    with the equivalent taxon in reference database against which the DNA barcodes are
    validated (i.e. GenBank, which uses the NCBI taxonomy).

    The resolver is the only place where taxonomic names are mapped to taxon IDs. Hence,
    it is the only place where various taxonomy data types are processed, such as
    NCBI taxdump format, DarwinCore archive format, and BOLD submissio spreadsheets. It
    is also the only place where taxonomic topologies are traversed and where taxonomic
    reconciliation is performed.

    Note that the resolver does not perform any taxonomic validation itself: it merely
    provides the necessary tooling for accessing taxonomic data. The taxonomic validation
    is otherwise orchestrated by the TaxonomicValidator class.

    Note also that the indexes of translation tables as used in bioinformatics (in this
    case to do amino acid translation) are decorations attached to the NCBI taxonomy. Hence,
    thes table indexes are also provided by the resolver.

    Examples:
        >>> resolver = TaxonomyResolver(config)
        >>> backbone_taxon = resolver.resolve_backbone("Homo sapiens", "bold")
        >>> backbone_valid, ncbi_valid = resolver.get_validation_taxon(backbone_taxon, "family")
    """


    def __init__(self, config: Config):
        """
        Initialize the taxonomy resolver.
        :param config: NBITK configuration object
        """
        Entrez.email = config.get('entrez_email', 'anonymous@example.org')
        self.config = config
        self.logger = get_formatted_logger(self.__class__.__name__, config)
        self.ncbi_tree = None
        self.backbone_tree = None
        self.tax_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    def download_ncbi_dump(self, destination: str) -> str:
        """
        Download NCBI taxdump.

        :param destination: Location to store the downloaded file.
        :return: Location of the downloaded file.
        """
        url = self.config.get("ncbi_taxonomy_url")
        if not url:
            self.logger.error("NCBI taxonomy URL is not set in the configuration.")
            sys.exit(1)

        self.logger.info(f"Downloading NCBI taxdump from {url}...")

        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()  # Raise an error for bad HTTP responses (4xx or 5xx)

            with open(destination, "wb") as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)

            self.logger.info(f"Download completed: {destination}")
            return destination

        except requests.exceptions.RequestException as e:
            self.logger.error(f"Failed to download NCBI taxdump: {e}")
            sys.exit(1)

    def download_nsr_dump(self, destination: str) -> str:
        """
        Download NSR dump.

        :param destination: Location to store the downloaded file.
        :return: Location of the downloaded file.
        """
        url = self.config.get("dwc_archive_url")
        if not url:
            self.logger.error("NSR taxonomy URL is not set in the configuration.")
            sys.exit(1)

        self.logger.info(f"Downloading NSR dump from {url}...")

        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()  # Raise an error for bad HTTP responses (4xx or 5xx)

            with open(destination, "wb") as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)

            self.logger.info(f"Download completed: {destination}")
            return destination

        except requests.exceptions.RequestException as e:
            self.logger.error(f"Failed to download NSR dump: {e}")
            sys.exit(1)

    def setup_taxonomy(self) -> None:
        """
        Setup taxonomy dumps using the command line arguments and update configuration.
        """
        # Check if the NCBI taxonomy is there. If not, download it. Then parse it.
        ncbi_path = self.config.get("ncbi_taxonomy")
        if ncbi_path is not None:
            path_obj = Path(ncbi_path)
            if not path_obj.exists():
                downloaded_path = self.download_ncbi_dump(ncbi_path)
                self.config.set("ncbi_taxonomy", downloaded_path)
            self.ncbi_tree = NCBIParser(tarfile.open(self.config.get("ncbi_taxonomy"), "r:gz")).parse()
        else:
            self.logger.error("NCBI taxonomy is not configured. Skipping setup.")
            sys.exit(1)

        # Check if the NSR taxonomy is there. If not, download it. Then parse it.
        nsr_path = self.config.get("dwc_archive")
        if nsr_path is not None:
            path_obj = Path(nsr_path)
            if not path_obj.exists():
                downloaded_path = self.download_nsr_dump(nsr_path)
                self.config.set("dwc_archive", downloaded_path)
            self.backbone_tree = self.backbone_tree = DwCParser(self.config.get("dwc_archive")).parse()

        # No NSR tree configured, assuming BOLD
        else:
            bold_file = self.config.get('bold_sheet_file')
            with open(bold_file, 'rb') as f:
                self.backbone_tree = BOLDParser(f).parse()

    def resolve_backbone(self, identification: str, backbone_type: str = "bold") -> Optional[Taxon]:
        """
        Resolve an identification against the backbone taxonomy.

        :param identification: Taxonomic name or record to resolve
        :param backbone_type: Type of backbone taxonomy ("bold" or "dwc")
        :return: Resolved taxon in backbone taxonomy or None if not found
        """
        if isinstance(identification, SeqRecord):
            if backbone_type == "bold":
                return self._resolve_from_bold(identification)
            else:  # dwc
                return self._resolve_from_dwc(identification)
        else:
            return self._resolve_from_name(identification)

    def _resolve_from_name(self, taxon_name: str) -> Optional[Taxon]:
        """
        Resolve a taxon name from the backbone taxonomy.

        :param taxon_name: Name of the taxon to resolve
        :return: Resolved taxon object or None if not found
        """
        if not taxon_name:
            return None

        # Look for exact match first
        for node in self.backbone_tree.find_clades():
            if node.name and node.name.lower() == taxon_name.lower():
                return node

        # If no exact match, try case-insensitive partial match
        matches = []
        for node in self.backbone_tree.find_clades():
            if node.name and taxon_name.lower() in node.name.lower():
                matches.append(node)

        if len(matches) == 1:
            return matches[0]
        elif len(matches) > 1:
            self.logger.warning(f"Multiple matches found for {taxon_name}")
            return None

        self.logger.warning(f"No matches found for {taxon_name}")
        return None

    def _resolve_from_bold(self, record: SeqRecord) -> Optional[Taxon]:
        """
        Resolve taxon from BOLD process ID in record.

        :param record: The DNA sequence record
        :return: Resolved taxon object or None if not found
        """
        process_id = record.annotations.get('bcdm_fields', {}).get('processid')
        if not process_id:
            self.logger.error(f"No process ID found in record {record.id}")
            return None

        # Find tip in backbone tree by process ID
        for node in self.backbone_tree.get_terminals():
            if process_id in node.guids:
                return node

        # If no match by process ID, try identification field
        if 'identification' in record.annotations.get('bcdm_fields', {}):
            return self._resolve_from_name(record.annotations['bcdm_fields']['identification'])

        return None

    def _resolve_from_dwc(self, record: SeqRecord) -> Optional[Taxon]:
        """
        Resolve taxon from DarwinCore fields in record.

        :param record: The DNA sequence record
        :return: Resolved taxon object or None if not found
        """
        # Try full taxonomy first
        taxonomy = record.annotations.get('taxonomy', [])
        if taxonomy:
            # Start from most specific (last) to least specific
            for name in reversed(taxonomy):
                if not name:
                    continue
                taxon = self._resolve_from_name(name)
                if taxon:
                    return taxon

        # Try description field if taxonomy failed
        if record.description:
            return self._resolve_from_name(record.description)

        self.logger.error(f"No taxonomy found in record {record.id}")
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
        :param level: Desired validation level (e.g., "family")
        :return: Tuple of (backbone_validation_taxon, ncbi_validation_taxon)
        """
        if not backbone_taxon:
            return None, None

        # Get validation level taxon in backbone
        backbone_valid = None
        for node in self.backbone_tree.root.get_path(backbone_taxon):
            if str(node.taxonomic_rank).lower() == level.lower():
                backbone_valid = node
                break

        if not backbone_valid:
            return None, None

        # Map to NCBI
        ncbi_valid = self.get_ncbi_taxon(backbone_valid)
        return backbone_valid, ncbi_valid

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
            if str(node.taxonomic_rank).lower() == constraint_level.lower():
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
