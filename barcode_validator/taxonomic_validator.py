from typing import Dict, Optional, Tuple
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.BaseTree import Tree
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from nbitk.Taxon import Taxon
from barcode_validator.blast_runner import BlastRunner
from barcode_validator.taxonomy_resolver import TaxonomyResolver
from enum import Enum


class TaxonomicBackbone(Enum):
    DWC = "dwc"  # DarwinCore Archive format (including NSR)
    BOLD = "bold"  # BOLD taxonomy from Excel


class TaxonomicValidator:
    """
    Validates taxonomic assignments by comparing expected taxa against BLAST results.

    Examples:
        >>> from nbitk.config import Config
        >>> config = Config()
        >>> config.load_config('/path/to/config.yaml')
        >>> validator = TaxonomicValidator(config, ncbi_tree, backbone_tree)
        >>> is_valid, details = validator.validate_taxonomy(record)

    :param config: Configuration object containing validation parameters
    :param ncbi_tree: NCBI taxonomy tree for reference lookups
    :param backbone_tree: Taxonomy tree from chosen backbone source
    :param taxonomy_resolver: Optional TaxonomyResolver instance (created if not provided)
    """

    def __init__(self, config: Config, ncbi_tree: Tree, backbone_tree: Tree,
                 taxonomy_resolver: Optional[TaxonomyResolver] = None):
        self.config = config
        self.logger = get_formatted_logger(self.__class__.__name__, config)
        self.ncbi_tree = ncbi_tree
        self.backbone_tree = backbone_tree
        self.backbone_type = TaxonomicBackbone(config.get('taxonomic_backbone', 'bold'))
        self.validation_rank = config.get('validation_rank', 'family')
        self.constraint_rank = config.get('constraint_rank', 'class')
        
        # Initialize taxonomy resolver if not provided
        self.taxonomy_resolver = taxonomy_resolver or TaxonomyResolver(
            config.get('entrez_email', 'bioinformatics@naturalis.nl'),
            self.logger,
            self.ncbi_tree
        )

        # Initialize BLAST runner
        self.blast_runner = BlastRunner(config)
        self.blast_runner.ncbi_tree = self.ncbi_tree

    def validate_taxonomy(self, record: SeqRecord, expected_taxon: Optional[str] = None) -> Tuple[bool, Dict]:
        """
        Validate the taxonomic assignment of a sequence.

        :param record: The DNA sequence record to validate
        :param expected_taxon: Optional explicit taxon name (for tabular input)
        :return: Tuple of (validation_success, validation_details)
        """
        # Get expected taxon
        resolved_taxon = self._resolve_expected_taxon(record, expected_taxon)
        if not resolved_taxon:
            return False, {'error': 'Could not resolve expected taxon'}

        # Get validation rank taxon
        rank_taxon = self._get_rank_taxon(resolved_taxon)
        if not rank_taxon:
            return False, {'error': f'Could not find {self.validation_rank} rank taxon'}

        # Get constraint rank taxon for BLAST search
        constraint_taxon = self._get_constraint_taxon(resolved_taxon)
        if not constraint_taxon:
            self.logger.warning(f"Using fallback constraint taxon Eukaryota (taxon:2759)")
            constraint_id = "2759"
        else:
            constraint_id = constraint_taxon.guids.get('taxon')
            if not constraint_id:
                self.logger.warning(f"No NCBI taxon ID for constraint taxon {constraint_taxon}")
                constraint_id = "2759"

        # Get observed taxa from BLAST
        observed_taxa = self.blast_runner.run_localblast(
            record, 
            constraint_id,
            self.validation_rank
        )
        if not observed_taxa:
            return False, {'error': 'BLAST search failed'}

        # Check if expected taxon is among observed
        is_valid = rank_taxon in observed_taxa
        
        details = {
            'expected_taxon': rank_taxon,
            'validation_rank': self.validation_rank,
            'constraint_rank': self.constraint_rank,
            'constraint_taxon': constraint_taxon,
            'observed_taxa': observed_taxa,
            'backbone_source': self.backbone_type.value
        }
        
        return is_valid, details

    def _resolve_expected_taxon(self, record: SeqRecord, explicit_taxon: Optional[str] = None) -> Optional[Taxon]:
        """
        Resolve the expected taxon from either record or explicit name.

        :param record: The DNA sequence record
        :param explicit_taxon: Optional explicit taxon name
        :return: Resolved taxon object or None if resolution fails
        """
        if explicit_taxon:
            return self._resolve_from_name(explicit_taxon)
            
        if self.backbone_type == TaxonomicBackbone.BOLD:
            return self._resolve_from_bold(record)
        else:  # DarwinCore (including NSR)
            return self._resolve_from_dwc(record)

    def _resolve_from_name(self, taxon_name: str) -> Optional[Taxon]:
        """
        Resolve a taxon name from the backbone taxonomy.

        :param taxon_name: Name of the taxon to resolve
        :return: Resolved taxon object or None if not found
        """
        for node in self.backbone_tree.find_clades():
            if node.name.lower() == taxon_name.lower():
                return node
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

        # Find tip in BOLD tree by process ID
        for node in self.backbone_tree.get_terminals():
            if process_id in node.guids:
                return node
        return None

    def _resolve_from_dwc(self, record: SeqRecord) -> Optional[Taxon]:
        """
        Resolve taxon from DarwinCore fields in record.

        :param record: The DNA sequence record
        :return: Resolved taxon object or None if not found
        """
        taxonomy = record.annotations.get('taxonomy', [])
        if not taxonomy:
            self.logger.error(f"No taxonomy found in record {record.id}")
            return None

        # Try to find most specific taxon in backbone tree
        for name in reversed(taxonomy):  # Start from most specific
            if not name:
                continue
            taxon = self._resolve_from_name(name)
            if taxon:
                return taxon
        return None

    def _get_rank_taxon(self, taxon: Taxon) -> Optional[Taxon]:
        """
        Get the taxon at the validation rank from a more specific taxon.

        :param taxon: The specific taxon to start from
        :return: Taxon at validation rank or None if not found
        """
        for node in self.backbone_tree.root.get_path(taxon):
            if str(node.taxonomic_rank).lower() == self.validation_rank.lower():
                return node
        return None

    def _get_constraint_taxon(self, taxon: Taxon) -> Optional[Taxon]:
        """
        Get the taxon at the constraint rank from a more specific taxon.

        :param taxon: The specific taxon to start from
        :return: Taxon at constraint rank or None if not found
        """
        for node in self.backbone_tree.root.get_path(taxon):
            if str(node.taxonomic_rank).lower() == self.constraint_rank.lower():
                return node
        return None