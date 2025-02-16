from typing import Optional
from Bio.SeqRecord import SeqRecord
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from barcode_validator.blast_runner import BlastRunner
from barcode_validator.taxonomy_resolver import TaxonomyResolver
from barcode_validator.dna_analysis_result import DNAAnalysisResult
from enum import Enum


class TaxonomicBackbone(Enum):
    DWC = "dwc"  # DarwinCore Archive format (including NSR)
    BOLD = "bold"  # BOLD taxonomy from Excel


class TaxonomicValidator:
    """
    Validator of DNA barcodes via BLAST-based reverse taxonomy.

    Validates taxonomic assignments by comparing expected taxa against those observed in
    BLAST results. Uses a backbone-first approach where identifications are first resolved
    against a configured backbone taxonomy before mapping to NCBI for validation.

    The class exposes a single method for validation (`validate_taxonomy`). This method
    is invoked by the overall orchestrator as part of its composed validation steps.
    In turn, this class delegates some of its work to the BLAST runner and TaxonomyResolver.
    The idea behind this is that the BLAST runner may change underlying implementation
    details (e.g. by switching databases, using a different search algorithm), and so may
    the TaxonomyResolver (e.g. by invoking web services).

    Examples:
        >>> from nbitk.config import Config
        >>> config = Config()
        >>> config.load_config('/path/to/config.yaml')
        >>> validator = TaxonomicValidator(config, ncbi_tree, backbone_tree)
        >>> result = DNAAnalysisResult("sequence_id")
        >>> validator.validate_taxonomy(record, result)

    :param config: Configuration object containing validation parameters
    :param taxonomy_resolver: TaxonomyResolver instance
    """

    def __init__(self, config: Config, taxonomy_resolver: TaxonomyResolver):
        self.config = config
        self.logger = get_formatted_logger(self.__class__.__name__, config)
        self.backbone_type = TaxonomicBackbone(config.get('taxonomic_backbone', 'bold'))
        self.validation_rank = config.get('validation_rank', 'family')
        self.constraint_rank = config.get('constraint_rank', 'class')
        self.taxonomy_resolver = taxonomy_resolver

        # Initialize BLAST runner
        self.blast_runner = BlastRunner(config)
        self.blast_runner.ncbi_tree = taxonomy_resolver.ncbi_tree

    # TODO: why is expected_taxon an optional(!) string(!)
    def validate_taxonomy(self, record: SeqRecord, result: DNAAnalysisResult,
                          expected_taxon: Optional[str] = None) -> None:
        """
        Validate the taxonomic assignment of a sequence using backbone-first approach.
        Populates the provided result object with validation outcomes.

        :param record: The DNA sequence record to validate
        :param result: DNAAnalysisResult object to populate with validation results
        :param expected_taxon: Optional explicit taxon name (for tabular input)
        """
        # Set the identification rank in result
        result.level = self.validation_rank

        # Get backbone taxon from TaxonomyResolver
        backbone_taxon = self.taxonomy_resolver.resolve_backbone(
            expected_taxon or record,
            self.backbone_type.value
        )
        if not backbone_taxon:
            result.error = 'Could not resolve expected taxon'
            return

        # Get taxa at validation rank in both taxonomies
        backbone_valid, ncbi_valid = self.taxonomy_resolver.get_validation_taxon(
            backbone_taxon,
            self.validation_rank
        )
        if not backbone_valid or not ncbi_valid:
            result.error = f'Could not find {self.validation_rank} rank taxon'
            return

        # Store expected taxon in result
        result.exp_taxon = backbone_valid

        # Get constraint taxon ID for BLAST
        constraint_id = self.taxonomy_resolver.get_constraint_taxon(
            ncbi_valid,
            self.constraint_rank
        )
        try:
            constraint_id_int = int(constraint_id)
        except ValueError:
            self.logger.warning(f"Invalid constraint taxon ID: {constraint_id}, using 2759")
            constraint_id_int = 2759

        # Run BLAST search
        observed_taxa = self.blast_runner.run_localblast(
            record,
            constraint_id_int,
            self.validation_rank
        )
        if not observed_taxa:
            result.error = 'BLAST search failed'
            return

        # Store observed taxa in result
        result.obs_taxon = observed_taxa

        # Add backbone source as ancillary data
        result.add_ancillary('backbone_source', self.backbone_type.value)