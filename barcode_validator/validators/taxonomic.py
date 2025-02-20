from typing import Optional
from Bio.SeqRecord import SeqRecord
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from barcode_validator.idservices.idservice import IDService
from barcode_validator.idservices.ncbi import NCBI
from barcode_validator.taxonomy_resolver import TaxonomicRank, TaxonomyResolver
from barcode_validator.dna_analysis_result import DNAAnalysisResult


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

    def __init__(self, config: Config, taxonomy_resolver: TaxonomyResolver, idservice: IDService):
        self.config = config
        self.logger = get_formatted_logger(self.__class__.__name__, config)
        self.taxonomy_resolver = taxonomy_resolver
        self.idservice = idservice

    def validate_taxonomy(self, record: SeqRecord, result: DNAAnalysisResult, constraint_rank: TaxonomicRank = TaxonomicRank.CLASS) -> None:
        """
        Validate the taxonomic assignment of a sequence using backbone-first approach.
        Populates the provided result object with validation outcomes. When this method
        is called, the result object should already have a result.exp_taxon field,
        which is the taxon at the validation rank that is expected for the record,
        and a result.level field, which is the taxonomic level at which the validation
        is being made.

        :param record: The DNA sequence record to validate
        :param result: DNAAnalysisResult object to populate with validation results
        """

        # Get taxon at validation rank in reflib taxonomy
        reflib_valid = self.taxonomy_resolver.find_taxon_at_level(
            result.exp_taxon,
            result.level
        )
        if reflib_valid is None:
            result.error = f'Could not reconcile expected taxon {result.exp_taxon} at rank {result.level} with reflib taxonomy'
            return

        # Get constraint taxon ID for BLAST
        reflib_constraint = self.taxonomy_resolver.get_constraint_taxon(
            reflib_valid,
            constraint_rank
        )
        if reflib_constraint is None:
            result.error = f'Could not get constraint taxon ID for {reflib_valid} at rank {self.constraint_rank}'
            return

        # Run BLAST search
        observed_taxa = self.idservice.identify_record(
            record,
            result.level,
            reflib_constraint
        )
        if not observed_taxa:
            result.error = 'BLAST search failed'
            return

        # Store observed taxa in result
        result.obs_taxon = observed_taxa

        # Add backbone source as ancillary data
        result.add_ancillary('backbone_source', self.taxonomy_resolver.backbone_type.value)