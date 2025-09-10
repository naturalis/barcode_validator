from typing import Tuple, List

from Bio.SeqRecord import SeqRecord
from nbitk.Taxon import Taxon
from nbitk.config import Config
from barcode_validator.dna_analysis_result import DNAAnalysisResult, DNAAnalysisResultSet
from barcode_validator.validators.validator import AbstractValidator
from barcode_validator.constants import TaxonomicRank


class TaxonomicValidator(AbstractValidator):
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
        >>> validator = TaxonomicValidator(config)
        >>> result = DNAAnalysisResult("sequence_id")
        >>> validator.validate_taxonomy(SeqRecord(), DNAAnalysisResult())

    :param config: Configuration object containing validation parameters
    """

    def __init__(self, config: Config):
        super().__init__(config)
        self.config = config

    def validate(self, resultset: DNAAnalysisResultSet, is_struct_validated: bool = False) -> None:
        """
        Validate a batch of DNA sequence records taxonomically. Populates the provided
        result object with validation outcomes.
        :param resultset: Populated DNA analysis result object
        :param is_struct_validated: Whether structural validation has already been performed
        """
        # Set up constraint rank
        constraint_rank = TaxonomicRank.CLASS
        if self.config.get('local_blast.extent') is not None:
            constraint_rank = TaxonomicRank(self.config.get('local_blast.extent'))

        # Iterate over records and results, aggregate structurally valid ones into batches
        batch = []
        for result in resultset.results:
            record = result.seq_record

            # Add to batch if no problems so far in structural validation
            if is_struct_validated:
                if not result.error and result.check_seq_quality():
                    batch.append((result, record))
                else:
                    self.logger.warning(f"Skipping {record.id} for due to prior errors: {result.error}")

            # Otherwise, if no structural validation done yet, add all records without errors
            else:
                if not result.error:

                    # Make sure that result.exp_taxon is assigned through tr.enrich_result() because this may not
                    # have been done if we aren't doing structural validation
                    if result.exp_taxon is None:
                        self.taxonomy_resolver.enrich_result(record, result)
                        if result.error:
                            continue
                    batch.append((result, record))
                else:
                    self.logger.warning(f"Skipping {record.id} due to prior errors: {result.error}")

        # Process batches
        self.logger.info(f"Going to taxonomically validate {len(batch)} records")
        for i in range(0, len(batch), 100):
            batch_slice = batch[i:i + 100]
            self.logger.info(f"Taxonomically validating batch {i // 100 + 1} containing {len(batch_slice)} records")
            self.validate_batch(batch_slice, constraint_rank)


    def validate_batch(self, batch: List[Tuple[DNAAnalysisResult,SeqRecord]], constraint: TaxonomicRank = TaxonomicRank.CLASS):
        """
        Validates taxonomic assignments by comparing expected taxa against observed in BLAST results.
        :param batch: A list of tuples of result objects, sequence records, and query extents
        """
        batch = self._annotate_batch(batch, constraint)
        self.idservice.identify_batch(batch)

    def _annotate_batch(self, batch: List[Tuple[DNAAnalysisResult,SeqRecord]], constraint: TaxonomicRank = TaxonomicRank.CLASS) -> List[Tuple[DNAAnalysisResult,SeqRecord,Taxon]]:
        """
        Attempt to find constraint Taxon at provided TaxonomicRank for each input Tuple.
        :param batch: A list of tuples containing (DNAAnalysisResult, SeqRecord) objects
        :param constraint: A TaxonomicRank.CLASS
        :return: A list of tuples containing (DNAAnalysisResult, SeqRecord, Taxon) objects
        """
        annotated_batch = []
        for item in batch:
            result = item[0]
            record = item[1]
            constraint = None
            annotated_batch.append((result, record, constraint))

            # # Get taxon at validation rank in reflib taxonomy
            # reflib_matches = self.taxonomy_resolver.match_nodes(result.exp_taxon, {"name", "taxonomic_rank"})
            # if len(reflib_matches) == 0:
            #     result.error = f'Could not reconcile expected taxon {result.exp_taxon} at rank {result.level} with reflib taxonomy'
            # elif len(reflib_matches) > 1:
            #     result.error = f'Multiple taxa found in reflib taxonomy for expected taxon {result.exp_taxon} at rank {result.level}: {reflib_matches}'
            # reflib_valid = reflib_matches[0]
            #
            # # Get constraint taxon ID for BLAST
            # if result.error is None:
            #     reflib_constraint = self.taxonomy_resolver.find_ancestor_at_rank(reflib_valid, constraint)
            #     if reflib_constraint is None:
            #         result.error = f'Could not get constraint taxon ID for {reflib_valid} at rank {constraint}'
            #
            #     if result.error is None:
            #         annotated_batch.append((result, record, reflib_constraint))
        return annotated_batch

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
        :param constraint_rank: Taxonomic rank to use as constraint for BLAST search
        """

        # Get taxon at validation rank in reflib taxonomy
        reflib_matches = self.taxonomy_resolver.match_nodes(result.exp_taxon, {"name", "taxonomic_rank"})
        if len(reflib_matches) == 0:
            result.error = f'Could not reconcile expected taxon {result.exp_taxon} at rank {result.level} with reflib taxonomy'
            return
        elif len(reflib_matches) > 1:
            result.error = f'Multiple taxa found in reflib taxonomy for expected taxon {result.exp_taxon} at rank {result.level}: {reflib_matches}'
            return
        reflib_valid = reflib_matches[0]

        # Get constraint taxon ID for BLAST
        reflib_constraint = self.taxonomy_resolver.find_ancestor_at_rank(
            reflib_valid,
            constraint_rank)
        if reflib_constraint is None:
            result.error = f'Could not get constraint taxon ID for {reflib_valid} at rank {constraint_rank}'
            return

        # Run BLAST search
        observed_taxa = self.idservice.identify_record(
            record,
            TaxonomicRank(result.level.lower()),
            reflib_constraint
        )
        if not observed_taxa:
            result.error = f'ID service {self.idservice} search failed'
            return

        # Store observed taxa in result
        result.obs_taxon = observed_taxa

        # Add backbone source as ancillary data
        result.add_ancillary('backbone_source', self.taxonomy_resolver.get_type().value)

    @staticmethod
    def requires_resolver() -> bool:
        """Override to declare if validator needs a taxonomy resolver."""
        return True

    @staticmethod
    def requires_idservice() -> bool:
        """Override to declare if validator needs an idservice."""
        return True