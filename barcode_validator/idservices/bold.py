from typing import Set, Dict, List, Tuple
from Bio.SeqRecord import SeqRecord

from nbitk.Taxon import Taxon
from nbitk.config import Config
from nbitk.Services.BOLD.IDService import IDService as BOLDIDService
from barcode_validator.constants import TaxonomicRank
from .idservice import IDService
from ..dna_analysis_result import DNAAnalysisResult


class BOLD(IDService):
    """
    BOLD identification service implementation.

    This service uses the BOLD Systems identification engine to identify
    taxonomic information from sequence records. This is a self-contained
    implementation that doesn't require the id_engine module.
    """

    def __init__(self, config: Config):
        super().__init__(config)

    def _build_taxonomic_trees(self, results: List[Dict]) -> List[Taxon]:
        """
        Build taxonomic trees from BOLD results, merging taxa with the same name.

        :param results: A list of dictionaries containing BOLD results
        :return: A list of root Taxon objects representing the taxonomic trees
        """
        # Dictionary to store unique taxa by (name, rank) pairs
        taxa_registry = {}

        # Define taxonomic levels in order from highest to lowest
        levels = [
            ("phylum", TaxonomicRank.PHYLUM),
            ("class", TaxonomicRank.CLASS),
            ("order", TaxonomicRank.ORDER),
            ("family", TaxonomicRank.FAMILY),
            ("subfamily", TaxonomicRank.SUBFAMILY),
            ("genus", TaxonomicRank.GENUS),
            ("species", TaxonomicRank.SPECIES)
        ]

        # Track root taxa (highest level taxa with no parents)
        root_taxa = set()

        for result in results[0]["results"]:
            previous_taxon = None

            # Build taxonomy chain from highest to lowest level
            for level_name, rank in levels:
                taxon_name = result.get(level_name)
                if not taxon_name or not taxon_name.strip():
                    continue

                taxon_name = taxon_name.strip()
                taxon_key = (taxon_name, rank.value)

                # Get or create taxon
                if taxon_key not in taxa_registry:
                    confidence = result.get("pct_identity", 0.0) / 100.0
                    taxon = Taxon(
                        name=taxon_name,
                        taxonomic_rank=rank.value,
                        confidence=confidence
                    )
                    taxa_registry[taxon_key] = taxon

                    # If this is the first (highest) level taxon, it's a root
                    if previous_taxon is None:
                        root_taxa.add(taxon)
                else:
                    taxon = taxa_registry[taxon_key]
                    # Update confidence if this result has higher confidence
                    current_confidence = result.get("pct_identity", 0.0) / 100.0
                    if current_confidence > taxon.confidence:
                        taxon.confidence = current_confidence

                # Link to parent if exists
                if previous_taxon is not None:
                    # Check if this taxon is already a child of the previous taxon
                    if taxon not in previous_taxon.clades:
                        previous_taxon.clades.append(taxon)

                previous_taxon = taxon

        return list(root_taxa)

    def _extract_taxa_at_level(self, trees: List[Taxon], level: TaxonomicRank) -> Set[Taxon]:
        """
        Extract all taxa at the specified taxonomic level from trees.

        :param trees: A list of root Taxon objects representing taxonomic trees
        :param level: The taxonomic rank at which to extract taxa
        :return: A set of Taxon objects at the specified level
        """
        taxa_at_level = set()

        def traverse_tree(taxon: Taxon, target_level: str):
            if taxon.taxonomic_rank == target_level:
                taxa_at_level.add(taxon)

            # Continue traversing children
            for child in taxon.clades:
                traverse_tree(child, target_level)

        for tree in trees:
            traverse_tree(tree, level.value)

        return taxa_at_level

    def identify_batch(self, batch: List[Tuple[DNAAnalysisResult, SeqRecord, Taxon]]) -> None:
        """
        Identify the taxonomic classification of a batch of sequence records.
        :param batch: List of tuples containing (DNAAnalysisResult, SeqRecord, Taxon)
        """
        try:
            self.logger.info(f"Identifying batch of {len(batch)} sequences using BOLD")

            # Initialize BOLD ID service with configuration
            bold_config = Config()
            bold_config.config_data = {
                'bold_database': 1,
                'bold_timeout': 300,
                'bold_params': {
                    'mi': self.min_identity,
                    'maxh': self.max_target_seqs
                }
            }
            bold_config.initialized = True
            service = BOLDIDService(bold_config)

            # Extract records from batch
            records = [item[1] for item in batch]

            # Wait for processing and get results
            results = service.identify_seqrecords(records)
            if not results:
                self.logger.warning("No BOLD results found for batch")
                return

            # Post-process results to build taxonomic trees and extract taxa
            for item in batch:
                result = item[0]
                record = item[1]

                # Filter results for the current record
                subset = [r for r in results if r['seqid'] == record.id]

                # Build taxonomic trees from the subset of results
                trees = self._build_taxonomic_trees(subset)

                # Extract taxa at requested level for each record in the batch
                taxa_at_level = self._extract_taxa_at_level(trees, TaxonomicRank(result.level))
                result.taxa = taxa_at_level

        except Exception as e:
            self.logger.error(f"Error identifying batch with BOLD: {e}")
            for item in batch:
                item[0].error = str(e)

    @staticmethod
    def requires_resolver():
        """BOLD service does not require a taxonomy resolver."""
        return False

    @staticmethod
    def requires_blastn():
        """BOLD service does not require BLAST."""
        return False