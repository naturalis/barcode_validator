from enum import Enum
from typing import Optional

from Bio import SeqRecord
from nbitk.config import Config
from nbitk.logger import get_formatted_logger

from .taxonomy_resolver import TaxonomyResolver, Marker, TaxonomicBackbone
from .dna_analysis_result import DNAAnalysisResult

class TaxonomyEnrichmentParser:
    """
    Base class for parsers that enrich result objects with taxonomic data.
    Works with TaxonomyResolver to map identifications to backbone.
    """

    def __init__(self, config: Config, taxonomy_resolver: TaxonomyResolver):
        self.config = config
        self.logger = None
        self.taxonomy: Optional[TaxonomicBackbone] = None
        self.taxonomy_resolver = taxonomy_resolver

    def get_query(self, record: SeqRecord) -> str:
        """
        Get whatever is going to be queried in the respective database.
        :return: A query string, e.g. a verbatim identification or a process ID
        """
        raise NotImplementedError

    def enrich_result(self, record: SeqRecord, result: DNAAnalysisResult, backbone: TaxonomicBackbone):
        """
        Enriches the result object with taxonomic data
        :param record: A SeqRecord object with process ID in the ID field.
        :param result: A DNAAnalysisResult object to populate with taxonomic data.
        :param backbone: The backbone taxonomy to use for the enrichment.
        :return:
        """
        query = self.get_query(record)

        # Use resolver to get backbone taxon
        species = self.taxonomy_resolver.resolve_backbone(query,backbone)
        if species is None:
            msg = f"No entry found for {query} in BOLD backbone taxonomy."
            self.logger.error(msg)
            result.error = msg
            return
        result.species = species

        # Get lineage from resolver's backbone tree
        for node in self.taxonomy_resolver.backbone_tree.root.get_path(result.species):
            if node.taxonomic_rank.lower() == result.level.lower():
                result.exp_taxon = node
                break

        # Check if there was a traversal problem
        if result.exp_taxon is None:
            msg = f"No {result.level} taxon found for {query} in BOLD backbone taxonomy."
            self.logger.error(msg)
            result.error = msg
            return

        # Infer the translation table
        ncbi_taxon = self.taxonomy_resolver.find_taxon_at_level(
            result.exp_taxon,
            result.level
        )
        if not ncbi_taxon:
            msg = f"No {result.level} taxon found for {query} in NCBI taxonomy."
            self.logger.error(msg)
            result.error = msg
            return
        marker_enum = Marker(result.ancillary.get('marker_code'))
        translation_table = self.taxonomy_resolver.get_translation_table(marker_enum, ncbi_taxon)
        result.add_ancillary('translation_table',str(translation_table))

class BOLDTaxonomyParser(TaxonomyEnrichmentParser):
    """Parser for BOLD-style FASTA records with process IDs."""
    def __init__(self, config: Config, taxonomy_resolver: TaxonomyResolver):
        super().__init__(config, taxonomy_resolver)
        self.taxonomy = TaxonomicBackbone.BOLD
        self.logger = get_formatted_logger(self.__class__.__name__, config)

    def get_query(self, record: SeqRecord) -> str:
        """
        Returns the process ID from the record, which is the part of the ID before the underscore.
        :param record: A SeqRecord object with process ID in the ID field.
        :return: A process ID string.
        """
        return record.id.split('_')[0]

class NSRTaxonomyParser(TaxonomyEnrichmentParser):
    """Parser for records with taxonomic names that are mapped to NSR."""
    def __init__(self, config: Config, taxonomy_resolver: TaxonomyResolver):
        super().__init__(config, taxonomy_resolver)
        self.taxonomy = TaxonomicBackbone.DWC
        self.logger = get_formatted_logger(self.__class__.__name__, config)

    def get_query(self, record: SeqRecord) -> str:
        """
        Returns the verbatim identification from the record, which is a bcdm_field
        :param record: A SeqRecord object the verbatim identification is stored in.
        :return: A verbatim identification string.
        """
        return record.annotations['bcdm_fields']['identification']