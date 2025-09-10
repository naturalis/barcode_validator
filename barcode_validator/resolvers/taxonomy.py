from pathlib import Path
from typing import Optional, List, Dict
from Bio.SeqRecord import SeqRecord
from nbitk.Taxon import Taxon
from nbitk.logger import get_formatted_logger
from nbitk.config import Config
from Bio.Phylo import BaseTree
from barcode_validator.constants import TaxonomicRank, TaxonomicBackbone, index_for_rank, rank_by_index

from barcode_validator.dna_analysis_result import DNAAnalysisResult

class TaxonResolver:
    """
    Base class for taxonomic name resolution and reconciliation. This can
    be done in a variety of ways, e.g. by using local taxonomy trees, or
    by web services.
    """

    def __init__(self, config: Config):
        self.logger = get_formatted_logger(self.__class__.__name__, config)
        self.tree: BaseTree = None
        self.parser = None

    def get_type(self) -> TaxonomicBackbone:
        """
        Return the type of taxonomy this resolver is working with.
        :return: TaxonomicBackbone enum value.
        """
        raise NotImplementedError

    def load_tree(self, file: Path) -> None:
        """
        Load a tree from a file.
        :param file: Location of the tree file.
        :return:
        """
        self.parser.file = file.as_posix()
        self.tree = self.parser.parse()

    def check_id(self, id_string: str, node: Taxon) -> bool:
        """
        Checks if the focal node is the one with the given ID.
        :param id_string: ID of the node to check.
        :param node: Focal node, e.g. in a traversal
        :return:
        """
        raise NotImplementedError

    def parse_id(self, record: SeqRecord) -> Optional[str]:
        """
        Parses the taxon ID from a SeqRecord. For example, by stripping
        any suffixes from a BOLD process ID, or by probing the BCDM
        annotations, by using a regular expression on the sequence
        description, and so on.
        :param record: The SeqRecord to parse.
        :return: A string representing the taxon ID.
        """
        raise NotImplementedError

    def enrich_result(self, record: SeqRecord, result: DNAAnalysisResult, rank: TaxonomicRank = TaxonomicRank.FAMILY) -> None:
        """
        Enriches the DNAAnalysisResult object with information from the taxonomy.
        This adds the following properties to the DNAAnalysisResult object:
        - species - Taxon object corresponding to the record's species (or lowest rank)
        - level - Taxonomic rank at which to validate the record
        - exp_taxon Taxon object corresponding to the species' ancestor at the level specified
        :param record: SeqRecord to probe
        :param result: DNAAnalysisResult to enrich
        :param rank: Taxonomic rank at which to validate the record
        """
        result.level = rank.value

        # Attempt to parse out the ID and look up the taxon
        tid = self.parse_id(record)
        if not tid:
            result.error = f"Could not parse ID from record {record.id}"
            self.logger.warning(result.error)
            return

        # Attempt to find the corresponding taxon for the ID
        nodes = self.find_nodes(tid)
        if len(nodes) == 0:
            result.error = f"Could not find nodes for ID {tid}"
            self.logger.warning(result.error)
            return
        elif len(nodes) > 1:

            # This is not ideal, e.g. if there are homonyms, but we can try to proceed
            self.logger.warning(f"Found multiple nodes for ID {tid}: {nodes}")
        result.species = nodes

        # Attempt to find the validation taxon
        exp_taxa = []
        for species in result.species:

            # There may be cases where the provided identification is at a higher rank than the validation rank.
            # For example, some specimens are identified as 'Lepidoptera (order)' only, but we may want to validate
            # at family rank. In such cases, we validate at the provided rank instead and issue a warning.
            validation_rank = rank
            provided_rank = TaxonomicRank(species.taxonomic_rank.lower())
            if index_for_rank(provided_rank) < index_for_rank(validation_rank):
                self.logger.warning(f"{species.name} is rank '{provided_rank.value}', i.e. higher than '{rank.value}'")
                validation_rank = provided_rank
                result.level = validation_rank

            # Now find the ancestor at the desired rank
            exp_taxon = self.find_ancestor_at_rank(species, validation_rank)
            if exp_taxon is None:
                self.logger.warning(f"Could not find ancestor of species {result.species} at rank {validation_rank}")
                continue
            exp_taxa.append(exp_taxon)

        result.exp_taxon = exp_taxa

    def find_nodes(self, id_string: str) -> List[Taxon]:
        """
        Find nodes in the tree by an ID. What constitutes an ID depends
        on the type of taxonomy, so matching this is delegated to the subclass.
        :param id_string: ID of the node to find.
        :return: List of Taxon objects.
        """
        nodes = []
        for node in self.tree.find_clades():
            if self.check_id(id_string, node):
                nodes.append(node)
        return nodes

    def match_nodes(self, node_to_match: Taxon, properties: set[str] = None) -> List[Taxon]:
        """
        Match the focal node, which comes from a different taxonomy, to nodes in the tree.
        In other words:  a type of reconciliation between two taxonomies.
        The matching is based on literal matches of the properties
        name and taxonomic_rank by default, but other criteria can be specified.
        These are the names of object properties of the Taxon class and any of
        its superclasses.
        :param node_to_match: Node to match. This node is from a different tree.
        :param properties: Properties to match.
        :return: Optional Taxon object.
        """
        all_properties = {"name", "taxonomic_rank"}
        if properties:
            all_properties.update(properties)

        nodes = []
        for node in self.tree.find_clades():
            matches = 0
            for prop in all_properties:

                # Either node doesn't have the property
                if not hasattr(node, prop) or not hasattr(node_to_match, prop):
                    break

                # The property value matches, otherwise break
                if getattr(node, prop) == getattr(node_to_match, prop):
                    matches += 1
                else:
                    break

            # All provided properties must match
            if matches == len(all_properties):
                nodes.append(node)
        return nodes

    def find_ancestor_at_rank(self, node: Taxon, rank: TaxonomicRank) -> Optional[Taxon]:
        """
        Find the ancestor of the focal node at the given rank.
        :param node: Node whose ancestor we want to find.
        :param rank: The rank of the ancestor we want to find.
        :return: Optional Taxon object.
        """
        for node in self.tree.root.get_path(node):
            if str(node.taxonomic_rank).lower() == rank.value.lower():
                return node
        return None

    def get_lineage_dict(self, nodes: List[Taxon], ranks: Optional[List[TaxonomicRank]] = None) -> Dict[str, str]:
        """
        Get a dictionary of rank:name pairs along the path to root.
        :param nodes: The nodes whose lineage to retrieve
        :param ranks: Optional list of ranks to include. If None, include all.
        :return: Dictionary mapping ranks to names along the path to root
        """
        lineage = {}

        for node in nodes:
            for ancestor in self.tree.root.get_path(node):
                rank = str(ancestor.taxonomic_rank).lower()
                if not ranks or rank in ranks:
                    if rank not in lineage:
                        lineage[rank] = []
                    lineage[rank].append(ancestor.name)
        return lineage