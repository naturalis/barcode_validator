import tempfile
import os
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Optional, Set
from nbitk.config import Config
from nbitk.Services.Galaxy.BLASTN import BLASTNClient, OutputFormat, TaxonomyMethod
from nbitk.Taxon import Taxon
from barcode_validator.resolvers.taxonomy import TaxonResolver
from barcode_validator.constants import TaxonomicRank
from .blast import BLAST
from copy import deepcopy


class GalaxyBLAST(BLAST):
    """
    Runs BLAST searches and parses results for taxonomic validation.

    This class wraps around the Galaxy BLAST service `nbitk.Services.Galaxy.BLASTN.BLASTNClient`, enhancing
    its functionality by some additional parsing and aggregation logic in dealing with tabular BLAST results.
    Specifically, it aggregates the taxon results returned across hits at a higher taxonomic level (usually family)
    to obtain a set of distinct Taxon objects. In the overall logic of barcode validation via reverse taxonomy,
    this set is checked to see if it contains the expected taxon (i.e. that which was provided by the collector).

    Requires the following environment variables to be set:

    GALAXY_DOMAIN: The domain of the Galaxy instance to use for BLAST searches, e.g. galaxy.naturalis.nl
    GALAXY_API_KEY: The API key for the Galaxy instance, used for authentication
    """

    def __init__(self, config: Config):
        super().__init__(config)
        self.taxonomy_resolver: TaxonResolver = Optional[TaxonResolver]
        self.blastn: BLASTNClient = BLASTNClient(config, self.logger)

    def run_localblast(self, sequence: SeqRecord, constraint: int) -> Optional[dict]:
        """
        Run local BLAST search constrained by taxonomy and collect results at specified rank.

        :param sequence: Bio.SeqRecord object to search
        :param constraint: NCBI taxon ID to constrain search
        :return: List of distinct taxa at specified rank, or None if search fails
        """
        if len(sequence.seq) == 0:
            return None

        # Create temporary files
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as temp_input:

            # Clone record and strip dashes from sequence
            cloned_record = deepcopy(sequence)
            cleaned_sequence = str(cloned_record.seq).replace('-', '').replace('~', '')
            cloned_record.seq = Seq(cleaned_sequence)

            SeqIO.write(cloned_record, temp_input, "fasta")
            temp_input_name = temp_input.name

        try:
            # Run BLAST using nbitk wrapper
            result = self.blastn.run_blast(
                input_file = temp_input_name,
                databases = ["Genbank CO1 (2023-11-15)"],
                max_target_seqs = self.max_target_seqs,
                output_format = OutputFormat.TABULAR,
                taxonomy_method = TaxonomyMethod.DEFAULT,
                coverage = 80.0,
                identity = self.min_identity * 100
            )
            return result

        except Exception as e:
            self.logger.error(f"Error running BLAST: {e}")
            return None

    def parse_blast_result(self, blast_result: str) -> Set[int]:
        """
        Parse BLAST results and collect distinct taxa at specified rank.

        :param blast_result: Path to BLAST result file
        :return: List of distinct taxa at specified rank
        """

        # Parse BLAST results
        distinct_taxids = set()
        with open(blast_result, 'r') as file:
            for line in file:
                columns = line.strip().split('\t')
                if columns:
                    taxid_field = columns[-1]
                    taxids = taxid_field.split(';')
                    distinct_taxids.update(int(taxid.strip()) for taxid in taxids if taxid.strip())

        self.logger.info(f'{len(distinct_taxids)} distinct taxids found in BLAST result')
        self.logger.debug(distinct_taxids)

        return distinct_taxids

    def collect_higher_taxa(self, taxids: Set[int], level: TaxonomicRank) -> Set[Taxon]:
        """
        Collect distinct higher taxa at specified rank from set of taxids.

        :param taxids: Set of NCBI taxonomy IDs
        :param level: Taxonomic level to collect
        :return: List of distinct taxa at specified rank
        """
        # Collect tips with matching taxids
        tips = []
        for taxon_id in taxids:
            taxon = self.taxonomy_resolver.find_nodes(str(taxon_id))
            if len(taxon) == 0:
                self.logger.warning(f"No taxon found for taxon ID '{taxon_id}'")
                continue
            elif len(taxon) > 1:
                self.logger.warning(f"Multiple taxons found for taxon ID '{taxon_id}': {taxon}")
                continue
            else:
                taxon = taxon[0]
                self.logger.debug(f"Found taxon '{taxon}' for taxon ID '{taxon_id}'")
                tips.append(taxon)
        self.logger.info(f'Found {len(tips)} tips for {len(taxids)} taxids in tree')

        # Collect distinct taxa at specified rank
        taxa = set()
        for tip in tips:
            taxon = self.taxonomy_resolver.find_ancestor_at_rank(tip, level)
            if taxon:
                taxa.add(taxon)
                self.logger.debug(f"Found ancestor '{taxon}' for '{tip}'")
            else:
                self.logger.warning(f"No {level} ancestor found for '{tip}'")
        self.logger.info(f'Collected {len(taxa)} higher taxa')
        return taxa

    def identify_record(self, record: SeqRecord, level: TaxonomicRank = TaxonomicRank.FAMILY, extent: Taxon = None) -> \
    Set[Taxon]:
        """
        Identify the taxonomic classification of a sequence record using BLAST.

        Implementation of the IDService interface that uses BLAST to identify sequences
        against the NCBI nucleotide database.

        :param record: A Bio.SeqRecord object containing the sequence to identify
        :param level: The taxonomic rank at which to return results (default: 'family')
        :param extent: A Taxon object representing the extent of the search (default: None)
        :return: A list of Taxon objects representing the observed higher taxa at the specified level
        """
        if extent:
            constraint = extent.guids.get('taxon')
        else:
            constraint = 33208  # Metazoa

        # Run BLAST
        results = self.run_localblast(record, constraint)
        blast_report = results['blast_output_fasta']
        distinct_taxids = self.parse_blast_result(f"{blast_report}.tsv")
        higher_taxa = self.collect_higher_taxa(distinct_taxids, level)

        # Clean up temporary files
        try:
            os.unlink(blast_report)
            os.unlink(f"{blast_report}.tsv")
        except OSError as e:
            self.logger.warning(f"Error cleaning up temporary files: {e}")

        return higher_taxa

    @staticmethod
    def requires_resolver():
        return True

    @staticmethod
    def requires_blastn():
        return False

