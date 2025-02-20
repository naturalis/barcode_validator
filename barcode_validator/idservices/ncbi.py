import tempfile
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.BaseTree import Tree
from typing import Optional, List
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from nbitk.Tools import Blastn
from nbitk.Taxon import Taxon
from barcode_validator.taxonomy_resolver import TaxonomicRank, TaxonomyResolver
from .idservice import IDService


class NCBI(IDService):
    """
    Runs BLAST searches and parses results for taxonomic validation.

    This class wraps around the NCBI BLAST+ wrapper `nbitk.Tools.Blastn`, enhancing
    its functionality by some additional parsing and aggregation logic in dealing
    with tabular BLAST results. Specifically, it aggregates the taxon results returned
    across hits at a higher taxonomic level (usually family) to obtain a set of distinct
    Taxon objects. In the overall logic of barcode validation via reverse taxonomy,
    this set is checked to see if it contains the expected taxon (i.e. that which was
    provided by the collector).
    """

    def __init__(self, config: Config, blastn: Blastn, taxonomy_resolver: TaxonomyResolver):
        self.logger = get_formatted_logger(self.__class__.__name__, config)

        # Lookup and traverse the NCBI taxonomy tree
        self.taxonomy_resolver = taxonomy_resolver

        # Configure BLASTN beyond what the orchestrator has done
        self.blastn = blastn
        self.blastn.set_task('megablast')
        self.blastn.set_outfmt("6 qseqid sseqid pident length qstart qend sstart send evalue bitscore staxids")

        # Immediately check for the environment variables
        blastdb = os.environ.get("BLASTDB")
        lmdb_map_size = os.environ.get("BLASTDB_LMDB_MAP_SIZE")

        if not blastdb:
            self.logger.warning("Environment variable 'BLASTDB' is not set.")
        if not lmdb_map_size:
            self.logger.warning("Environment variable 'BLASTDB_LMDB_MAP_SIZE' is not set.")


    def run_localblast(self, sequence: SeqRecord, constraint: int, level: TaxonomicRank = TaxonomicRank.FAMILY) -> str:
        """
        Run local BLAST search constrained by taxonomy and collect results at specified rank.

        :param sequence: Bio.SeqRecord object to search
        :param constraint: NCBI taxon ID to constrain search
        :param level: Taxonomic level for collecting results
        :return: List of distinct taxa at specified rank, or None if search fails
        """
        if len(sequence.seq) == 0:
            return None

        # Create temporary files
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as temp_input:
            SeqIO.write(sequence, temp_input, "fasta")
            temp_input_name = temp_input.name

        try:
            # Run BLAST using nbitk wrapper
            self.blastn.set_query(temp_input_name)
            self.blastn.set_taxids([str(constraint)])
            self.blastn.set_out(f"{temp_input_name}.tsv")

            return_code = self.blastn.run()
            if return_code != 0:
                self.logger.error(f"BLAST search failed with return code {return_code}")
                return None
            return temp_input_name

        except Exception as e:
            self.logger.error(f"Error running BLAST: {e}")
            return None

    def parse_blast_result(self, blast_result: str, level: TaxonomicRank) -> set:
        """
        Parse BLAST results and collect distinct taxa at specified rank.

        :param blast_result: Path to BLAST result file
        :param level: Taxonomic level to collect
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
                    distinct_taxids.update(taxid.strip() for taxid in taxids if taxid.strip())

        self.logger.info(f'{len(distinct_taxids)} distinct taxids found in BLAST result')
        self.logger.debug(distinct_taxids)

        return distinct_taxids

    def collect_higher_taxa(self, taxids: set, level: TaxonomicRank) -> List:
        """
        Collect distinct higher taxa at specified rank from set of taxids.

        :param taxids: Set of NCBI taxonomy IDs
        :param level: Taxonomic level to collect
        :return: List of distinct taxa at specified rank
        """
        # Collect tips with matching taxids
        tips = []
        for taxon_id in taxids:
            taxon = self.taxonomy_resolver.get_taxon_by_id(taxon_id)
            if taxon:
                tips.append(taxon)
        self.logger.info(f'Found {len(tips)} tips for {len(taxids)} taxids in tree')

        # Collect distinct taxa at specified rank
        taxa = set()
        for tip in tips:
            taxon = self.taxonomy_resolver.get_constraint_taxon(tip, level)
            if taxon:
                taxa.add(taxon)
                self.logger.debug(f"Found ancestor '{taxon}' for '{tip}'")

        self.logger.info(f'Collected {len(taxa)} higher taxa')
        return taxa

    def identify_record(self, record: SeqRecord, level: TaxonomicRank = TaxonomicRank.FAMILY, extent: Taxon = None) -> List[Taxon]:
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
        blast_report = self.run_localblast(record, extent.guids.get('taxon', 33208), level)
        distinct_taxids = self.parse_blast_result(f"{blast_report}.tsv", level)
        higher_taxa = self.collect_higher_taxa(distinct_taxids, level)

        # Clean up temporary files
        try:
            os.unlink(blast_report)
            os.unlink(f"{blast_report}.tsv")
        except OSError as e:
            self.logger.warning(f"Error cleaning up temporary files: {e}")

        return higher_taxa
