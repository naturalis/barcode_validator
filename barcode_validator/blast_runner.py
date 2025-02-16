import tempfile
import os
from Bio import SeqIO
from Bio.Phylo.BaseTree import Tree
from typing import Optional, List
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from nbitk.Tools import Blastn


class BlastRunner:
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

    def __init__(self, config: Config):
        self.logger = get_formatted_logger(self.__class__.__name__, config)

        # TODO: it is probably better if this is a reference to a TaxonomyResolver instance
        self.ncbi_tree: Optional[Tree] = None

        # Initialize BLASTN using nbitk wrapper
        self.blastn = Blastn(config)
        self.blastn.set_db(config.get('blast_db'))
        self.blastn.set_num_threads(config.get('num_threads'))
        self.blastn.set_evalue(config.get('evalue'))
        self.blastn.set_max_target_seqs(config.get('max_target_seqs'))
        self.blastn.set_word_size(config.get('word_size'))
        self.blastn.set_task('megablast')
        self.blastn.set_outfmt("6 qseqid sseqid pident length qstart qend sstart send evalue bitscore staxids")

        # Immediately check for the environment variables
        blastdb = os.environ.get("BLASTDB")
        lmdb_map_size = os.environ.get("BLASTDB_LMDB_MAP_SIZE")

        if not blastdb:
            self.logger.warning("Environment variable 'BLASTDB' is not set.")
        if not lmdb_map_size:
            self.logger.warning("Environment variable 'BLASTDB_LMDB_MAP_SIZE' is not set.")


    def run_localblast(self, sequence, constraint: int, level='family') -> Optional[List]:
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

            return self.parse_blast_result(f"{temp_input_name}.tsv", level)

        except Exception as e:
            self.logger.error(f"Error running BLAST: {e}")
            return None

        finally:
            # Clean up temporary files
            try:
                os.unlink(temp_input_name)
                os.unlink(f"{temp_input_name}.tsv")
            except OSError as e:
                self.logger.warning(f"Error cleaning up temporary files: {e}")

    def parse_blast_result(self, blast_result: str, level: str) -> Optional[List]:
        """
        Parse BLAST results and collect distinct taxa at specified rank.

        :param blast_result: Path to BLAST result file
        :param level: Taxonomic level to collect
        :return: List of distinct taxa at specified rank
        """
        if not self.ncbi_tree:
            self.logger.error("NCBI taxonomy tree not initialized")
            return None

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

        return self.collect_higher_taxa(distinct_taxids, level)

    def collect_higher_taxa(self, taxids: set, level: str) -> List:
        """
        Collect distinct higher taxa at specified rank from set of taxids.

        :param taxids: Set of NCBI taxonomy IDs
        :param level: Taxonomic level to collect
        :return: List of distinct taxa at specified rank
        """
        # Collect tips with matching taxids
        tips = [tip for tip in self.ncbi_tree.get_terminals()
                if tip.guids.get('taxon') in taxids]

        self.logger.info(f'Found {len(tips)} tips for {len(taxids)} taxids in tree')

        # Collect distinct taxa at specified rank
        taxa = []
        for tip in tips:
            for node in self.ncbi_tree.root.get_path(tip):
                if node.taxonomic_rank == level and node not in taxa:
                    taxa.append(node)
                    self.logger.info(f"Found ancestor '{node}' for '{tip}'")

        self.logger.info(f'Collected {len(taxa)} higher taxa')
        return taxa