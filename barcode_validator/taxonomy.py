import tempfile
import subprocess
import os
from Bio import SeqIO
from Bio.Phylo.BaseTree import Tree
from typing import Optional
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from nbitk.Tools import Blastn


class BlastRunner:

    def __init__(self, config: Config):
        """
        Initialize the BarcodeValidator object.
        """
        class_name = self.__class__.__name__
        self.logger = get_formatted_logger(class_name, config)
        self.ncbi_tree: Optional[Tree] = Optional[Tree]

        # Initialize the BLASTN tool
        self.blastn = Blastn(config)
        self.blastn.set_db(config.get('blast_db'))
        self.blastn.set_num_threads(config.get('num_threads'))
        self.blastn.set_evalue(config.get('evalue'))
        self.blastn.set_max_target_seqs(config.get('max_target_seqs'))
        self.blastn.set_word_size(config.get('word_size'))
        self.blastn.set_task('megablast')
        self.blastn.set_outfmt("6 qseqid sseqid pident length qstart qend sstart send evalue bitscore staxids")

        # Initialize environment variables
        #self.BLASTDB_LMDB_MAP_SIZE: Optional[int] = config.get('BLASTDB_LMDB_MAP_SIZE')
        self.BLASTDB: Optional[str] = config.get('BLASTDB')

    def run_localblast(self, sequence, constraint, level='family'):
        """
        Run a local BLASTN search against a local database and return the taxonomic lineages of the hits.
        :param sequence: A Bio.SeqRecord object
        :param constraint: An NCBI taxon ID to constrain the search to, e.g. only search within class Insecta's taxon ID
        :param level: The taxonomic level at which to collect higher taxa, e.g. only return distinct families
        :return: A list of distinct higher taxa at the specified rank
        """
        self.logger.info("Running local BLASTN...")

        # Check if there is a sequence
        if len(sequence.seq) == 0:
            return None

        # Create a temporary file for the input sequence
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as temp_input:
            SeqIO.write(sequence, temp_input, "fasta")
            temp_input_name = temp_input.name

        # Run local BLASTN
        #os.environ['BLASTDB_LMDB_MAP_SIZE'] = str(self.BLASTDB_LMDB_MAP_SIZE)
        os.environ['BLASTDB'] = str(self.BLASTDB)
        try:
            self.blastn.set_query(temp_input_name)
            self.blastn.set_taxids(constraint)
            self.blastn.set_out(f"{temp_input_name}.tsv")
            return_code = self.blastn.run()
            if return_code != 0:
                raise subprocess.CalledProcessError(return_code, 'blastn')
            return self.parse_blast_result(f"{temp_input_name}.tsv", level)

        # Handle exception
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Error running local BLASTN: {e}")
            raise

    def parse_blast_result(self, blast_result, level):
        """
        Parse the BLAST result file and return the distinct higher taxa that sit at the configured taxonomic level.
        :param blast_result: A path to the BLAST result file
        :param level: The taxonomic level at which to collect higher taxa
        :return: A list of distinct higher taxa at the specified rank
        """
        # Parse BLAST result
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

    def collect_higher_taxa(self, taxids, level):
        """
        Collect the distinct higher taxa that sit at the configured taxonomic level for an input set of NCBI taxids. For
        example, if the taxonomic level is 'family', the function will return a list of distinct families that are
        encountered in the traversal of the tree from the tips with the specified taxids to the root.
        :param taxids: A set of NCBI taxids
        :param level: The taxonomic level at which to collect higher taxa
        :return: A list of distinct higher taxa at the specified rank
        """
        tips = []

        # Iterate over all tips in the NCBI taxonomy tree
        # to collect all tips with the specified taxids
        for tip in self.ncbi_tree.get_terminals():
            taxid = tip.guids['taxon']

            # Focal tip is annotated with an NCBI taxon ID in the provided lists
            if taxid in taxids:
                tips.append(tip)
                self.logger.debug(f'Found tip {tip.name} with taxid {taxid}')
        self.logger.info(f'Found {len(tips)} tips for {len(taxids)} in the tree')

        # Iterate over the collected tips to find their lineages to build a set
        # of distinct higher taxa with the specified rank
        taxa = []
        for tip in tips:
            for node in self.ncbi_tree.root.get_path(tip):
                self.logger.debug(f'Traversing {node} from lineage {tip}')
                if node.taxonomic_rank == level:
                    if node not in taxa:  # Check for uniqueness
                        taxa.append(node)
                        self.logger.info(f"Found ancestor '{node}' for '{tip}'")
        self.logger.info(f'Collected {len(taxa)} higher taxa')
        return list(taxa)
