import tempfile
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Optional, Set, List, Tuple
from nbitk.config import Config
from nbitk.Services.Galaxy.BLASTN import BLASTNClient
try:
    # Try importing from module level (older nbitk versions)
    from nbitk.Services.Galaxy.BLASTN import OutputFormat, TaxonomyMethod
except ImportError:
    # Fall back to class attributes (newer nbitk versions)
    OutputFormat = BLASTNClient.OutputFormat
    TaxonomyMethod = BLASTNClient.TaxonomyMethod
from nbitk.Taxon import Taxon
from barcode_validator.resolvers.taxonomy import TaxonResolver
from copy import deepcopy

from .idservice import IDService
from ..dna_analysis_result import DNAAnalysisResult


class GalaxyBLAST(IDService):
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
        self.database = config.get('galaxy_blast.database', 'BOLD species only no duplicates')

    def run_galaxy_blast(self, sequences: List[SeqRecord]) -> Optional[dict]:
        """
        Run local BLAST search constrained by taxonomy and collect results at specified rank.

        :param sequences: List of Bio.SeqRecord objects to search
        :return: Galaxy BLAST search results
        """
        if len(sequences) == 0:
            return None

        # Create temporary files
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as temp_input:

            # Clone records, skip empty ones, and strip dashes from sequences
            cleaned = []
            for sequence in sequences:
                if len(sequence.seq) == 0:
                    continue
                cloned_record = deepcopy(sequence)
                cleaned_sequence = str(cloned_record.seq).replace('-', '').replace('~', '')
                cloned_record.seq = Seq(cleaned_sequence)
                cleaned.append(cloned_record)

            # Prepare input FASTA file
            SeqIO.write(cleaned, temp_input.name, "fasta")
            temp_input_name = temp_input.name

            try:
                # Run BLAST using nbitk wrapper
                result = self.blastn.run_blast(
                    input_file = temp_input_name,
                    databases = [self.database],
                    max_target_seqs = self.max_target_seqs,
                    output_format = OutputFormat.CUSTOM_TAXONOMY,
                    taxonomy_method = TaxonomyMethod.DEFAULT,
                    coverage = 80.0,
                    identity = self.min_identity * 100
                )
                return result

            except Exception as e:
                self.logger.error(f"Error running BLAST: {e}")
                return None

    def parse_blast_result(self, blast_result: str) -> dict[str,Set[Taxon]]:
        """
        Parse BLAST results and collect distinct taxa at specified rank.

        :param blast_result: Path to BLAST result file
        :return: Dictionary of sequence IDs to taxon names
        """

        # Parse BLAST results
        mapping = {}
        with open(blast_result, 'r') as file:
            for line in file:
                columns = line.strip().split('\t')
                if columns:
                    lineage_field = columns[-1]
                    if lineage_field.startswith('#'):

                        # Skip header lines
                        continue

                    # Parse out family and seq ID
                    lineage = lineage_field.split(' / ')
                    family = lineage[-3]
                    seq_id = columns[0].strip()

                    # Instantiate set under focal key
                    if seq_id not in mapping:
                        mapping[seq_id] = set()

                    # Add to set
                    mapping[seq_id].add(Taxon(name=family, taxonomic_rank="family"))

        self.logger.info(f'{len(mapping)} sequences found in BLAST result')
        return mapping


    def identify_batch(self, batch: List[Tuple[DNAAnalysisResult,SeqRecord,Taxon]]) -> None:
        """
        Identify the taxonomic classification of a batch of sequence records. The third item
        in the tuple is a Taxon that limits the query extent. This is only applicable to a local
        BLAST install that is indexed against the NCBI taxonomy. Here it is simply ignored.
        :param batch: List of (DNAAnalysisResult,SeqRecord,Taxon)
        """

        # Create list of sequences to query and run BLAST
        sequences = []
        for item in batch:
            sequences.append(item[1])
        galaxy_output = self.run_galaxy_blast(sequences)

        # Process BLAST result
        mapping = self.parse_blast_result(galaxy_output['blast_output_fasta'])
        for item in batch:
            seq_id = item[1].id
            if seq_id in mapping:
                result = item[0]
                result.obs_taxon = mapping[seq_id]

        # Clean up temporary files
        try:
            os.unlink(galaxy_output['blast_output_fasta'])
            os.unlink(galaxy_output['log_output'])
        except OSError as e:
            self.logger.warning(f"Error cleaning up temporary files: {e}")

    @staticmethod
    def requires_resolver():
        return False

    @staticmethod
    def requires_blastn():
        return False

