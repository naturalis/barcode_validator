import tempfile
import os
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Optional, Set
from nbitk.config import Config
from nbitk.Tools import Blastn
from nbitk.Taxon import Taxon
from barcode_validator.resolvers.taxonomy import TaxonResolver
from barcode_validator.constants import TaxonomicRank
from .idservice import IDService
from copy import deepcopy


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

    def __init__(self, config: Config):
        super().__init__(config)
        self.taxonomy_resolver: TaxonResolver = Optional[TaxonResolver]
        self.blastn: Blastn = Optional[Blastn]

        # Set the BLASTDB environment variable if not already set and if blast_db is in config
        blast_db = config.get('reference_library.database_path', None)

        if blast_db is not None:
            # Get the directory containing the blast database files
            blast_db_dir = str(Path(blast_db).parent)

            # Check if BLASTDB environment variable is set
            if 'BLASTDB' not in os.environ:
                # Set BLASTDB to the database directory
                os.environ['BLASTDB'] = blast_db_dir

        if blast_db is None:
            self.logger.warning("Config variable `blast_db` is not set. BLAST searches may fail.")

    def run_localblast(self, sequence: SeqRecord, constraint: int) -> Optional[str]:
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

    def identify_record(self, record: SeqRecord, level: TaxonomicRank = TaxonomicRank.FAMILY, extent: Taxon = None) -> Set[Taxon]:
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
        blast_report = self.run_localblast(record, constraint)
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
        return True

    def set_blastn(self, blastn) -> None:
        """
        Assign BLASTN object to the IDService. Configure it further. Test environment variables.
        :param blastn: an instance of Blastn class.
        """
        self.blastn = blastn

        # Configure BLASTN beyond what the orchestrator has done
        self.blastn.set_task('megablast')
        self.blastn.set_outfmt("6 qseqid sseqid pident length qstart qend sstart send evalue bitscore staxids")

        # Immediately check for the environment variables
        blastdb = os.environ.get("BLASTDB")
        lmdb_map_size = os.environ.get("BLASTDB_LMDB_MAP_SIZE")
        if not blastdb:
            self.logger.warning("Environment variable 'BLASTDB' is not set.")
        if not lmdb_map_size:
            self.logger.warning("Environment variable 'BLASTDB_LMDB_MAP_SIZE' is not set.")
