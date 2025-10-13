import tempfile
import os
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Optional, Set, List, Tuple
from nbitk.config import Config
from nbitk.Tools import Blastn
from nbitk.Taxon import Taxon
from barcode_validator.resolvers.taxonomy import TaxonResolver
from barcode_validator.constants import TaxonomicRank
from .idservice import IDService
from ..dna_analysis_result import DNAAnalysisResult
from copy import deepcopy


class BLAST(IDService):
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

    def run_localblast(self, sequences: List[SeqRecord], constraint: int) -> Optional[str]:
        """
        Run local BLAST search constrained by taxonomy and collect results at specified rank.

        :param sequences: List of Bio.SeqRecord object to search
        :param constraint: NCBI taxon ID to constrain search
        :return: List of distinct taxa at specified rank, or None if search fails
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
            self.blastn.set_query(temp_input_name)
            if constraint != 0:
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

    def parse_blast_result(self, blast_result: str) -> dict[str,Set[Taxon]]:
        """
        Parse BLAST results and collect distinct taxa at specified rank.

        :param blast_result: Path to BLAST result file
        :return: List of distinct taxa at specified rank
        """

        # Parse BLAST results
        matches = {}
        seen = {}
        with open(blast_result, 'r') as file:
            for line in file:
                columns = line.strip().split('\t')
                if columns:
                    seq_id = columns[0]

                    # process and/or cache match
                    match_id = columns[1]
                    bin_uri = match_id.split('|')[1]
                    if bin_uri not in seen:
                        seen[bin_uri] = self.taxonomy_resolver.find_nodes(bin_uri)

                    matches[seq_id] = seen[bin_uri]


        self.logger.info(f'{len(matches)} distinct matches found in BLAST result')

        return matches

    def collect_higher_taxa(self, matches: Set[Taxon], level: TaxonomicRank) -> Set[Taxon]:
        """
        Collect distinct higher taxa at specified rank from set of taxids.

        :param matches: Set of Taxon objects
        :param level: Taxonomic level to collect
        :return: List of distinct taxa at specified rank
        """

        # Collect distinct taxa at specified rank
        taxa = set()
        for tip in matches:
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
        blast_report = self.run_localblast([record], constraint)
        distinct_taxids = self.parse_blast_result(f"{blast_report}.tsv")
        higher_taxa = self.collect_higher_taxa(distinct_taxids, level)

        # Clean up temporary files
        try:
            os.unlink(blast_report)
            os.unlink(f"{blast_report}.tsv")
        except OSError as e:
            self.logger.warning(f"Error cleaning up temporary files: {e}")

        return higher_taxa

    def identify_batch(self, batch: List[Tuple[DNAAnalysisResult, SeqRecord, Taxon]]) -> None:
        """
        Identify the taxonomic classification of a batch of sequence records. The third item
        in the tuple is a Taxon that limits the query extent. This is only applicable to a local
        BLAST install that is indexed against the NCBI taxonomy. Here it is simply ignored.
        """

        # Create list of sequences to query and run BLAST
        sequences = []
        for item in batch:
            sequences.append(item[1])

        # Run local BLAST
        blast_output = self.run_localblast(sequences, 0)
        matches = self.parse_blast_result(f"{blast_output}.tsv")

        # Process BLAST result
        for item in batch:
            result = item[0]
            seq_id = item[1].id
            if seq_id in matches:
                level = result.level if result.level else TaxonomicRank.FAMILY
                result.obs_taxon = self.collect_higher_taxa(matches[seq_id], level)

        # Clean up temporary files
        try:
            os.unlink(blast_output)
            os.unlink(f"{blast_output}.tsv")
        except OSError as e:
            self.logger.warning(f"Error cleaning up temporary files: {e}")

    @staticmethod
    def requires_resolver():
        return True

    @staticmethod
    def requires_blastn():
        return True

    def set_blastn(self, blastn) -> None:
        """
        Assign BLASTN object to the IDService. This must have been configured by the orchestrator.
        :param blastn: A fully configured instance of Blastn class.
        """
        self.blastn = blastn
