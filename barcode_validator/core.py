import logging
from typing import List, Optional
from barcode_validator.config import Config
from barcode_validator.taxonomy import read_ncbi_taxonomy, read_bold_taxonomy, get_tip_by_processid, run_localblast
from barcode_validator.alignment import parse_fasta, align_to_hmm, translate_sequence, get_stop_codons, \
    marker_seqlength, num_ambiguous
from barcode_validator.result import DNAAnalysisResult
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.BaseTree import Tree


class BarcodeValidator:
    def __init__(self, config: Config):
        """
        Initialize the BarcodeValidator object.
        :param config: An instance of the Config class
        """
        self.config = config
        self.ncbi_tree: Optional[Tree] = None
        self.bold_tree: Optional[Tree] = None

    def initialize(self) -> None:
        """
        Initialize the taxonomy trees. This is separate from __init__ to allow for lazy loading.
        :return: None
        """
        logging.info("Initializing taxonomy trees...")
        self.ncbi_tree = read_ncbi_taxonomy(self.config.get('ncbi_taxonomy'))
        self.bold_tree = read_bold_taxonomy(self.config.get('bold_sheet_file'))
        logging.info("Initialization complete.")

    def validate_fasta(self, fasta_file_path: str) -> List[DNAAnalysisResult]:
        """
        Validate a FASTA file of DNA sequences.
        :param fasta_file_path: A path to a FASTA file
        :return: A list of DNAAnalysisResult objects
        """
        results = []
        for process_id, record in parse_fasta(fasta_file_path):
            result = self.validate_record(process_id, record)
            results.append(result)
        return results

    def validate_record(self, process_id: str, record: SeqRecord) -> DNAAnalysisResult:
        """
        Validate a single DNA sequence record.
        :param process_id: A process ID
        :param record: A Bio.SeqRecord object
        :return: A DNAAnalysisResult object
        """

        # Instantiate result object with process ID and calculate full sequence stats
        result = DNAAnalysisResult(process_id)
        result.full_length = len(record.seq)
        result.full_ambiguities = num_ambiguous(record)

        # Lookup expected species and higher taxon, infer observed higher taxa from blast
        result.species = get_tip_by_processid(process_id, self.bold_tree)
        logging.info(f"Species: {result.species}")
        for node in self.bold_tree.root.get_path(result.species):
            if node.taxonomic_rank == self.config.get('level'):
                result.exp_taxon = node
                break
        result.obs_taxon = run_localblast(record, self.ncbi_tree, self.bold_tree)

        # Compute marker quality metrics
        aligned_sequence = align_to_hmm(record)
        amino_acid_sequence = translate_sequence(aligned_sequence)
        result.stop_codons = get_stop_codons(amino_acid_sequence)
        result.seq_length = marker_seqlength(aligned_sequence)
        result.ambiguities = num_ambiguous(aligned_sequence)

        # Return the result object
        return result
