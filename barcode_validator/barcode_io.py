from pathlib import Path
from typing import Iterator
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from barcode_validator.dna_analysis_result import DNAAnalysisResult, DNAAnalysisResultSet


class BarcodeIO:
    """
    Handle file I/O operations for barcode validation.
    
    This class centralizes all file operations for the barcode validator,
    including parsing different input formats and writing results.
    
    Examples:
        >>> from nbitk.config import Config
        >>> config = Config()
        >>> io = BarcodeIO(config)
        >>> for record, metadata in io.parse_input('sequences.fasta'):
        ...     print(record.id, metadata)
    """
    
    def __init__(self, config: Config):
        """Initialize with configuration."""
        self.config = config
        self.logger = get_formatted_logger(self.__class__.__name__, config)

    def parse_input(self, file_path: Path) -> Iterator[SeqRecord]:
        """
        Parse input file in either FASTA or tabular format.

        :param file_path: Path to input file
        :return: Iterator of SeqRecord objects
        """
        suffix = Path(file_path).suffix.lower()
        if suffix in ['.fa', '.fasta', '.fna']:
            yield from SeqIO.parse(file_path, 'fasta')
        elif suffix in ['.tsv', '.txt']:
            yield from self._parse_tsv(file_path)
        else:
            raise ValueError(f"Unsupported file format: {suffix}")

    def _parse_tsv(self, file_path: Path) -> Iterator[SeqRecord]:
        """
        Parse CSC-style tabular format with sequence data and metadata.

        Expected columns include:
        - local_id: Primary key for tracking records
        - nuc: The sequence data
        - marker_code: Barcode marker name (e.g. COI-5P)
        - verbatim_identification: Taxonomic identification
        - verbatim_kingdom: Kingdom classification
        - verbatim_rank: Taxonomic rank of identification
        - originial_source: Source of the metadata (i.e. ada or sheets)

        :param file_path: Path to TSV file
        :return: Iterator of SeqRecord objects
        """
        self.logger.info(f"Parsing TSV file: {file_path}")
        with open(file_path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # Extract required fields
                local_id = row.get('local_id')
                sequence = row.get('nuc')

                if not all([local_id, sequence]):
                    self.logger.warning(f"Skipping row with missing required fields")
                    continue

                # Create sequence record
                record = SeqRecord(
                    seq=Seq(sequence),
                    id=local_id,
                    name=local_id,
                    description=row.get('verbatim_identification', '')
                )

                # Add fields to bcdm_fields
                record.annotations['bcdm_fields'] = {
                    'marker_code': row.get('marker_code'),
                    'identification': row.get('verbatim_identification'),
                    'kingdom': row.get('verbatim_kingdom'),
                    'rank': row.get('verbatim_rank'),
                    'source': row.get('originial_source', 'unknown')
                }

                yield record

    def write_results(self, results: DNAAnalysisResultSet, output_path: Path) -> None:
        """
        Write validation results to TSV file.
        
        :param results: Set of validation results
        :param output_path: Path to output file
        """
        self.logger.info(f"Writing results to: {output_path}")
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            
            # Write header
            writer.writerow(DNAAnalysisResult.result_fields())
            
            # Write results
            for result in results.results:
                writer.writerow(result.get_values())

    def write_filtered_fasta(self, results: DNAAnalysisResultSet, 
                           output_path: Path, valid_only: bool = True) -> None:
        """
        Write sequences to FASTA based on validation results.
        
        :param results: Set of validation results
        :param output_path: Path to output FASTA file
        :param valid_only: If True, only write valid sequences
        """
        self.logger.info(f"Writing filtered FASTA to: {output_path}")
        valid_sequences = []
        for result in results.results:
            if not valid_only or result.passes_all_checks():
                if hasattr(result, 'sequence'):  # Need to store sequence in result
                    valid_sequences.append(result.sequence)
                else:
                    self.logger.warning(f"No sequence stored for {result.sequence_id}")
        
        SeqIO.write(valid_sequences, output_path, "fasta")