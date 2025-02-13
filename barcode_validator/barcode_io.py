from pathlib import Path
import json
from typing import Iterator, Tuple, Optional, Dict
from Bio import SeqIO
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

    def parse_input(self, file_path: Path) -> Iterator[Tuple[SeqRecord, Optional[Dict]]]:
        """
        Parse input file in either FASTA or tabular format.
        
        :param file_path: Path to input file
        :return: Iterator of (sequence_record, metadata) tuples
        """
        suffix = Path(file_path).suffix.lower()
        if suffix in ['.fa', '.fasta', '.fna']:
            yield from self._parse_fasta(file_path)
        elif suffix in ['.tsv', '.txt']:
            yield from self._parse_tsv(file_path)
        else:
            raise ValueError(f"Unsupported file format: {suffix}")

    def _parse_fasta(self, file_path: Path) -> Iterator[Tuple[SeqRecord, Optional[Dict]]]:
        """
        Parse FASTA file with optional JSON configuration in descriptions.
        
        :param file_path: Path to FASTA file
        :return: Iterator of (sequence_record, config) tuples
        """
        self.logger.info(f"Parsing FASTA file: {file_path}")
        for record in SeqIO.parse(file_path, 'fasta'):
            config = None
            json_start = record.description.find('{')
            if json_start != -1:
                try:
                    json_str = record.description[json_start:]
                    config = json.loads(json_str)
                    # Remove JSON from description
                    record.description = record.description[:json_start].strip()
                except json.JSONDecodeError as e:
                    self.logger.warning(f"Failed to parse JSON for {record.id}: {e}")
            
            # Add process ID to record annotations
            if 'bcdm_fields' not in record.annotations:
                record.annotations['bcdm_fields'] = {}
            record.annotations['bcdm_fields']['processid'] = record.id.split('_')[0]
            
            yield record, config

    def _parse_tsv(self, file_path: Path) -> Iterator[Tuple[SeqRecord, Dict]]:
        """
        Parse tabular format with sequence data and metadata.
        
        :param file_path: Path to TSV file
        :return: Iterator of (sequence_record, metadata) tuples
        """
        self.logger.info(f"Parsing TSV file: {file_path}")
        with open(file_path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if 'sequence_id' not in row or 'nuc' not in row:
                    raise ValueError("TSV must contain 'sequence_id' and 'nuc' columns")
                
                # Create sequence record
                record = SeqRecord(
                    seq=row['nuc'],
                    id=row['sequence_id'],
                    name=row['sequence_id'],
                    description=row.get('identification', '')
                )
                
                # Create metadata dict from remaining columns
                metadata = {k: v for k, v in row.items() 
                          if k not in ['sequence_id', 'nuc']}
                
                yield record, metadata

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