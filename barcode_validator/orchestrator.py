import csv
from pathlib import Path
from typing import Optional, Iterator

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from nbitk.SeqIO.BCDM import BCDMIterator
from .dna_analysis_result import DNAAnalysisResult, DNAAnalysisResultSet
from .taxonomic_enrichment import NSRTaxonomyParser, BOLDTaxonomyParser
from .taxonomy_resolver import Marker, TaxonomyResolver, TaxonomicBackbone
from .validators.non_coding import NonCodingValidator
from .validators.taxonomic import TaxonomicValidator
from .validators.protein_coding import ProteinCodingValidator

class ValidationOrchestrator:
    """
    Main orchestrator for DNA barcode validation.

    This class coordinates the complete validation process including:
    - Input parsing (FASTA/TSV)
    - Validator initialization and coordination
    - Taxonomic resolution
    - Result set management
    - Output generation

    The orchestrator maintains responsibility for the entire validation
    lifecycle to ensure consistent data flow and error handling.

    Examples:
        >>> from nbitk.config import Config
        >>> config = Config()
        >>> config.load_config('/path/to/config.yaml')
        >>> orchestrator = ValidationOrchestrator(config)
        >>> results = orchestrator.validate_file('sequences.fasta')
        >>> orchestrator.write_results(results, 'output.tsv')
    """

    def __init__(self, config: Config):
        """
        Initialize the orchestrator.

        :param config: Configuration object containing validation parameters
        """
        self.config = config
        self.logger = get_formatted_logger(self.__class__.__name__, config)
        self.taxonomy_resolver = TaxonomyResolver(config)
        self.taxonomy_parser = None
        self.structural_validator = None
        self.taxonomic_validator = None

    def validate_file(self, input_path: Path,
                     csv_path: Optional[Path] = None,
                     yaml_path: Optional[Path] = None) -> DNAAnalysisResultSet:
        """
        Validate sequences from an input file.

        :param input_path: Path to input file (FASTA/TSV)
        :param csv_path: Optional path to CSV with record-level analytics
        :param yaml_path: Optional path to YAML with analysis-level config
        :return: Set of validation results
        :raises RuntimeError: If orchestrator not initialized
        :raises ValueError: If input file format not supported
        """
        # We defer initialization until we definitely need it, because some of the validators
        # need assets and setup time that are quite costly.
        backbone_type = TaxonomicBackbone(self.config.get('taxonomic_backbone'))
        marker_type = Marker(self.config.get('marker'))
        self._initialize(backbone_type, marker_type)

        # Validate records and create result set
        results = []
        for record in self._parse_input(input_path):
            result = self._validate_record(record, str(input_path), backbone_type, marker_type)
            results.append(result)
        result_set = DNAAnalysisResultSet(results)

        # Add additional data if provided
        if csv_path:
            result_set.add_csv_file(str(csv_path))
        if yaml_path:
            result_set.add_yaml_file(str(yaml_path))

        return result_set

    def _initialize(self, backbone_type: TaxonomicBackbone, marker_type: Marker) -> None:
        """
        Initialize validators and other resources.
        """
        # Determine if we need taxonomy. We defer initialization of the TR because it's costly.
        needs_taxonomy = False

        # Taxonomic validation needs taxonomy for observed/expected resolution & reconciliation
        if self.config.get('validate_taxonomy', True):
            needs_taxonomy = True
            self.logger.info("Will need taxonomy for taxonomic validation")

        # Structural validation of protein-coding markers needs taxonomy to infer translation table
        do_structural = self.config.get('validate_structure', True)
        if do_structural:
            if marker_type in [Marker.COI_5P, Marker.MATK, Marker.RBCL]:
                needs_taxonomy = True
                self.logger.info("Will need taxonomy for protein-coding marker translation")

        # Setup taxonomy resolver, parser, and possibly taxonomic validator
        if needs_taxonomy:
            self.logger.info("Initializing taxonomy. This may take a while...")
            self.taxonomy_resolver.setup_taxonomy(backbone_type)
            self.logger.debug("Taxonomy initialization complete")

            # Create taxonomy parser
            if backbone_type == TaxonomicBackbone.DWC:
                self.taxonomy_parser = NSRTaxonomyParser(self.config, self.taxonomy_resolver)
            elif backbone_type == TaxonomicBackbone.BOLD:
                self.taxonomy_parser = BOLDTaxonomyParser(self.config, self.taxonomy_resolver)

            # Create taxonomic validator
            if self.config.get('validate_taxonomy', True):
                self.taxonomic_validator = TaxonomicValidator(
                    self.config,
                    self.taxonomy_resolver
                )

        # Create structural validator if needed
        if do_structural:
            if marker_type in [Marker.COI_5P, Marker.MATK, Marker.RBCL]:
                hmm_dir = self.config.get('hmm_profile_dir')
                self.structural_validator = ProteinCodingValidator(self.config, hmm_dir)
            else:
                self.structural_validator = NonCodingValidator(self.config)

    def _parse_input(self, file_path: Path) -> Iterator[SeqRecord]:
        """
        Parse input file in FASTA or TSV format.

        :param file_path: Path to input file
        :return: Iterator of SeqRecord objects
        :raises ValueError: If file format not supported
        """
        suffix = file_path.suffix.lower()

        if suffix in ['.fa', '.fasta', '.fna']:
            self.logger.info(f"Parsing FASTA file: {file_path}")
            yield from SeqIO.parse(file_path, 'fasta')

        elif suffix in ['.tsv', '.txt']:
            self.logger.info(f"Parsing TSV file: {file_path}")
            with open(file_path) as handle:
                yield from SeqIO.parse(handle, 'bcdm-tsv')

        else:
            raise ValueError(f"Unsupported file format: {suffix}")

    def _validate_record(self, record: SeqRecord, dataset: str, backbone_type: TaxonomicBackbone, marker_type: Marker) -> DNAAnalysisResult:
        """
        Validate a single sequence record.

        :param record: The sequence record to validate
        :param dataset: Dataset identifier (e.g., source file name)
        :param backbone_type: Taxonomic backbone to use for validation
        :param marker_type: Marker type to use for validation
        :return: DNAAnalysisResult object containing validation results
        """
        # Instantiate result object, annotate target marker at sequence level (CSC) or config level, set validation level
        result = DNAAnalysisResult(record.id, dataset)
        marker_code = record.annotations.get('bcdm_fields', {}).get('marker_code') # may be in CSC/BCDM
        if marker_code is not None:
            marker_type = Marker(marker_code)
        result.add_ancillary('marker_code', marker_type.value)
        result.level = self.config.get('validation_rank')

        # Enrich taxonomic data and translation table
        self.taxonomy_parser.enrich_result(record, result, backbone_type)
        if result.error:
            return result

        # Perform validations
        if self.structural_validator:
            self.structural_validator.validate_sequence(record, result)
        if self.taxonomic_validator and not result.error:
            self.taxonomic_validator.validate_taxonomy(record, result)

        return result

    def write_results(self, results: DNAAnalysisResultSet,
                     output_path: Path,
                     write_valid_fasta: bool = False,
                     fasta_path: Optional[Path] = None) -> None:
        """
        Write validation results.

        :param results: Validation result set
        :param output_path: Path for TSV output
        :param write_valid_fasta: Whether to output valid sequences as FASTA
        :param fasta_path: Optional path for FASTA output
        """
        # Write TSV results
        self.logger.info(f"Writing results to {output_path}")
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(DNAAnalysisResult.result_fields())
            for result in results.results:
                writer.writerow(result.get_values())

        # Write valid sequences if requested
        if write_valid_fasta:
            self._write_valid_sequences(results, fasta_path)

    def _write_valid_sequences(self, results: DNAAnalysisResultSet,
                             fasta_path: Optional[Path]) -> None:
        """
        Write valid sequences to FASTA.

        :param results: Set of validation results
        :param fasta_path: Optional path for FASTA output
        """
        valid_sequences = []
        for result in results.results:
            if result.passes_all_checks():
                if hasattr(result, 'sequence'):
                    valid_sequences.append(result.sequence)
                else:
                    self.logger.warning(f"No sequence for {result.sequence_id}")

        if valid_sequences:
            output_path = fasta_path or Path('valid_sequences.fasta')
            self.logger.info(f"Writing {len(valid_sequences)} valid sequences to {output_path}")
            SeqIO.write(valid_sequences, output_path, "fasta")
