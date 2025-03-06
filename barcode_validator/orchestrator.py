import csv
from pathlib import Path
from typing import Optional, Iterator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from nbitk.SeqIO.BCDM import BCDMIterator
from nbitk.Tools import Blastn, Hmmalign
from .dna_analysis_result import DNAAnalysisResult, DNAAnalysisResultSet
from .resolvers.factory import ResolverFactory
from .idservices.factory import IDServiceFactory
from .validators.factory import StructureValidatorFactory
from .validators.taxonomic import TaxonomicValidator
from .constants import Marker, TaxonomicRank, TaxonomicBackbone

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
        >>> results = orchestrator.validate_file(Path('sequences.fasta'))
        >>> orchestrator.write_results(results, Path('output.tsv'))
    """

    def __init__(self, config: Config):
        """
        Initialize the orchestrator.

        :param config: Configuration object containing validation parameters
        """
        self.config = config
        self.logger = get_formatted_logger(self.__class__.__name__, config)
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
        marker_type = Marker(self.config.get('marker'))
        self._initialize(marker_type)

        # Validate records and create result set
        results = []
        for record in self._parse_input(input_path):
            result = self._validate_record(record, str(input_path), marker_type)
            results.append(result)
        result_set = DNAAnalysisResultSet(results)

        # Add additional data if provided
        if csv_path:
            result_set.add_csv_file(str(csv_path))
        if yaml_path:
            result_set.add_yaml_file(str(yaml_path))

        return result_set

    def _initialize(self, marker_type: Marker) -> None:
        """
        Initialize validators and other resources.
        """

        # Instantiate taxonomic validator if reverse taxonomy is required
        if self.config.get('validate_taxonomy', True):
            self._initialize_tv()

        # Instantiate structural validator subclass if structural validation is required
        if self.config.get('validate_structure', True):
            self._initialize_sv(marker_type)

    def _initialize_sv(self, marker_type) -> None:
        """
        Initialize structural validator.
        :param marker_type: Marker type to use for validation
        """
        sv = StructureValidatorFactory.create_validator(self.config, marker_type)
        self.structural_validator = sv

        # Set up hmmalign if needed
        if sv.requires_hmmalign():
            sv.set_hmm_profile_dir(self.config.get('hmm_profile_dir'))
            sv.set_hmmalign(Hmmalign(self.config))

        # Set up resolver if needed
        if sv.requires_resolver():
            svbb = TaxonomicBackbone(self.config.get('input_taxonomy'))
            svpath = self.config.get(svbb.value + '_file')
            sr = ResolverFactory.create_resolver(self.config, svbb)
            self.logger.info(f"Loading input taxonomy from {svpath}")
            sr.load_tree(Path(svpath))
            sv.set_taxonomy_resolver(sr)

    def _initialize_tv(self) -> None:
        """
        Initialize taxonomic validator.
        """
        self.taxonomic_validator = TaxonomicValidator(self.config)

        # Set up the resolver for taxonomic validation
        if self.taxonomic_validator.requires_resolver():

            # Initialize the TaxonResolver for taxonomic validation
            tvbb = TaxonomicBackbone(self.config.get('reference_taxonomy'))
            tvpath = self.config.get(tvbb.value + '_file')
            tr = ResolverFactory.create_resolver(self.config, tvbb)
            self.logger.info(f"Loading reference taxonomy from {tvpath}")
            tr.load_tree(Path(tvpath))
            self.taxonomic_validator.set_taxonomy_resolver(tr)

            # Initialize an IDService if needed
            if self.taxonomic_validator.requires_idservice():
                ids = IDServiceFactory.create_idservice(self.config, self.config.get('reference_taxonomy'))

                # There is no way right now that the IDService could have a different
                # TaxonResolver than the TaxonomicValidator
                if ids.requires_resolver():
                    ids.set_taxonomy_resolver(tr)

                # If the IDService needs BLASTN, all externally configurable settings
                # are injected here
                if ids.requires_blastn():
                    blastn = Blastn(self.config)
                    blastn.set_db(self.config.get('blast_db'))
                    blastn.set_num_threads(self.config.get('num_threads'))
                    blastn.set_evalue(self.config.get('evalue'))
                    blastn.set_outfmt(self.config.get('outfmt'))
                    blastn.set_max_target_seqs(self.config.get('max_target_seqs'))
                    blastn.set_word_size(self.config.get('word_size'))
                    ids.set_blastn(blastn)
                self.taxonomic_validator.set_idservice(ids)

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

    def _validate_record(self, record: SeqRecord, dataset: str, marker_type: Marker) -> Optional[DNAAnalysisResult]:
        """
        Validate a single sequence record.

        :param record: The sequence record to validate
        :param dataset: Dataset identifier (e.g., source file name)
        :param marker_type: Marker type to use for validation
        :return: DNAAnalysisResult object containing validation results
        """
        # Extract group ID if present
        group_id = None
        if self.config.get('group_id_separator'):
            group_id = record.id.split(self.config.get('group_id_separator'))[0]

        # Instantiate result object, annotate target marker at sequence level (CSC) or config level, set validation level
        result = DNAAnalysisResult(record.id, dataset, group_id=group_id)
        marker_code = record.annotations.get('bcdm_fields', {}).get('marker_code') # may be in CSC/BCDM
        if marker_code is not None:
            marker_type = Marker(marker_code)
        result.add_ancillary('marker_code', marker_type.value)

        # Perform validations
        if self.structural_validator:
            if self.structural_validator.marker != marker_type:
                result.error = f"Marker type mismatch: expected {self.structural_validator.marker_type.value}, got {marker_type.value}"
                return None
            else:
                self.structural_validator.validate(record, result)
        if self.taxonomic_validator and not result.error:
            constraint_rank = TaxonomicRank(self.config.get('constraint_rank', 'class'))
            self.taxonomic_validator.validate(record, result, constraint_rank)

        return result

    def write_results(self, results: DNAAnalysisResultSet,
                     output_format: str = 'tsv', triage: bool = False) -> None:
        """
        Write validation results.

        :param results: Validation result set
        :param output_format: Output format (tsv or fasta, default: tsv)
        :param triage: Perform triage on the result set (default: false)
        """
        # Write TSV results
        self.logger.info(f"Writing results")

        # Triage the results
        if triage:
            results = results.triage()

        # Write the results
        print(results.to_string(output_format))


