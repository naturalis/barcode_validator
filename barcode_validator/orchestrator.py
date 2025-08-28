import csv
import os
from pathlib import Path
from typing import Optional, Iterator, List
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
from .constants import Marker, TaxonomicRank, TaxonomicBackbone, RefDB, ValidationMode
from .criteria import MarkerCriteriaFactory

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
    lifecycle to ensure consistent data flow and error handling. It is
    therefore solely responsible for translating configuration settings
    from the CLI to internal variables (ideally Enums or other boxed types).
    If you find yourself fiddling with config options in deeper classes
    this should be considered a major red flag that will impede maintainability.
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
        mode = ValidationMode(self.config.get('mode'))
        self._initialize(marker_type, mode)

        # Validate records and create result set
        result_set = DNAAnalysisResultSet([])
        records = self._parse_input(input_path)
        self.logger.info(f"Validating {len(records)} records from {input_path}")
        self._validate_records(records, str(input_path), marker_type, result_set)

        # Add additional data if provided
        if csv_path:
            result_set.add_csv_file(str(csv_path))
        if yaml_path:
            result_set.add_yaml_file(str(yaml_path))

        return result_set

    def _initialize(self, marker_type: Marker, mode: ValidationMode) -> None:
        """
        Initialize validators and other resources.
        """

        # Instantiate taxonomic validator if reverse taxonomy is required
        if mode == ValidationMode.TAXONOMIC or mode == ValidationMode.BOTH:
            self._initialize_tv()

        # Instantiate structural validator subclass if structural validation is required
        if mode == ValidationMode.STRUCTURAL or mode == ValidationMode.BOTH:
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

            # Initialize the TaxonResolver for the input. For every input sequence, this resolver
            # is able to figure out its intended taxon (e.g. from the BOLD process ID) and higher
            # lineage. This is both to compare higher taxa in this lineage against those returned
            # by the ID service, and to infer the translation table.
            sr = self._prepare_resolver(self.config.get('input_resolver.format'),
                                        self.config.get('input_resolver.file'))

            sv.set_taxonomy_resolver(sr)

    def _initialize_tv(self) -> None:
        """
        Initialize taxonomic validator.
        """
        self.taxonomic_validator = TaxonomicValidator(self.config)

        # Set up the resolver for taxonomic validation
        if self.taxonomic_validator.requires_resolver():

            # Initialize the TaxonResolver for the input. For every input sequence, this resolver
            # is able to figure out its intended taxon (e.g. from the BOLD process ID) and higher
            # lineage. This is both to compare higher taxa in this lineage against those returned
            # by the ID service, and to infer the translation table.
            tr = self._prepare_resolver(self.config.get('input_resolver.format'),
                                        self.config.get('input_resolver.file'))
            self.taxonomic_validator.set_taxonomy_resolver(tr)

            # Initialize an IDService if needed
            if self.taxonomic_validator.requires_idservice():
                db = RefDB(self.config.get('taxon_validation.method'))
                ids = IDServiceFactory.create_idservice(self.config, db)
                ids.set_min_identity(self.config.get('taxon_validation.min_identity'))
                ids.set_max_target_seqs(self.config.get('taxon_validation.max_target_seqs'))

                # Initialize the TaxonResolver for the ID Service. Local BLAST needs this to
                # construct the higher lineage for the hits.
                if ids.requires_resolver():
                    rdbr = self._prepare_resolver(self.config.get('reflib_resolver.format'),
                                                  self.config.get('reflib_resolver.file'))
                    ids.set_taxonomy_resolver(rdbr)

                # If the IDService needs BLASTN, all externally configurable settings
                # are injected here
                if ids.requires_blastn():
                    ids.set_blastn(self._configure_blastn())
                self.taxonomic_validator.set_idservice(ids)

    def _prepare_resolver(self, file_format: str, file: Path):
        """
        Prepares a taxon resolver that reads a particular file format
        :param file_format: A value of TaxonomicBackbone
        :param file: Path to a taxonomy dump
        :return:
        """
        tvbb = TaxonomicBackbone(file_format)
        tr = ResolverFactory.create_resolver(self.config, tvbb)
        self.logger.info(f"Loading {file_format} reference taxonomy from {file}")
        tr.load_tree(file)
        return tr

    def _configure_blastn(self):
        """
        Configures the blastn instance with user provided variables
        :return: A configured blastn instance
        """

        # These options are system specific for when users run a local blast database, which needs
        # a location. The system also will have system specific, optimal thread numbers and user
        # settable values for the maximum E and the starting word size.
        blastn = Blastn(self.config)
        blastn.set_db(self.config.get('local_blast.db'))
        blastn.set_num_threads(self.config.get('local_blast.threads'))
        blastn.set_evalue(self.config.get('local_blast.max_evalue'))
        blastn.set_word_size(self.config.get('local_blast.word_size'))

        # Set the BLASTDB environment variable if not already set and if blast_db is in config
        blast_db = self.config.get('local_blast.db')
        if blast_db is not None:

            # Get the directory containing the blast database files
            blast_db_dir = str(Path(blast_db).parent)

            # Check if BLASTDB environment variable is set
            if 'BLASTDB' not in os.environ:
                # Set BLASTDB to the database directory
                os.environ['BLASTDB'] = blast_db_dir
        if blast_db is None:
            self.logger.warning(
                "Command line variable `--local-blast db=<path>` is not set. Local BLAST searches will fail.")

        return blastn

    def _parse_input(self, file_path: Path) -> List[SeqRecord]:
        """
        Parse input file in FASTA or TSV format.

        :param file_path: Path to input file
        :return: List of SeqRecord objects
        :raises ValueError: If file format not supported
        """
        suffix = file_path.suffix.lower()
        records = []

        if suffix in ['.fa', '.fasta', '.fna', '.fas']:
            self.logger.info(f"Parsing FASTA file: {file_path}")
            for record in SeqIO.parse(file_path, 'fasta'):
                records.append(record)

        elif suffix in ['.tsv', '.txt']:
            self.logger.info(f"Parsing TSV file: {file_path}")
            with open(file_path) as handle:
                for record in SeqIO.parse(handle, 'bcdm-tsv'):
                    records.append(record)
        else:
            raise ValueError(f"Unsupported file format: {suffix}")

        return records

    def _validate_records(self, records: List[SeqRecord], dataset: str, marker: Marker, rs: DNAAnalysisResultSet):
        """
        Validate a single sequence record.

        :param records: A list of sequence record to validate
        :param dataset: Dataset identifier (e.g., source file name)
        :param marker: Marker type to use for validation
        :param rs: DNAAnalysisResultSet object to be populated with validation results
        """
        self._populate_resultset(dataset, marker, records, rs)

        # Perform structural validation if requested. This is quick so can be done sequentially.
        if self.structural_validator:
            for result, record in zip(rs.results, records):
                if result.criteria.marker_type !=  self.structural_validator.marker:
                    result.error = f"Marker type mismatch: expected {self.structural_validator.marker_type.value}, got {marker.value}"
                else:
                    self.logger.info(f"Structurally validating record {record.id}")
                    self.structural_validator.validate(record, result)

        # If taxonomic validation is requested and there were no errors thus far, perform the validation.
        # To speed things up, we don't do this on errors or if the sequence was structurally invalid.
        if self.taxonomic_validator:

            # Set up constraint rank
            constraint_rank = TaxonomicRank.CLASS
            if self.config.get('local_blast.extent') is not None:
                constraint_rank = TaxonomicRank(self.config.get('local_blast.extent'))

            # Iterate over records and results, aggregate structurally valid ones into batches
            batch = []
            for result, record in zip(rs.results, records):

                # Add to batch if no problems so far in structural validation
                if self.structural_validator:
                    if not result.error and result.check_seq_quality():
                        batch.append((result,record))
                    else:
                        self.logger.warning(f"Skipping {record.id} for taxonomic validation due to prior errors: {result.error}")

                # Otherwise, if no structural validation done yet, add all records without errors
                elif not self.structural_validator:
                    if not result.error:
                        # Make sure that result.exp_taxon is assigned through tr.enrich_result() because this may not have been done
                        # if we aren't doing taxonomic validation
                        if result.exp_taxon is None:
                            self.taxonomic_validator.taxonomy_resolver.enrich_result(record, result)
                            if result.error:
                                continue
                            result.add_ancillary('nuc', str(record.seq))
                        batch.append((result,record))
                    else:
                        self.logger.warning(f"Skipping {record.id} for taxonomic validation due to prior errors: {result.error}")

            # Process batches
            self.logger.info(f"Going to taxonomically validate {len(batch)} records")
            for i in range(0, len(batch), 100):
                batch_slice = batch[i:i + 100]
                self.logger.info(f"Taxonomically validating batch {i//100 + 1} containing {len(batch_slice)} records")
                self.taxonomic_validator.validate_batch(batch_slice, constraint_rank)

    def _populate_resultset(self, dataset: str, marker: Marker, records: List[SeqRecord], rs: DNAAnalysisResultSet):
        """
        Populate resultset with initial sequence metadata
        :param dataset: Dataset identifier (e.g., source file name)
        :param marker: Marker type to use for validation
        :param records: List of SeqRecord objects
        :param rs: DNAAnalysisResultSet object to be populated with validation results
        """

        for record in records:

            # Parse out group ID if needed
            group_id = None
            if self.config.get('triage_config.group_by_sample'):
                sep = self.config.get('triage_config.group_id_separator')
                group_id = record.id.split(sep)[0]

            # Overwrite marker type at record level if provided
            marker_code = record.annotations.get('bcdm_fields', {}).get('marker_code')  # may be in CSC/BCDM
            if marker_code is not None:
                marker = Marker(marker_code)

            # Instantiate result
            criteria = MarkerCriteriaFactory.get_criteria(marker, self.config)
            result = DNAAnalysisResult(record.id, dataset=dataset, group_id=group_id, criteria=criteria)
            result.add_ancillary('marker_code', str(marker.value))

            # Add to result set
            rs.results.append(result)

    def write_results(self, results: DNAAnalysisResultSet,
                     output_fasta: Path, output_tsv: Path, mode: ValidationMode) -> None:
        """
        Write validation results.

        :param results: Validation result set
        :param output_fasta: Output FASTA file path for valid sequences
        :param output_tsv: Output TSV file path for results
        :param mode: Validation mode (e.g., 'both', 'taxonomic', 'structural')
        """
        # We do the triage first so that the TSV output will contain any applicable descriptive error
        # Write FASTA for valid sequences
        triaged = results.triage(mode, self.config.get('triage_config.group_by_sample', False))
        with open(output_fasta, 'w') as fasta_file:
            fasta_file.write(triaged.to_string('fasta'))

        # Write TSV results
        self.logger.info(f"Writing results")
        with open(output_tsv, 'w', newline='\n') as tsv_file:
            tsv_file.write(results.to_string('tsv'))




