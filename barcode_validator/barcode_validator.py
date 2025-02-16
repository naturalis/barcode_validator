import tarfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from nbitk.Phylo.NCBITaxdmp import Parser as NCBIParser
from nbitk.Phylo.BOLDXLSXIO import Parser as BOLDParser
from nbitk.Phylo.DwCATaxonomyIO import Parser as DwCParser
from barcode_validator.dna_analysis_result import DNAAnalysisResult, DNAAnalysisResultSet
from barcode_validator.protein_coding_validator import ProteinCodingValidator
from barcode_validator.non_coding_validator import NonCodingValidator
from barcode_validator.taxonomic_validator import TaxonomicValidator, TaxonomicBackbone
from barcode_validator.taxonomy_resolver import Marker, TaxonomyResolver


class BarcodeValidator:
    """
    Main validator class for DNA barcodes.

    This class orchestrates the validation of DNA barcodes by combining structural
    and taxonomic validation. It can handle both FASTA and tabular input formats,
    supports different taxonomic backbones, and can be configured to perform
    structural validation, taxonomic validation, or both.

    Examples:
        >>> from nbitk.config import Config
        >>> config = Config()
        >>> config.load_config('/path/to/config.yaml')
        >>> validator = BarcodeValidator(config)
        >>> validator.initialize()
        >>> results = validator.validate_fasta('sequences.fasta')

    :param config: Configuration object containing validation parameters
    """

    def __init__(self, config: Config, taxonomy_resolver: TaxonomyResolver):
        """
        Initialize the barcode validator.

        :param config: Configuration object containing validation parameters
        :param taxonomy_resolver: TaxonomyResolver instance for handling taxonomic operations
        """
        self.config = config
        self.logger = get_formatted_logger(self.__class__.__name__, config)
        self.taxonomy_resolver = taxonomy_resolver
        self.structural_validator = None
        self.taxonomic_validator = None

    def initialize(self) -> None:
        """
        Initialize validator components.

        Creates appropriate validator instances based on configuration.
        Must be called before validation can be performed.
        """
        # Create structural validator based on marker type
        marker = Marker(self.config.get('marker', 'COI-5P'))
        hmm_dir = Path(self.config.get('hmm_profile_dir'))

        if marker in [Marker.COI_5P, Marker.MATK, Marker.RBCL]:
            self.structural_validator = ProteinCodingValidator(self.config, hmm_dir)
        else:
            self.structural_validator = NonCodingValidator(self.config)

        # Create taxonomic validator if needed
        if self.config.get('validate_taxonomy', True):
            self.taxonomic_validator = TaxonomicValidator(
                self.config,
                self.taxonomy_resolver
            )

    def validate_fasta(self, fasta_file: str) -> DNAAnalysisResultSet:
        """
        Validate sequences from a FASTA file.

        :param fasta_file: Path to FASTA file
        :return: Set of validation results
        """
        results = []
        for record in SeqIO.parse(fasta_file, 'fasta'):
            result = DNAAnalysisResult(record.id, fasta_file)
            self.validate_record(record, result)
            results.append(result)
        return DNAAnalysisResultSet(results)

    def validate_table(self, table_file: str) -> DNAAnalysisResultSet:
        """
        Validate sequences from a tabular file.

        :param table_file: Path to tabular file
        :return: Set of validation results
        """
        results = []
        for record in SeqIO.parse(table_file, 'bcdm-tsv'):
            result = DNAAnalysisResult(record.id, table_file)
            self.validate_record(record, result)
            results.append(result)
        return DNAAnalysisResultSet(results)

    def validate_record(self, record: SeqRecord, result: DNAAnalysisResult) -> None:
        """
        Validate a single sequence record by orchestrating taxonomic resolution
        and delegating structural and taxonomic validation to their respective validators.

        The validators update the result object directly; no overall validity is determined here.

        :param record: The sequence record to validate.
        :param result: Result object to store validation outcomes.
        """
        # Extract identification and rank from bcdm_fields
        bcdm_fields = record.annotations.get('bcdm_fields', {})
        identification = bcdm_fields.get('identification')
        rank = bcdm_fields.get('rank', 'null')

        if not identification:
            result.error = "No taxonomic identification provided"
            return

        # Resolve backbone taxonomy and get validation taxa
        backbone_taxon = self.taxonomy_resolver.resolve_backbone(identification, rank)
        if not backbone_taxon:
            result.error = f"Could not resolve taxon in backbone: {identification}"
            return

        validation_level = self.config.get('level', 'family')
        backbone_validation, ncbi_validation = self.taxonomy_resolver.get_validation_taxon(
            backbone_taxon, validation_level
        )
        if not backbone_validation:
            result.error = f"Could not find {validation_level} rank in backbone for {identification}"
            return
        if not ncbi_validation:
            result.error = f"Could not map {validation_level} {backbone_validation.name} to NCBI taxonomy"
            return

        # Store expected taxon and level in the result
        result.level = validation_level
        result.exp_taxon = backbone_validation

        # Configure structural validation (translation table and constraint)
        if self.structural_validator:
            marker = bcdm_fields.get('marker_code', 'COI-5P')
            try:
                marker_enum = Marker(marker)
            except ValueError:
                self.logger.warning(f"Unknown marker {marker}, using COI-5P")
                marker_enum = Marker.COI_5P

            trans_table = self.taxonomy_resolver.get_translation_table(marker_enum, ncbi_validation)
            self.config.set('translation_table', trans_table)
            constraint_level = self.config.get('constrain', 'class')
            constraint_id = self.taxonomy_resolver.get_constraint_taxon(ncbi_validation, constraint_level)
            self.config.set('constraint_taxid', constraint_id)

        # Delegate validations to the respective validators.
        if self.structural_validator:
            self.structural_validator.validate_sequence(record, result)
        if self.taxonomic_validator:
            self.taxonomic_validator.validate_taxonomy(record, result)

    def _load_taxonomy_trees(self) -> None:
        """Load NCBI and backbone taxonomy trees."""
        # Load NCBI taxonomy
        ncbi_tax_file = self.config.get('ncbi_taxonomy')
        if ncbi_tax_file:
            self.logger.info(f"Loading NCBI taxonomy from {ncbi_tax_file}")
            self.ncbi_tree = NCBIParser(tarfile.open(ncbi_tax_file, "r:gz")).parse()

        # Load appropriate backbone taxonomy
        backbone_type = TaxonomicBackbone(self.config.get('taxonomic_backbone', 'bold'))
        if backbone_type == TaxonomicBackbone.BOLD:
            bold_file = self.config.get('bold_sheet_file')
            if bold_file:
                self.logger.info(f"Loading BOLD taxonomy from {bold_file}")
                with open(bold_file, 'rb') as f:
                    self.backbone_tree = BOLDParser(f).parse()
        else:  # DarwinCore
            dwc_file = self.config.get('dwc_archive')
            if dwc_file:
                self.logger.info(f"Loading DarwinCore taxonomy from {dwc_file}")
                self.backbone_tree = DwCParser(dwc_file).parse()
