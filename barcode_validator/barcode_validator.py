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
from barcode_validator.taxonomy_resolver import Marker


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

    def __init__(self, config: Config):
        """Initialize the barcode validator."""
        self.config = config
        self.logger = get_formatted_logger(self.__class__.__name__, config)
        self.structural_validator = None
        self.taxonomic_validator = None
        self.ncbi_tree = None
        self.backbone_tree = None

    def initialize(self) -> None:
        """
        Initialize validator components.

        Creates appropriate validator instances based on configuration.
        Must be called before validation can be performed.
        """
        # Load taxonomy trees
        self._load_taxonomy_trees()

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
                self.ncbi_tree,
                self.backbone_tree
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
        Validate a single sequence record.

        :param record: The sequence record to validate
        :param result: Result object to store validation outcomes
        """
        # Structural validation
        if self.structural_validator:
            structural_valid, details = self.structural_validator.validate_sequence(record)
            result.add_ancillary('structural_valid', str(structural_valid))
            for key, value in details.items():
                result.add_ancillary(f'structural_{key}', str(value))

        # Taxonomic validation
        if self.taxonomic_validator:
            taxonomic_valid, details = self.taxonomic_validator.validate_taxonomy(record)
            result.add_ancillary('taxonomic_valid', str(taxonomic_valid))
            for key, value in details.items():
                result.add_ancillary(f'taxonomic_{key}', str(value))

            # Store expected taxon for the validation rank
            if details.get('rank_taxon'):
                result.exp_taxon = details['rank_taxon']

        # Overall validation requires both checks to pass if both are enabled
        is_valid = True
        if self.structural_validator:
            is_valid = is_valid and structural_valid
        if self.taxonomic_validator:
            is_valid = is_valid and taxonomic_valid
        result.add_ancillary('is_valid', str(is_valid))

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
