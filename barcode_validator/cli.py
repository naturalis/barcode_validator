import argparse
import sys
import logging
from pathlib import Path
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from .orchestrator import ValidationOrchestrator
from .dna_analysis_result import DNAAnalysisResultSet


class BarcodeValidatorCLI:
    """
    Command Line Interface for DNA Barcode Validator.

    This class exposes the functionality of the barcode_validator package to the command line.
    The overall program flow is as follows:
    1. Parse command line arguments
    2. Load configuration from YAML and override with command line arguments
    3. Instantiate the orchestrator and initialize it with the loaded configuration
    4. Run the orchestrator and obtain a result set object
    5. Merge ancillary data from CSV and YAML files into the result set
    6. Optionally perform triage on the result set, possibly at the sequence group level
    7. Output the result set as TSV or FASTA to stdout
    """
    def __init__(self) -> None:
        """
        Initialize the BarcodeValidatorCLI.
        """
        self.args: argparse.Namespace = self.parse_args()
        self.config: Config = self.load_and_override_config()
        self.logger: logging.Logger = get_formatted_logger(__name__, self.config)
        self.logger.info("Starting DNA Barcode Validation CLI")

    @staticmethod
    def parse_args() -> argparse.Namespace:
        """
        Parse command line arguments.

        :return: Parsed command line arguments.
        """
        parser = argparse.ArgumentParser(
            description="""DNA Barcode Validator
This tool validates DNA barcode sequence data in FASTA or BCDM/TSV format and validates it. Validation consists of
either or both of the following steps:

1. Structural validation: Check the sequence for minimum length, number of ambiguous bases, and, in the case of
   protein-coding sequences, for the presence of unexpected stop codons.
2. Taxonomic validation: Compare the sequence against a reference library and an expected taxonomy to determine
   whether the sequence is taxonomically valid in the sense that the provided taxon is observed among the results
   of the reference library search. This check can be performed at different taxonomic ranks.

Optionally, the tool performs a triage on the result set that filters out invalid sequences based on the validation
results. If input sequences are grouped (e.g., by sample), the triage can be performed at the group level, in which
case the longest valid sequence in each group is retained.

Optionally, the tool can merge CSV analytics and YAML configuration into each result, which are then returned in the
output TSV.
""",
            formatter_class=argparse.RawTextHelpFormatter
        )

        # Input file options.
        parser.add_argument("--input_file", required=True,
                            help="Path to the input sequence file (FASTA or BCDM/TSV).")
        parser.add_argument("--csv_file",
                            help="Optional CSV file with record-level analytics to merge as ancillary data.")
        parser.add_argument("--yaml_file",
                            help="Optional YAML file with analysis-level configuration to merge into every result.")

        # Taxonomy dump options. The locations are updated the config, which is passed into the system.
        parser.add_argument("--reflib_taxonomy",
                            help="Path to the taxonomy of the reference library, i.e. NCBI taxdump (tar.gz).")
        parser.add_argument("--exp_taxonomy",
                            help="Path to the expected taxonomy, e.g. an NSR dump or a BOLD dump spreadsheet.")
        parser.add_argument("--exp_taxonomy_type", choices=["nsr", "bold"], default="nsr",
                            help="Type of the expected taxonomy dump (nsr or bold).")

        # Validation and output options
        parser.add_argument("--mode", choices=["structural", "taxonomic", "both"], default="both",
                            help="Validation mode: structural, taxonomic, or both (default: both).")
        parser.add_argument("--output_format", choices=["tsv", "fasta"], default="tsv",
                            help="Whether to output results as TSV or FASTA (default: tsv).")
        parser.add_argument("--triage", action="store_true", help="Perform triage on the result set.")

        # Logging and configuration. The log level is updated in the config, and the config is passed into the system.
        parser.add_argument("--log_level", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                            help="Set logging verbosity (overrides config log_level).")
        parser.add_argument("--config", required=True,
                            help="Path to the configuration YAML file (e.g., config.yml).")
        return parser.parse_args()

    def load_and_override_config(self) -> Config:
        """
        Load and override configuration from command line arguments.

        :return: Config object.
        """
        config = Config()
        try:
            config.load_config(self.args.config)
        except Exception as e:
            sys.exit(f"Error loading config file: {e}")
        if self.args.reflib_taxonomy:
            config.set("reflib_taxonomy", self.args.reflib_taxonomy)
        if self.args.exp_taxonomy_type == 'nsr':
            config.set("exp_taxonomy", self.args.exp_taxonomy)
            config.set("exp_taxonomy_type", self.args.exp_taxonomy_type)
        elif self.args.exp_taxonomy_type == 'bold':
            config.set("exp_taxonomy", self.args.exp_taxonomy)
            config.set("exp_taxonomy_type", self.args.exp_taxonomy_type)
        if self.args.log_level:
            config.set("log_level", self.args.log_level)
        return config

    def run(self) -> DNAAnalysisResultSet:
        """Run the validation process."""
        try:
            # Initialize orchestrator
            orchestrator = ValidationOrchestrator(self.config)

            # Validate sequences
            results = orchestrator.validate_file(
                Path(self.args.input_file),
                Path(self.args.csv_file) if self.args.csv_file else None,
                Path(self.args.yaml_file) if self.args.yaml_file else None
            )

            # Write results
            orchestrator.write_results(
                results,
                self.args.output_format,
                self.args.triage
            )
            return results

        except Exception as e:
            self.logger.error(f"Validation failed: {e}")
            sys.exit(1)