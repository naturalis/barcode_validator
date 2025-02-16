import argparse
import os
import sys
import logging
import tempfile
from pathlib import Path

from nbitk.config import Config
from nbitk.logger import get_formatted_logger

from .orchestrator import ValidationOrchestrator
from .dna_analysis_result import DNAAnalysisResultSet, DNAAnalysisResult


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
    6. Output the result set as TSV to stdout
    7. Optionally emit only valid sequences as FASTA to stdout
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
            description="""DNA Barcode Validator CLI.
Ingests FASTA/TSV sequence data, merges CSV analytics and YAML config,
performs validations, and outputs results as TSV. Optionally emits valid sequences as FASTA.""",
            formatter_class=argparse.RawTextHelpFormatter
        )

        # TODO: ensure that these are valid flags as object properties (probably need snake_case)
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
                            help="Path to the NSR dump (DarwinCore .zip) or the the BOLD dump (XLSX).")
        parser.add_argument("--exp_taxonomy_type", choices=["nsr", "bold"], default="nsr",
                            help="Type of the expected taxonomy dump (nsr or bold).")

        # Validation and output options
        parser.add_argument("--mode", choices=["structural", "taxonomic", "both"], default="both",
                            help="Validation mode: structural, taxonomic, or both (default: both).")
        parser.add_argument("--emit_valid_fasta", action="store_true",
                            help="If set, emit only valid sequences as FASTA (in addition to TSV output).")
        parser.add_argument("--output_fasta",
                            help="Path to output FASTA file for valid sequences (if not specified, valid FASTA is printed to stdout).")
        parser.add_argument("--output_file", required=True,  # This is used in run() but not defined
                            help="Output TSV file for validation results")

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
            config.set("ncbi_taxonomy", self.args.reflib_taxonomy)
        if self.args.exp_taxonomy_type == 'nsr':
            config.set("dwc_archive", self.args.exp_taxonomy)
        elif self.args.exp_taxonomy_type == 'bold':
            config.set("bold_sheet_file", self.args.exp_taxonomy)
        if self.args.log_level:
            config.set("log_level", self.args.log_level)
        return config

    def run(self):
        """Run the validation process."""
        try:
            # Initialize orchestrator
            orchestrator = ValidationOrchestrator(self.config)
            orchestrator.initialize()

            # Validate sequences
            results = orchestrator.validate_file(
                Path(self.args.input_file),
                Path(self.args.csv_file) if self.args.csv_file else None,
                Path(self.args.yaml_file) if self.args.yaml_file else None
            )

            # Write results
            orchestrator.write_results(
                results,
                Path(self.args.output_file),
                self.args.emit_valid_fasta,
                Path(self.args.output_fasta) if self.args.output_fasta else None
            )

        except Exception as e:
            self.logger.error(f"Validation failed: {e}")
            sys.exit(1)