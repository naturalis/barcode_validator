import argparse
import os
import sys
import logging
import tempfile
import requests

from nbitk.config import Config
from nbitk.logger import get_formatted_logger

from barcode_validator.barcode_io import BarcodeIO
from barcode_validator.barcode_validator import BarcodeValidator
from barcode_validator.dna_analysis_result import DNAAnalysisResultSet, DNAAnalysisResult


class BarcodeValidatorCLI:
    def __init__(self) -> None:
        """
        Initialize the BarcodeValidatorCLI.
        Stores command line arguments, configuration, and logger for later use.
        """
        self.args: argparse.Namespace = self.parse_args()
        self.config: Config = self.load_and_override_config()
        self.logger: logging.Logger = get_formatted_logger(__name__, self.config)
        self.logger.info("Starting DNA Barcode Validation CLI")
        self.setup_taxonomy()

    def download_ncbi_dump(self, destination: str) -> str:
        """
        Download NCBI taxdump.

        :param destination: Location to store the downloaded file.
        :return: Location of the downloaded file.
        """
        url = self.config.get("ncbi_taxonomy_url")
        if not url:
            self.logger.error("NCBI taxonomy URL is not set in the configuration.")
            sys.exit(1)

        self.logger.info(f"Downloading NCBI taxdump from {url}...")

        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()  # Raise an error for bad HTTP responses (4xx or 5xx)

            with open(destination, "wb") as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)

            self.logger.info(f"Download completed: {destination}")
            return destination

        except requests.exceptions.RequestException as e:
            self.logger.error(f"Failed to download NCBI taxdump: {e}")
            sys.exit(1)

    def download_nsr_dump(self, destination: str) -> str:
        """
        Download NSR dump.

        :param destination: Location to store the downloaded file.
        :return: Location of the downloaded file.
        """
        url = self.config.get("dwc_archive_url")
        if not url:
            self.logger.error("NSR taxonomy URL is not set in the configuration.")
            sys.exit(1)

        self.logger.info(f"Downloading NSR dump from {url}...")

        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()  # Raise an error for bad HTTP responses (4xx or 5xx)

            with open(destination, "wb") as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)

            self.logger.info(f"Download completed: {destination}")
            return destination

        except requests.exceptions.RequestException as e:
            self.logger.error(f"Failed to download NSR dump: {e}")
            sys.exit(1)

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
        # Input options
        parser.add_argument("--input-file", required=True,
                            help="Path to the input sequence file (FASTA or TSV).")
        parser.add_argument("--tsv-file",
                            help="Optional TSV file with additional metadata (if using FASTA input).")
        parser.add_argument("--csv-file",
                            help="Optional CSV file with record-level analytics to merge as ancillary data.")
        parser.add_argument("--yaml-file",
                            help="Optional YAML file with analysis-level configuration to merge into every result.")
        # Taxonomy dump options
        parser.add_argument("--ncbi-dump",
                            help="Path to the NCBI taxdump (tar.gz).")
        parser.add_argument("--nsr-dump",
                            help="Path to the NSR dump (DarwinCore .zip).")
        parser.add_argument("--bold-dump",
                            help="Path to the BOLD dump (XLSX).")
        parser.add_argument("--download-ncbi", action="store_true",
                            help="Download the NCBI taxdump if not available locally.")
        parser.add_argument("--download-nsr", action="store_true",
                            help="Download the NSR dump if not available locally.")
        # Validation and output options
        parser.add_argument("--mode", choices=["structural", "taxonomic", "both"], default="both",
                            help="Validation mode: structural, taxonomic, or both (default: both).")
        parser.add_argument("--emit-valid-fasta", action="store_true",
                            help="If set, emit only valid sequences as FASTA (in addition to TSV output).")
        parser.add_argument("--output-fasta",
                            help="Path to output FASTA file for valid sequences (if not specified, valid FASTA is printed to stdout).")
        # Logging and configuration
        parser.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                            help="Set logging verbosity (overrides config log_level).")
        parser.add_argument("--config-file", required=True,
                            help="Path to the configuration YAML file (e.g., config.yml).")
        return parser.parse_args()

    def load_and_override_config(self) -> Config:
        """
        Load and override configuration from command line arguments.

        :return: Config object.
        """
        config = Config()
        try:
            config.load_config(self.args.config_file)
        except Exception as e:
            sys.exit(f"Error loading config file: {e}")
        if self.args.ncbi_dump:
            config.set("ncbi_taxonomy", self.args.ncbi_dump)
        if self.args.nsr_dump:
            config.set("dwc_archive", self.args.nsr_dump)
        if self.args.bold_dump:
            config.set("bold_sheet_file", self.args.bold_dump)
        if self.args.log_level:
            config.set("log_level", self.args.log_level)
        return config

    def setup_taxonomy(self) -> None:
        """
        Setup taxonomy dumps using the command line arguments and update configuration.
        """
        if self.args.download_ncbi:
            ncbi_path = self.config.get("ncbi_taxonomy")
            downloaded_path = self.download_ncbi_dump(ncbi_path)
            self.config.set("ncbi_taxonomy", downloaded_path)
        if self.args.download_nsr:
            nsr_path = self.config.get("dwc_archive")
            downloaded_path = self.download_nsr_dump(nsr_path)
            self.config.set("dwc_archive", downloaded_path)

    def perform_validation(self) -> DNAAnalysisResultSet:
        """
        Perform validation on input sequences.

        :return: Set of validation results
        """
        # Initialize validator
        validator = BarcodeValidator(self.config)
        validator.initialize()

        # Setup I/O handler
        io_handler = BarcodeIO(self.config)
        results = []

        # Process each record
        for record in io_handler.parse_input(self.args.input_file):
            result = DNAAnalysisResult(record.id, self.args.input_file)
            validator.validate_record(record, result)
            results.append(result)

        return DNAAnalysisResultSet(results)

    def merge_ancillary_data(self, result_set: DNAAnalysisResultSet) -> DNAAnalysisResultSet:
        """
        Merge ancillary data from CSV and YAML files into the result set.

        :param result_set: DNAAnalysisResultSet object.
        :return: Updated DNAAnalysisResultSet object.
        """
        if self.args.csv_file:
            self.logger.info(f"Merging CSV analytics from: {self.args.csv_file}")
            result_set.add_csv_file(self.args.csv_file)
        if self.args.yaml_file:
            self.logger.info(f"Merging YAML configuration from: {self.args.yaml_file}")
            result_set.add_yaml_file(self.args.yaml_file)
        return result_set

    @staticmethod
    def output_tsv_results(result_set: DNAAnalysisResultSet) -> None:
        """
        Output the TSV report to stdout.
        """
        print(result_set)

    def emit_valid_fasta(self, result_set: DNAAnalysisResultSet) -> None:
        """
        Emit only valid sequences as FASTA.
        """
        if self.args.emit_valid_fasta:
            self.logger.info("Emitting valid sequences as FASTA")
            valid_sequences = []
            for result in result_set.results:
                if result.passes_all_checks():
                    if hasattr(result, 'sequence'):
                        valid_sequences.append(result.sequence)
                    else:
                        self.logger.warning(f"No sequence stored for result {result.sequence_id}")
            io_handler = BarcodeIO(self.config)
            if self.args.output_fasta:
                output_path = self.args.output_fasta
                io_handler.write_filtered_fasta(result_set, output_path, valid_only=True)
                self.logger.info(f"Valid sequences written to: {output_path}")
            else:
                with tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta", delete=False) as tmp_file:
                    tmp_path = tmp_file.name
                io_handler.write_filtered_fasta(result_set, tmp_path, valid_only=True)
                self.logger.info(f"Valid sequences written to temporary file: {tmp_path}")
                with open(tmp_path, "r") as f:
                    sys.stdout.write(f.read())
                os.unlink(tmp_path)

    def run(self) -> None:
        """
        Run the CLI tool. This method orchestrates the entire validation process:
        parsing arguments, loading configuration, setting up taxonomy, performing validation,
        merging ancillary data, and outputting results.
        """
        result_set = self.perform_validation()
        result_set = self.merge_ancillary_data(result_set)
        self.output_tsv_results(result_set)
        self.emit_valid_fasta(result_set)