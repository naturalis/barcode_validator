import argparse
import sys
import logging
import traceback
from pathlib import Path
from nbitk.logger import get_formatted_logger

from barcode_validator.constants import ValidationMode
from .orchestrator import ValidationOrchestrator
from .dna_analysis_result import DNAAnalysisResultSet
from barcode_validator.config.schema_config import SchemaConfig


class BarcodeValidatorCLI:
    """
    Command Line Interface for DNA Barcode Validator.

    This class exposes the functionality of the barcode_validator package to the command line.
    The overall program flow is as follows:
    1. Initialize schema-driven configuration
    2. Parse command line arguments using schema-generated parser
    3. Update configuration with command line arguments
    4. Validate required fields
    5. Instantiate the orchestrator and initialize it with the loaded configuration
    6. Run the orchestrator and obtain a result set object
    7. Merge ancillary data from CSV and YAML files into the result set
    8. Optionally perform triage on the result set, possibly at the sequence group level
    9. Output the result set as TSV or FASTA to stdout
    """

    def __init__(self) -> None:
        """Initialize the BarcodeValidatorCLI."""
        self.config: SchemaConfig = SchemaConfig()
        self.args: argparse.Namespace = self.parse_args()
        self.logger: logging.Logger = get_formatted_logger(__name__, self.config)
        self.logger.info("Starting DNA Barcode Validation CLI")

    def create_parser(self) -> argparse.ArgumentParser:
        """
        Create and configure the argument parser using the schema.

        :return: Configured ArgumentParser
        """
        parser = argparse.ArgumentParser(
            description="""DNA Barcode Validator

This tool validates DNA barcode sequence data in FASTA or BCDM/TSV format. Validation consists of
either or both of the following steps:

1. Structural validation: Check the sequence for minimum length, number of ambiguous bases, and, in the case of
   protein-coding sequences, for the presence of unexpected stop codons.

2. Taxonomic validation: Compare the sequence against a reference library and an expected taxonomy to determine
   whether the sequence is taxonomically valid in the sense that the provided taxon is observed among the results
   of the reference library search. This check can be performed at different taxonomic ranks.

Optionally, the tool performs a triage on the result set that filters out invalid sequences based on the validation
results. If input sequences are grouped (e.g., by sample), the triage can be performed at the group level, in which
case the longest valid sequence in each group is retained.

Configuration parameters can be set using command line arguments. For nested parameters, use the format:
--section-name key=value

""",
            formatter_class=argparse.RawTextHelpFormatter
        )

        # Let the schema populate the parser
        self.config.populate_argparse(parser)

        return parser

    def parse_args(self) -> argparse.Namespace:
        """
        Parse command line arguments using schema-generated parser.

        :return: Parsed command line arguments.
        """
        parser = self.create_parser()
        args = parser.parse_args()

        # Update configuration with parsed arguments
        self.config.update_from_args(args)

        # Validate that all required fields are present
        try:
            self.config.validate_required_fields()
        except Exception as e:
            parser.error(f"Configuration validation failed: {e}")

        return args

    def run(self) -> DNAAnalysisResultSet:
        """Run the validation process."""
        try:
            # Initialize orchestrator
            orchestrator = ValidationOrchestrator(self.config)

            # Get file paths from configuration
            input_file = self.config.get('input_file')
            csv_file = self.config.get('csv_file')
            yaml_file = self.config.get('yaml_file')

            # Validate sequences
            results = orchestrator.validate_file(
                Path(input_file),
                Path(csv_file) if csv_file else None,
                Path(yaml_file) if yaml_file else None
            )

            # Write results
            output_fasta = self.config.get('output_fasta')
            output_tsv = self.config.get('output_tsv')
            mode = ValidationMode(self.config.get('mode', 'both'))
            orchestrator.write_results(
                results,
                output_fasta,
                output_tsv,
                mode
            )
            return results

        except Exception as e:
            self.logger.error(f"Validation failed: {e}")
            stack_trace = traceback.format_exc()
            self.logger.error(stack_trace)
            sys.exit(1)


def main():
    """Main entry point for the CLI."""
    cli = BarcodeValidatorCLI()
    cli.run()


if __name__ == "__main__":
    main()