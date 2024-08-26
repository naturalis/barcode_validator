import argparse
import logging
from barcode_validator.config import Config
from barcode_validator.core import BarcodeValidator
import barcode_validator.result


def main(fasta_file_path):

    # Initialize BarcodeValidator
    validator = BarcodeValidator(config)
    validator.initialize()

    # Print header
    logging.info(f"Starting analysis for file: {fasta_file_path}")
    header = barcode_validator.result.result_fields()
    header.append('fasta_file')
    print('\t'.join(header))  # print TSV header

    # Validate the FASTA file
    results = validator.validate_fasta(fasta_file_path)

    # Print results
    for result in results:
        values = result.get_values()
        values.append(fasta_file_path)
        print('\t'.join(map(str, values)))

    logging.info("Analysis completed")


if __name__ == "__main__":
    # Process command line arguments
    parser = argparse.ArgumentParser(description="Analyze DNA sequences from a FASTA file.")
    parser.add_argument("-f", "--fasta_file", required=True, help="Path to the input FASTA file")
    parser.add_argument("-c", "--config_file", required=True, help="Path to the configuration YAML file")
    parser.add_argument("-v", "--verbosity", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help="Set the logging verbosity (default: INFO)")
    args = parser.parse_args()

    # Setup logging
    config = Config()
    try:
        config.load_config(args.config_file)
        config.setup_logging(args.verbosity)
    except ValueError as e:
        print(f"Error setting up logging: {e}")
        exit(1)

    # Run the main analysis
    main(args.fasta_file)
