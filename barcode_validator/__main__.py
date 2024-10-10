import argparse
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from barcode_validator.core import BarcodeValidator
from barcode_validator.result import DNAAnalysisResult, DNAAnalysisResultSet


def main(fasta_file_path, yaml_file_path, analytics_file_path, logger, config):
    # Initialize BarcodeValidator
    validator = BarcodeValidator(config)
    validator.initialize()

    # Print header
    logger.info(f"Starting analysis for file: {fasta_file_path}")

    # Validate the FASTA file
    results = validator.validate_fasta(fasta_file_path, config)
    rs = DNAAnalysisResultSet(results)

    # Add YAML and analytics files
    if yaml_file_path is not None:
        rs.add_yaml_file(yaml_file_path)
    if analytics_file_path is not None:
        rs.add_csv_file(analytics_file_path)

    # Add fasta file path to each result
    for result in rs.results:
        result.add_ancillary('fasta_file', fasta_file_path)

    print(rs)  # print TSV results

    logger.info("Analysis completed")


if __name__ == "__main__":
    # Process command line arguments
    parser = argparse.ArgumentParser(description="Analyze DNA sequences from a FASTA file.")
    parser.add_argument("-f", "--fasta_file", required=True, help="Path to the input FASTA file")
    parser.add_argument("-y", "--yaml_file", required=True, help="Path to the assembly param YAML file")
    parser.add_argument("-a", "--analytics_file", required=True, help="Path to the analytics CSV file")
    parser.add_argument("-c", "--config_file", required=True, help="Path to the configuration YAML file")
    parser.add_argument("-v", "--verbosity", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help="Set the logging verbosity (default: INFO)")
    args = parser.parse_args()

    # Setup logging
    main_config = Config()
    try:
        main_config.load_config(args.config_file)
        main_logger = get_formatted_logger(__name__, main_config)
    except ValueError as e:
        print(f"Error setting up logging: {e}")
        exit(1)

    # Run the main analysis
    main(args.fasta_file, args.yaml_file, args.analytics_file, main_logger, main_config)
