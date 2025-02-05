import argparse
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from barcode_validator.barcode_validator import BarcodeValidator
from barcode_validator.result import DNAAnalysisResult, DNAAnalysisResultSet
from barcode_validator.sequence_handler import SequenceHandler

def main(fasta_file_path, logger, config):

    # Process the FASTA file
    sh = SequenceHandler(config)
    validator = BarcodeValidator(config)
    results = []
    logger.info(f"Starting analysis for file: {fasta_file_path}")
    for record, json_config in sh.parse_fasta(fasta_file_path):
        scoped_config = config.local_clone(json_config)
        result = DNAAnalysisResult(record.id, fasta_file_path)
        validator.validate_sequence_quality(scoped_config, record, result)
        results.append(result)

    # Print results
    print(DNAAnalysisResultSet(results))


if __name__ == "__main__":
    # Process command line arguments
    parser = argparse.ArgumentParser(description="Analyze DNA sequences from a FASTA file.")
    parser.add_argument("-f", "--fasta_file", required=True, help="Path to the input FASTA file")
    parser.add_argument("-c", "--config_file", required=True, help="Path to the configuration YAML file")
    parser.add_argument("-v", "--verbosity", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help="Set the logging verbosity (default: WARNING)")
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
    main(args.fasta_file, main_logger, main_config)