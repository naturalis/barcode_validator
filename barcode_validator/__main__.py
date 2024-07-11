import argparse
import yaml
from config import Config
from barcode_validator.alignment import *
from barcode_validator.taxonomy import *
from barcode_validator.persistence import *
from barcode_validator.result import DNAAnalysisResult


def main(fasta_file_path, bold_sheet, actions):
    logging.info(f"Starting analysis for file: {fasta_file_path}")

    # Pre-fetch process ID => taxonomy mapping
    if 'id' in actions:
        if bold_sheet:
            mapping = fetch_bold_taxonomy(bold_sheet)
        else:
            raise ValueError("ID action requested but no BOLD spreadsheet provided")

    # Start reading the FASTA file
    with open(fasta_file_path, 'r') as file:
        record = next(SeqIO.parse(file, 'fasta'))
        process_id = record.id.split('_')[0]
        logging.info(f'Processing FASTA record {process_id}')

        # Instantiate result object, unalign sequence
        result = DNAAnalysisResult(process_id)
        record = unalign_sequence(record)

        # Do the ID check
        if 'id' in actions:
            level = config.get('level')
            result.exp_taxon = mapping[process_id][level]
            result.obs_taxon = run_seqid(record)
            result.species = mapping[process_id]['species']

        # Check if sequence has no early stop codons
        if 'stops' in actions:
            aligned_sequence = align_to_hmm(record)
            logging.debug(f"Sequence aligned to HMM: {aligned_sequence.seq}")
            amino_acid_sequence = translate_sequence(aligned_sequence)
            result.stop_codons = get_stop_codons(amino_acid_sequence)
            result.seq_length = marker_seqlength(aligned_sequence)

        # Persist the results to a database
        if 'persist' in actions:
            persist_result(result)

        logging.debug(result)

    logging.info("Analysis completed")


if __name__ == "__main__":

    # Process command line arguments
    parser = argparse.ArgumentParser(description="Analyze DNA sequences from a FASTA file.")
    parser.add_argument("-f", "--fasta_file", help="Path to the input FASTA file")
    parser.add_argument("-c", "--config_file", help="Path to the configuration YAML file")
    parser.add_argument("-b", "--bold_sheet", help="BOLD XLSX spreadsheet with Lab Sheet and Taxonomy tabs")
    parser.add_argument("-a", "--action", help="Comma-separated list of actions, e.g. id,stops,persist")
    parser.add_argument("-v", "--verbosity", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO', help="Set the logging verbosity (default: INFO)")
    args = parser.parse_args()

    # Setup configuration
    config = Config()
    try:
        config.load_config(args.config_file)
    except (FileNotFoundError, yaml.YAMLError) as e:
        print(f"Error loading configuration: {e}")
        exit(1)

    # Setup logging
    try:
        config.setup_logging(args.verbosity)
    except ValueError as e:
        print(f"Error setting up logging: {e}")
        exit(1)

    # Parse the actions
    actions = args.action.split(',')

    # Run the main analysis
    main(args.fasta_file, args.bold_sheet, actions)
