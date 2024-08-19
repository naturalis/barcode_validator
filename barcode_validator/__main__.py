import argparse
import yaml
from barcode_validator.config import Config
from barcode_validator.alignment import *
from barcode_validator.taxonomy import *
from barcode_validator.result import DNAAnalysisResult


def get_tip_by_processid(process_id, tree):
    for tip in tree.get_terminals():
        if tip.guids['processid'] == process_id:
            return tip
    return None


def main(fasta_file_path, bold_sheet):
    logging.info(f"Starting analysis for file: {fasta_file_path}")

    # Create taxonomy trees
    ncbi_tree = read_ncbi_taxonomy(config.get('ncbi_taxonomy'))
    bold_tree = read_bold_taxonomy(bold_sheet)

    # Open the FASTA file and process each record
    with open(fasta_file_path, 'r') as file:
        for process_id, record in parse_fasta(file):

            # Instantiate result object with process ID
            logging.info(f'Processing FASTA record {process_id}')
            result = DNAAnalysisResult(process_id)

            # Lookup species name from process_id
            tip = get_tip_by_processid(process_id, bold_tree)
            result.species = tip.name
            logging.info(f"Species: {result.species}")

            # Lookup expected and observed higher taxon
            for node in bold_tree.root.get_path(tip):
                if node.taxonomic_rank == config.get('level'):
                    result.exp_taxon = node.name
                    break
            result.obs_taxon = run_seqid(record, ncbi_tree)

            # Compute sequence quality metrics
            aligned_sequence = align_to_hmm(record)
            logging.debug(f"Sequence aligned to HMM: {aligned_sequence.seq}")
            amino_acid_sequence = translate_sequence(aligned_sequence)
            result.stop_codons = get_stop_codons(amino_acid_sequence)
            result.seq_length = marker_seqlength(aligned_sequence)
            result.ambiguities = num_ambiguous(record)

            # Stringify the result
            print(result)

    logging.info("Analysis completed")


if __name__ == "__main__":

    # Process command line arguments
    parser = argparse.ArgumentParser(description="Analyze DNA sequences from a FASTA file.")
    parser.add_argument("-f", "--fasta_file", help="Path to the input FASTA file")
    parser.add_argument("-c", "--config_file", help="Path to the configuration YAML file")
    parser.add_argument("-b", "--bold_sheet", help="BOLD XLSX spreadsheet with Lab Sheet and Taxonomy tabs")
    parser.add_argument("-v", "--verbosity", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help="Set the logging verbosity (default: INFO)")
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

    # Run the main analysis
    main(args.fasta_file, args.bold_sheet)
