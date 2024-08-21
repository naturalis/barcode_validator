import argparse
import yaml
import logging
import barcode_validator.result
from barcode_validator.config import Config
from barcode_validator.alignment import parse_fasta, align_to_hmm, translate_sequence, get_stop_codons, \
    marker_seqlength, num_ambiguous
from barcode_validator.taxonomy import read_ncbi_taxonomy, read_bold_taxonomy, get_tip_by_processid, run_localblast
from barcode_validator.result import DNAAnalysisResult


def main(fasta_file_path, bold_sheet):
    logging.info(f"Starting analysis for file: {fasta_file_path}")
    print(barcode_validator.result.result_fields())  # print TSV header

    # Create taxonomy trees
    ncbi_tree = read_ncbi_taxonomy(config.get('ncbi_taxonomy'))
    bold_tree = read_bold_taxonomy(bold_sheet)

    # Open the FASTA file and process each record
    for process_id, record in parse_fasta(fasta_file_path):

        # Validate the record. Here we are only printing the result as a TSV row
        # but we can do more with the result object. In particular, we can use the
        # result.calculate_ranks() method to get the barcode and full sequence ranks
        # as well as a list of messages. These we would attach to the pull request
        # as part of the review process.
        result = validate_record(bold_tree, ncbi_tree, process_id, record)
        print(result)

    logging.info("Analysis completed")


def validate_record(bold_tree, ncbi_tree, process_id, record):

    # Instantiate result object with process ID and calculate full sequence stats
    logging.info(f'Processing FASTA record {process_id}')
    result = DNAAnalysisResult(process_id)
    result.full_length = len(record.seq)
    result.full_ambiguities = num_ambiguous(record)

    # Lookup species Taxon from process_id
    result.species = get_tip_by_processid(process_id, bold_tree)
    logging.info(f"Species: {result.species}")

    # Lookup expected and observed higher taxon
    for node in bold_tree.root.get_path(result.species):
        if node.taxonomic_rank == config.get('level'):
            result.exp_taxon = node
            break
    result.obs_taxon = run_localblast(record, ncbi_tree, bold_tree)

    # Compute marker quality metrics
    aligned_sequence = align_to_hmm(record)
    logging.debug(f"Sequence aligned to HMM: {aligned_sequence.seq}")
    amino_acid_sequence = translate_sequence(aligned_sequence)
    result.stop_codons = get_stop_codons(amino_acid_sequence)
    result.seq_length = marker_seqlength(aligned_sequence)
    result.ambiguities = num_ambiguous(aligned_sequence)

    # Done
    return result


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
