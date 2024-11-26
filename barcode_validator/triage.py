import argparse
import logging
import os
import csv
import sys
from Bio import SeqIO
from pathlib import Path


def setup_logging(verbosity):
    """
    Configure logging format and level based on verbosity flag count.

    Args:
        verbosity (int): Number of verbose flags (-v) provided:
            0: WARNING level (default)
            1: INFO level
            2: DEBUG level
    """
    levels = {
        0: logging.WARNING,
        1: logging.INFO,
        2: logging.DEBUG
    }
    # If verbosity is higher than 2, cap it at DEBUG level
    level = levels.get(min(verbosity, 2), logging.DEBUG)

    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logging.debug(f"Logging level set to: {logging.getLevelName(level)}")


def parse_args():
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace: Object containing the following attributes:
            - tsv (str): Path to input TSV file with sequence metadata
            - fasta (str): Path to input FASTA file with sequences
            - good_dir (str): Path to output directory for good sequences
            - bad_dir (str): Path to output directory for bad sequences
            - verbose (int): Count of verbose flags (0-2)
    """
    parser = argparse.ArgumentParser(description='Triage FASTA sequences based on quality criteria.')
    parser.add_argument('--tsv', required=True, help='Input TSV file with sequence metadata')
    parser.add_argument('--fasta', required=True, help='Input FASTA file with sequences')
    parser.add_argument('--good_dir', required=True, help='Output directory for good sequences')
    parser.add_argument('--bad_dir', required=True, help='Output directory for bad sequences')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='Increase output verbosity (use -v for INFO, -v -v for DEBUG)')
    return parser.parse_args()


def read_tsv_data(tsv_file):
    """
    Read TSV file and return data as dictionary keyed by process_id.

    Args:
        tsv_file (str): Path to the input TSV file containing sequence metadata.

    Returns:
        dict: Dictionary where keys are process_ids and values are dictionaries
             containing the corresponding row data from the TSV file.
    """
    data = {}
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        logging.debug(f"Reading TSV file: {tsv_file}")
        for row in reader:
            try:
                data[row['sequence_id']] = row
            except KeyError:
                print(f"Error: Missing process_id in file {tsv_file}")
                sys.exit(1)
    logging.debug(f"Read {len(data)} records from TSV file")
    return data


def check_sequence(seq_id, metadata):
    """
    Check if sequence meets quality criteria.

    Args:
        seq_id (str): Identifier of the sequence being checked.
        metadata (dict): Dictionary containing sequence metadata with the following keys:
            - nuc_basecount (str): Number of nucleotide bases
            - stop_codons (str): Number of stop codons
            - ambig_basecount (str): Number of ambiguous bases
            - obs_taxon (str): Comma-separated list of observed taxa
            - identification (str): Expected taxon identification

    Returns:
        tuple: A pair containing:
            - bool: True if sequence passes all criteria, False otherwise
            - list: List of strings indicating which criteria failed (empty if all passed)
    """
    failed_criteria = []
    logging.debug(f"Checking sequence {seq_id}")

    # Check nuc_basecount >= 500
    nuc_count = 0
    if metadata['nuc_basecount'] != 'None':
        nuc_count = int(metadata['nuc_basecount'])
    if nuc_count >= 500:
        logging.info(f"Sequence {seq_id} passes nucleotide base count check ({nuc_count} bases)")
    else:
        logging.warning(f"Sequence {seq_id} fails nucleotide base count check ({nuc_count} bases)")
        failed_criteria.append('nuc_basecount')

    # Check stop_codons == 0
    stop_count = int(metadata['stop_codons'])
    if stop_count == 0:
        logging.info(f"Sequence {seq_id} passes stop codons check (no stop codons)")
    else:
        logging.warning(f"Sequence {seq_id} fails stop codons check ({stop_count} stop codons)")
        failed_criteria.append('stop_codons')

    # Check ambig_basecount <= 6
    ambig_count = 7
    if metadata['ambig_basecount'] != 'None':
        ambig_count = int(metadata['ambig_basecount'])
    if ambig_count <= 6:
        logging.info(f"Sequence {seq_id} passes ambiguous base count check ({ambig_count} ambiguous bases)")
    else:
        logging.warning(f"Sequence {seq_id} fails ambiguous base count check ({ambig_count} ambiguous bases)")
        failed_criteria.append('ambig_basecount')

    # Check if identification is in obs_taxon
    obs_taxa = {taxon.strip() for taxon in metadata['obs_taxon'].split(',')}
    if metadata['identification'] in obs_taxa:
        logging.info(f"Sequence {seq_id} passes taxonomy check (identification: {metadata['identification']})")
    else:
        logging.warning(
            f"Sequence {seq_id} fails taxonomy check (identification '{metadata['identification']}' "
            f"not found in observed taxa: {', '.join(obs_taxa)})"
        )
        failed_criteria.append('taxonomy')

    logging.debug(
        f"Sequence {seq_id} check complete. Failed criteria: {', '.join(failed_criteria) if failed_criteria else 'none'}")
    return len(failed_criteria) == 0, failed_criteria


def process_sequences(args):
    """
    Process sequences and write to appropriate output directories.

    Args:
        args (argparse.Namespace): Command line arguments containing:
            - tsv (str): Path to input TSV file
            - fasta (str): Path to input FASTA file
            - good_dir (str): Path to output directory for good sequences
            - bad_dir (str): Path to output directory for bad sequences

    Side effects:
        - Creates output directories if they don't exist
        - Writes filtered FASTA and TSV files to good_dir and bad_dir
        - Logs processing results using the logging module
    """
    # Create output directories if they don't exist
    Path(args.good_dir).mkdir(parents=True, exist_ok=True)
    Path(args.bad_dir).mkdir(parents=True, exist_ok=True)
    logging.debug(f"Created output directories: {args.good_dir}, {args.bad_dir}")

    # Read TSV data
    metadata = read_tsv_data(args.tsv)

    # Prepare output files
    good_fasta = os.path.join(args.good_dir, os.path.basename(args.fasta))
    bad_fasta = os.path.join(args.bad_dir, os.path.basename(args.fasta))
    good_tsv = os.path.join(args.good_dir, os.path.basename(args.tsv))
    bad_tsv = os.path.join(args.bad_dir, os.path.basename(args.tsv))

    logging.debug(f"Output files prepared:\n"
                  f"Good FASTA: {good_fasta}\n"
                  f"Bad FASTA: {bad_fasta}\n"
                  f"Good TSV: {good_tsv}\n"
                  f"Bad TSV: {bad_tsv}")

    # Process sequences
    good_records = []
    bad_records = []
    good_metadata = []
    bad_metadata = []

    logging.info(f"Starting sequence processing from {args.fasta}")
    for record in SeqIO.parse(args.fasta, 'fasta'):
        seq_id = record.id
        logging.debug(f"Processing sequence: {seq_id}")

        if seq_id not in metadata:
            logging.error(f"No metadata found for sequence {seq_id}")
            continue
        try:
            passes_check, failed_criteria = check_sequence(seq_id, metadata[seq_id])
        except KeyError as e:
            logging.error(f"Error processing sequence {seq_id}: {e} in {args.tsv}")
            continue

        if passes_check:
            good_records.append(record)
            good_metadata.append(metadata[seq_id])
        else:
            bad_records.append(record)
            bad_metadata.append(metadata[seq_id])

    # Write FASTA files
    SeqIO.write(good_records, good_fasta, 'fasta')
    SeqIO.write(bad_records, bad_fasta, 'fasta')
    logging.debug(f"Written FASTA files:\n"
                  f"Good sequences: {good_fasta} ({len(good_records)} records)\n"
                  f"Bad sequences: {bad_fasta} ({len(bad_records)} records)")

    # Write TSV files
    fieldnames = next(csv.DictReader(open(args.tsv), delimiter='\t')).keys()

    def write_tsv(filename, data):
        with open(filename, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(data)
        logging.debug(f"Written TSV file: {filename} ({len(data)} records)")

    write_tsv(good_tsv, good_metadata)
    write_tsv(bad_tsv, bad_metadata)

    logging.info(
        f"Processing complete. Processed {len(good_records)} good sequences and {len(bad_records)} bad sequences")


def main():
    """
    Main function to run the script.

    This function:
    1. Parses command line arguments
    2. Sets up logging configuration based on verbosity level
    3. Processes and triages sequences based on quality criteria
    """
    args = parse_args()
    setup_logging(args.verbose)
    process_sequences(args)


if __name__ == "__main__":
    main()