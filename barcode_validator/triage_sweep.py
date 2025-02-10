import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: Namespace with parsed arguments
    """
    parser = argparse.ArgumentParser(
        description='Select best DNA barcode sequences based on validation criteria'
    )
    parser.add_argument('-i', '--input-fasta', required=True, type=Path,
                        help='Input FASTA file with sequence variants')
    parser.add_argument('-v', '--validation', required=True, type=Path,
                        help='TSV file with validation results')
    parser.add_argument('-o', '--output-fasta', required=True, type=Path,
                        help='Output FASTA file with best sequences')
    parser.add_argument('-f', '--output-fails', required=True, type=Path,
                        help='Output TSV file with failed sequences')
    return parser.parse_args()


def setup_logging():
    """Configure logging to write to stderr."""
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def extract_process_id(sequence_id: str) -> str:
    """Extract the BOLD process ID from a sequence identifier.

    :param sequence_id: Full sequence identifier (e.g., 'BHNHM001-24_r_1_s_50')
    :return: Process ID portion (e.g., 'BHNHM001-24')
    """
    return sequence_id.split('_')[0]


def meets_basic_criteria(row: pd.Series) -> bool:
    """Check if a sequence meets the basic validation criteria.

    :param row: Pandas Series with validation data
    :return: True if sequence meets all basic criteria, False otherwise
    """
    return (
            row['ambig_basecount'] == 0 and
            row['error'] == 'None' and
            row['stop_codons'] == 0 and
            row['nuc_basecount'] >= 500 and
            row['nuc_full_basecount'] >= 500
    )


def select_best_sequences(validation_df: pd.DataFrame) -> pd.DataFrame:
    """Select the best sequence variant for each process ID.

    :param validation_df: DataFrame with validation results
    :return: DataFrame with only the best sequences
    """
    # Remove negative controls
    df = validation_df[~validation_df['sequence_id'].str.contains('-NC')]

    # Add process_id column
    df['process_id'] = df['sequence_id'].apply(extract_process_id)

    # Filter sequences that meet basic criteria
    valid_df = df[df.apply(meets_basic_criteria, axis=1)]

    if valid_df.empty:
        logging.warning("No sequences meet basic criteria")
        return pd.DataFrame()

    # Sort by criteria (first by process_id to group them, then by quality metrics)
    sorted_df = valid_df.sort_values(
        by=['process_id', 'ambig_full_basecount', 'nuc_full_basecount'],
        ascending=[True, True, False]  # True for ambig (lower better), False for nuc (higher better)
    )

    # Select first (best) sequence for each process_id
    return sorted_df.groupby('process_id').first().reset_index()


def write_fasta(sequences: Dict[str, str], best_df: pd.DataFrame, output_path: Path):
    """Write selected sequences to FASTA format.

    :param sequences: Dictionary mapping sequence IDs to sequences
    :param best_df: DataFrame with selected best sequences
    :param output_path: Path to output FASTA file
    """
    with open(output_path, 'w') as fh:
        for _, row in best_df.iterrows():
            seq_id = row['sequence_id']
            process_id = extract_process_id(seq_id)
            if seq_id in sequences:
                fh.write(f">{process_id}\n{sequences[seq_id]}\n")


def write_fails(validation_df: pd.DataFrame, best_df: pd.DataFrame, output_path: Path):
    """Write failed sequences to TSV format.

    :param validation_df: Original validation DataFrame
    :param best_df: DataFrame with selected best sequences
    :param output_path: Path to output TSV file
    """
    # Get process IDs that made it to best sequences
    successful_ids = set(best_df['process_id'])

    # Add process_id column if not already present
    if 'process_id' not in validation_df.columns:
        validation_df['process_id'] = validation_df['sequence_id'].apply(extract_process_id)

    # Filter for failed sequences
    failed_df = validation_df[~validation_df['process_id'].isin(successful_ids)]
    failed_df.to_csv(output_path, sep='\t', index=False)


def main():
    """Main function to run the sequence selection pipeline."""
    args = parse_args()
    setup_logging()

    # Read validation data
    logging.info(f"Reading validation data from {args.validation}")
    validation_df = pd.read_csv(args.validation, sep='\t')

    # Read sequences into dictionary
    logging.info(f"Reading sequences from {args.input_fasta}")
    sequences = {record.id: str(record.seq)
                 for record in SeqIO.parse(args.input_fasta, 'fasta')}

    # Select best sequences
    logging.info("Selecting best sequences")
    best_df = select_best_sequences(validation_df)

    # Write outputs
    logging.info(f"Writing best sequences to {args.output_fasta}")
    write_fasta(sequences, best_df, args.output_fasta)

    logging.info(f"Writing failed sequences to {args.output_fails}")
    write_fails(validation_df, best_df, args.output_fails)

    logging.info("Done")


if __name__ == '__main__':
    main()