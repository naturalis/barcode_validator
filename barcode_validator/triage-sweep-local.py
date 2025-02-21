import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
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
    """Extract the BOLD process ID from a sequence identifier."""
    return sequence_id.split('_')[0]


def check_validation_criteria(row: pd.Series) -> List[str]:
    """Check each validation criterion separately and return list of failures.
    
    :param row: Pandas Series with validation data
    :return: List of failed criteria descriptions
    """
    failures = []
    
    # Convert "None" string or NaN to None for proper comparison
    error_value = row['error']
    if pd.isna(error_value):
        error_value = "None"
    
    if error_value.lower() != "none":  # Case-insensitive comparison
        failures.append(f"Has error: {error_value}")
    
    if row['ambig_full_basecount'] > 6:
        failures.append(f"Too many full ambiguous bases: {row['ambig_full_basecount']} > 6")
        
    if row['ambig_basecount'] > 6:
        failures.append(f"Too many ambiguous bases: {row['ambig_basecount']} > 6")
        
    if row['stop_codons'] > 0:
        failures.append(f"Contains stop codons: {row['stop_codons']}")
        
    if row['nuc_basecount'] < 500:
        failures.append(f"Sequence too short: {row['nuc_basecount']} < 500")
        
    if row['identification'] not in str(row['obs_taxon']).split(','):
        failures.append(f"Identification mismatch: {row['identification']} not in {row['obs_taxon']}")
    
    return failures


def meets_basic_criteria(row: pd.Series) -> bool:
    """Check if a sequence meets the basic validation criteria."""
    return len(check_validation_criteria(row)) == 0


def select_best_sequences(validation_df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, List[str]]]:
    """Select the best sequence variant for each process ID and track failures.
    
    :param validation_df: DataFrame with validation results
    :return: Tuple of (DataFrame with best sequences, Dict mapping process IDs to failure reasons)
    """
    # Remove negative controls and create a copy
    df = validation_df[~validation_df['sequence_id'].str.contains('-NC')].copy()

    # Convert numeric columns
    numeric_cols = ['ambig_basecount', 'ambig_full_basecount', 'nuc_basecount', 'nuc_full_basecount', 'stop_codons']
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Add process_id column
    df.loc[:, 'process_id'] = df['sequence_id'].apply(extract_process_id)
    
    # Track failures for each sequence
    failure_reasons = {}
    for _, row in df.iterrows():
        failures = check_validation_criteria(row)
        if failures:
            if row['process_id'] not in failure_reasons:
                failure_reasons[row['process_id']] = []
            failure_reasons[row['process_id']].extend(failures)

    # Filter sequences that meet basic criteria
    valid_df = df[df.apply(meets_basic_criteria, axis=1)]

    if valid_df.empty:
        logging.warning("No sequences meet basic criteria")
        # Log specific failures
        for process_id, reasons in failure_reasons.items():
            unique_reasons = list(set(reasons))  # Remove duplicates
            logging.warning(f"Process ID {process_id} failed due to:")
            for reason in unique_reasons:
                logging.warning(f"  - {reason}")
        return pd.DataFrame(), failure_reasons

    # Sort by criteria
    sorted_df = valid_df.sort_values(
        by=['process_id', 'ambig_full_basecount', 'nuc_full_basecount'],
        ascending=[True, True, False]
    )

    # Select first (best) sequence for each process_id
    return sorted_df.groupby('process_id').first().reset_index(), failure_reasons


def write_fasta(sequences: Dict[str, str], best_df: pd.DataFrame, output_path: Path):
    """Write selected sequences to FASTA format."""
    with open(output_path, 'w') as fh:
        for _, row in best_df.iterrows():
            seq_id = row['sequence_id']
            process_id = extract_process_id(seq_id)
            if seq_id in sequences:
                fh.write(f">{process_id}\n{sequences[seq_id]}\n")


def write_fails(validation_df: pd.DataFrame, best_df: pd.DataFrame, failure_reasons: Dict[str, List[str]], output_path: Path):
    """Write failed sequences to TSV format with failure reasons."""
    # Create a copy and add process_id column
    failed_df = validation_df.copy()
    failed_df.loc[:, 'process_id'] = failed_df['sequence_id'].apply(extract_process_id)

    if not best_df.empty:
        successful_ids = set(best_df['process_id'])
        failed_df = failed_df[~failed_df['process_id'].isin(successful_ids)]
    
    # Add failure reasons column
    failed_df['failure_reasons'] = failed_df['process_id'].map(
        lambda x: '; '.join(set(failure_reasons.get(x, [])))
    )

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
    best_df, failure_reasons = select_best_sequences(validation_df)

    # Write outputs
    logging.info(f"Writing best sequences to {args.output_fasta}")
    write_fasta(sequences, best_df, args.output_fasta)

    logging.info(f"Writing failed sequences to {args.output_fails}")
    write_fails(validation_df, best_df, failure_reasons, args.output_fails)

    # Summary statistics
    total_sequences = len(validation_df)
    failed_sequences = len(validation_df) - len(best_df)
    logging.info(f"\nSummary:")
    logging.info(f"Total sequences processed: {total_sequences}")
    logging.info(f"Failed sequences: {failed_sequences}")
    logging.info(f"Successful sequences: {len(best_df)}")

    logging.info("Done")


if __name__ == '__main__':
    main()
