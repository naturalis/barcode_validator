import argparse
import logging
import sys
from pathlib import Path
from typing import List, Set, Dict

import pandas as pd

# Required columns that must be present in each TSV file
REQUIRED_COLUMNS = {
    'sequence_id',
    'error',
    'ambig_full_basecount',
    'ambig_basecount',
    'stop_codons',
    'nuc_basecount',
    'identification',
    'obs_taxon'
}


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: Namespace with parsed arguments
    """
    parser = argparse.ArgumentParser(
        description='Merge MGE TSV files with appropriate sequence_id modifications'
    )
    parser.add_argument('-i', '--input-dir', required=True, type=Path,
                        help='Input directory containing MGE TSV files')
    parser.add_argument('-o', '--output-file', required=True, type=Path,
                        help='Output TSV file path')
    parser.add_argument('--dry-run', action='store_true',
                        help='Perform validation without writing output')
    return parser.parse_args()


def setup_logging():
    """Configure logging to write to stderr."""
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )


def get_mge_files(input_dir: Path) -> List[Path]:
    """Get all MGE TSV files from the input directory.

    :param input_dir: Path to directory containing TSV files
    :return: List of paths to MGE TSV files
    """
    return sorted(
        path for path in input_dir.glob('*.tsv')
        if path.name.startswith(('mge_fastp_', 'mge_standard_'))
    )


def validate_columns(df: pd.DataFrame, file_path: Path) -> bool:
    """Validate that all required columns are present in the DataFrame.

    :param df: DataFrame to validate
    :param file_path: Path to source file (for logging)
    :return: True if all required columns are present, False otherwise
    """
    missing_columns = REQUIRED_COLUMNS - set(df.columns)
    if missing_columns:
        logging.error(f"File {file_path} is missing required columns: {missing_columns}")
        return False
    return True


def modify_sequence_id(df: pd.DataFrame, file_path: Path) -> pd.DataFrame:
    """Modify sequence_id values based on file prefix.

    :param df: DataFrame containing sequence_id column
    :param file_path: Path to source file
    :return: DataFrame with modified sequence_id values
    """
    # Determine suffix based on filename
    if file_path.name.startswith('mge_fastp_'):
        suffix = '_fastp'
    else:
        suffix = '_standard'

    # Create a copy to avoid modifying the original
    df = df.copy()

    # Strip whitespace and add suffix
    df['sequence_id'] = df['sequence_id'].astype(str).str.strip() + suffix

    return df


def process_files(files: List[Path], dry_run: bool = False) -> pd.DataFrame:
    """Process all MGE TSV files and combine them into a single DataFrame.

    :param files: List of paths to TSV files
    :param dry_run: If True, only validate files without combining
    :return: Combined DataFrame or empty DataFrame if dry_run is True
    """
    # Initialize set to track all columns
    all_columns: Set[str] = set()

    # First pass: validate files and collect all column names
    valid_files: List[Path] = []
    for file_path in files:
        logging.info(f"Validating {file_path}")
        try:
            df = pd.read_csv(file_path, sep='\t')
            if validate_columns(df, file_path):
                valid_files.append(file_path)
                all_columns.update(df.columns)
            else:
                logging.error(f"Skipping {file_path} due to missing columns")
        except Exception as e:
            logging.error(f"Error reading {file_path}: {str(e)}")

    if dry_run:
        logging.info("Dry run completed")
        return pd.DataFrame()

    # Second pass: read and combine valid files
    dfs: List[pd.DataFrame] = []
    for file_path in valid_files:
        logging.info(f"Processing {file_path}")
        try:
            # Read file
            df = pd.read_csv(file_path, sep='\t')

            # Modify sequence_ids
            df = modify_sequence_id(df, file_path)

            # Add missing columns with 'None' values
            missing_cols = all_columns - set(df.columns)
            for col in missing_cols:
                df[col] = 'None'

            dfs.append(df)

        except Exception as e:
            logging.error(f"Error processing {file_path}: {str(e)}")

    if not dfs:
        logging.error("No valid files to process")
        return pd.DataFrame()

    # Combine all DataFrames
    return pd.concat(dfs, ignore_index=True)


def main():
    """Main function to run the MGE file merger pipeline."""
    args = parse_args()
    setup_logging()

    # Validate input directory
    if not args.input_dir.is_dir():
        logging.error(f"Input directory does not exist: {args.input_dir}")
        sys.exit(1)

    # Get list of MGE files
    files = get_mge_files(args.input_dir)
    if not files:
        logging.error(f"No MGE TSV files found in {args.input_dir}")
        sys.exit(1)

    logging.info(f"Found {len(files)} MGE TSV files")

    # Process files
    result_df = process_files(files, args.dry_run)

    if args.dry_run:
        logging.info("Dry run completed, no output file written")
        return

    if result_df.empty:
        logging.error("No data to write")
        sys.exit(1)

    # Write output
    logging.info(f"Writing combined data to {args.output_file}")
    result_df.to_csv(args.output_file, sep='\t', index=False)
    logging.info("Done")


if __name__ == '__main__':
    main()