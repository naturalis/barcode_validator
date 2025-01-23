import argparse
import csv

from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from barcode_validator.core import BarcodeValidator
from barcode_validator.result import DNAAnalysisResult, DNAAnalysisResultSet
from barcode_validator.translation_tables import TaxonomyResolver, get_translation_table, Marker

# Always tell NCBI who you are
Entrez.email = "bioinformatics@naturalis.nl"

def main(table_file_path, logger, config):
    # Initialize BarcodeValidator
    validator = BarcodeValidator(config)
    validator.initialize()

    # Initialize TaxonResolver
    tr = TaxonomyResolver(Entrez.email)

    # Print header
    logger.info(f"Starting analysis for file: {table_file_path}")

    # Read the table file as TSV, iterate over the records as dictionaries
    results = []
    with open(table_file_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for r in reader:

            # Raw sequence and locally unique ID, which will go into the SeqRecord for v.validate_record()
            # and in the result object so that the output can be joined with the input
            seq_str = r['nuc']
            seq_id = r['local_id']
            seq_obj = SeqRecord(Seq(seq_str), id=seq_id)
            res = DNAAnalysisResult(seq_id)

            # Try to instantiate the Marker enum from the marker code. If it fails, log the error and skip.
            try:
                marker = Marker(r['marker_code'])
            except ValueError as e:
                logger.error(f"Invalid marker code: {r['marker_code']} ({e})")
                res.error = f"Invalid marker code: {r['marker_code']} ({e})"
                continue

            # Attempt to resolve the taxonomic lineage at the specified ranks. If it fails, log the error and skip.
            specific_taxonomy = preprocess_taxa(config, r, res, tr)
            if specific_taxonomy is None:
                logger.error(f"Taxon not found: {r['verbatim_identification']}")
                res.error = f"Taxon not found: {r['verbatim_identification']}"
                continue

            # Update the translation table using the specific taxonomy and marker.
            config['translation_table'].update(get_translation_table(marker, specific_taxonomy))

            # Indicate in the result object what the expected taxon is at the specified level, which is, by
            # default, the family level but otherwise possibly the order level.
            if config['level'] in specific_taxonomy:
                res.exp_taxon = specific_taxonomy[config['level']]
            else:
                logger.error(f"Expected taxon not found: {config['level']} in {specific_taxonomy}")
                res.error = f"Expected taxon not found: {config['level']} in {specific_taxonomy}"
                continue

            # Do the validation
            validator.validate_record(seq_obj, config, res)
            res.add_ancillary('is_valid', str(res.passes_all_checks()))
            results.append(res)

    # Validate the FASTA file
    rs = DNAAnalysisResultSet(results)

    print(rs)  # print TSV results

    logger.info("Analysis completed")


def preprocess_taxa(config, r, result, tr):

    # Handle the rank variation: CSC has identifications at various ranks, including above the
    # family level. If so, this needs to be indicated in the result object so that the expected
    # and observed taxa are at the same level.
    taxon_rank = r['verbatim_rank']
    if taxon_rank not in ['family', 'genus', 'species']:
        config['level'] = str(taxon_rank).lower()  # Probably 'order'
    else:
        config['level'] = 'family'
    result.level = config['level']  # Will be identification_rank in the result object

    # Try to resolve the taxonomic lineage at the specified ranks. If it fails, log the error and skip.
    ranks = ['phylum', 'class', 'order', 'family']
    specific_taxonomy = tr.get_lineage_at_ranks(r['verbatim_identification'], ranks, r['verbatim_kingdom'])

    return specific_taxonomy


if __name__ == "__main__":
    # Process command line arguments
    parser = argparse.ArgumentParser(description="Analyze DNA sequences from a CSC dump file.")
    parser.add_argument("-t", "--table_file", required=True, help="Path to the input dump file")
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
    main(args.table_file, main_logger, main_config)
