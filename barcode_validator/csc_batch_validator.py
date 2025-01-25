import argparse
import csv
import re

from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from nbitk.Taxon import Taxon
from barcode_validator.core import BarcodeValidator
from barcode_validator.result import DNAAnalysisResult, DNAAnalysisResultSet
from barcode_validator.translation_tables import TaxonomyResolver, get_translation_table, Marker

# Always tell NCBI who you are
Entrez.email = "bioinformatics@naturalis.nl"

def main(table_file_path, logger, global_config):
    # Initialize BarcodeValidator
    validator = BarcodeValidator(global_config)
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

            # Clone the global config for each record because it may be modified
            config = global_config.local_clone()

            # Raw sequence and locally unique ID, which will go into the SeqRecord for v.validate_record()
            # and in the result object so that the output can be joined with the input
            seq_str = r['nuc']
            seq_id = r['local_id']
            seq_obj = SeqRecord(Seq(seq_str), id=seq_id)
            seq_obj.annotations['bcdm_fields'] = {'processid': seq_id}
            res = DNAAnalysisResult(seq_id, table_file_path)

            # Try to instantiate the Marker enum from the marker code. If it fails, log the error and skip.
            try:
                marker = Marker(r['marker_code'])

                # If the marker is anything other than COI-5P, this will throw exceptions later on. That's
                # fine, because the marker is not supported by the current implementation - but we do want
                # to log this so that the user knows that the marker is not supported.
                hmm_file = 'examples/' + marker.value + '.hmm'
                config.set('hmm_file', hmm_file)
            except ValueError as e:
                logger.error(f"Invalid marker code: {r['marker_code']} ({e})")
                res.error = f"Invalid marker code: {r['marker_code']} ({e})"
                continue

            # Attempt to resolve the taxonomic lineage at the specified ranks. If it fails, log the error and skip.
            specific_taxonomy = preprocess_taxa(config, r, res, tr, logger)
            if specific_taxonomy is None:
                logger.error(f"Taxon not found: {r['verbatim_identification']}")
                res.error = f"Taxon not found: {r['verbatim_identification']}"
                continue

            # Update the translation table using the specific taxonomy and marker.
            config.set('translation_table', get_translation_table(marker, specific_taxonomy))
            logger.info(f"Updated translation table: {config['translation_table']}")

            # Indicate in the result object what the expected taxon is at the specified level, which is, by
            # default, the family level but otherwise possibly the order level.
            if config['level'] in specific_taxonomy:
                res.exp_taxon = Taxon(specific_taxonomy[config['level']])
            else:
                logger.error(f"Expected taxon not found: {config['level']} in {specific_taxonomy}")
                res.error = f"Expected taxon not found: {config['level']} in {specific_taxonomy}"
                continue

            # Do the validation
            validator.validate_record(seq_obj, config, res)
            res.add_ancillary('is_valid', str(res.passes_all_checks()))
            results.append(res)

    # Create a DNAAnalysisResultSet object from the list of DNAAnalysisResult objects and print it
    rs = DNAAnalysisResultSet(results)
    print(rs)  # prints TSV
    logger.info("Analysis completed")


def preprocess_taxa(config, r, result, tr, logger):

    # Handle the rank variation: CSC has identifications at various ranks, including above the
    # family level. If so, this needs to be indicated in the result object so that the expected
    # and observed taxa are at the same level.
    taxon_rank = r['verbatim_rank']

    # Attempt to identify possible species binomials in records that have no rank information.
    pattern = r'^[A-Z][a-z]+ [a-z]+$'
    if taxon_rank == 'null' and re.match(pattern, r['verbatim_identification']):
        taxon_rank = 'species'
        logger.warning(f"Assuming {r['verbatim_identification']} is a species")

    # If the rank is in none of these, it is probably an order, so far as we've been able to see.
    # We can then only validate at that higher level instead of the default family level. For
    # this, both the config object needs to be updated and the result object needs to be updated.
    if taxon_rank not in ['Family', 'Genus', 'Species', 'species', 'null']:
        config.set('level', str(taxon_rank).lower())  # Probably 'order'
        logger.warning(f"Taxon {r['verbatim_identification']} can only be validated at the config['level'] level")
    else:
        config.set('level', 'family')
    result.level = config.get('level')  # Will be identification_rank in the result object

    # Try to resolve the taxonomic lineage at the specified ranks. This may return None, in which
    # case the record will be skipped. The webservice optionally takes a kingdom name to avoid
    # homonyms. If the kingdom is not specified, it will be set to None.
    if r['verbatim_kingdom'] == 'null':
        r['verbatim_kingdom'] = None
        logger.warning(f"{r['verbatim_identification']} has no kingdom, this may lead to homonyms")
    ranks = ['kingdom', 'phylum', 'class', 'order', 'family']
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
