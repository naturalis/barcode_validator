import argparse
import csv
import logging
import re

from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from nbitk.config import Config
from nbitk.logger import get_formatted_logger
from barcode_validator.core import BarcodeValidator
from barcode_validator.result import DNAAnalysisResult, DNAAnalysisResultSet
from barcode_validator.translation_tables import TaxonomyResolver, Marker

# Always tell NCBI who you are
Entrez.email = "bioinformatics@naturalis.nl"

def main(table_file_path, logger, global_config):
    # Initialize BarcodeValidator and TaxonomyResolver objects
    validator = BarcodeValidator(global_config)
    validator.initialize()
    tr = TaxonomyResolver(Entrez.email, logger, validator.ncbi_tree)

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
            seq_obj = SeqRecord(Seq(seq_str), id=seq_id, description='')
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
            taxon = preprocess_taxa(config, r, res, logger, tr)
            if taxon is None:
                logger.error(f"Taxon not found: {r['verbatim_identification']}")
                res.error = f"Taxon not found: {r['verbatim_identification']}"
                continue
            res.species = taxon

            # Update the translation table using the specific taxonomy and marker.
            config.set('translation_table', tr.get_translation_table(marker, taxon))
            logger.info(f"Updated translation table: {config['translation_table']}")

            # Indicate in the result object what the expected taxon is at the specified level, which is, by
            # default, the family level but otherwise possibly the order level.
            for node in validator.ncbi_tree.root.get_path(taxon):
                if str(node.taxonomic_rank).lower() == str(config['level']).lower():
                    res.exp_taxon = node
                    break
            if res.exp_taxon is None:
                logger.error(f"Expected taxon not found: {config['level']} for {taxon}")
                res.error = f"Expected taxon not found: {config['level']} in {taxon}"
                continue

            # Fetch the taxon ID for the ancestor of the expected taxon at the specified level
            for node in validator.ncbi_tree.root.get_path(taxon):
                if str(node.taxonomic_rank).lower() == str(config['constrain']).lower():
                    config.set('constraint_taxid', node.guids['taxon'])
                    break
            if config.get('constraint_taxid') is None:
                logger.warning(f"Constraint taxon not found: {config['constrain']} for {taxon}")
                logger.warning(f"Using Eukaryota (taxon:2759) instead. This will be a lot slower.")
                config.set('constraint_taxid', 2759)

            # Do the validation
            validator.validate_record(seq_obj, config, res)
            res.add_ancillary('is_valid', str(res.passes_all_checks()))
            results.append(res)

    # Create a DNAAnalysisResultSet object from the list of DNAAnalysisResult objects and print it
    rs = DNAAnalysisResultSet(results)
    print(rs)  # prints TSV
    logger.info("Analysis completed")


def preprocess_taxa(config: Config, record: dict, result: DNAAnalysisResult, logger: logging.Logger, tr: TaxonomyResolver):

    # Handle the rank variation: CSC has identifications at various ranks, including above the
    # family level. If so, this needs to be indicated in the result object so that the expected
    # and observed taxa are at the same level.
    taxon_rank = record['verbatim_rank']

    # Attempt to identify possible species binomials in records that have no rank information.
    pattern = r'^[A-Z][a-z]+ [a-z]+$'
    if taxon_rank == 'null' and re.match(pattern, record['verbatim_identification']):
        taxon_rank = 'species'
        logger.warning(f"Assuming {record['verbatim_identification']} is a species")

    # If the rank is in none of these, it is probably an order, so far as we've been able to see.
    # We can then only validate at that higher level instead of the default family level. For
    # this, both the config object needs to be updated and the result object needs to be updated.
    if taxon_rank not in ['Family', 'Genus', 'Species', 'species', 'null']:
        config.set('level', str(taxon_rank).lower())  # Probably 'order'
        logger.warning(f"Taxon {record['verbatim_identification']} can only be validated at the config['level'] level")
    else:
        config.set('level', 'family')
    result.level = config.get('level')  # Will be identification_rank in the result object

    # Try to resolve the verbatim_identification wrt the NCBI taxonomy via the Entrez service. This
    # may return None, in which case the record will be skipped.
    taxon = tr.get_taxon(record['verbatim_identification'])
    return taxon


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
        main_config.set('log_level', args.verbosity or 'INFO')
        main_logger = get_formatted_logger(__name__, main_config)
    except ValueError as e:
        print(f"Error setting up logging: {e}")
        exit(1)

    # Run the main analysis
    main(args.table_file, main_logger, main_config)
