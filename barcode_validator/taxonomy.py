import logging
import tempfile
import subprocess
import os
import time
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
from config import Config


def fetch_bold_taxonomy(spreadsheet):
    # Check if the file exists
    if not os.path.exists(spreadsheet):
        raise FileNotFoundError(f"The file {spreadsheet} does not exist.")

    # Read the Excel file
    xl = pd.ExcelFile(spreadsheet)

    # Check if required sheets exist
    required_sheets = ['Lab Sheet', 'Taxonomy']
    if not all(sheet in xl.sheet_names for sheet in required_sheets):
        raise ValueError("The Excel file must contain 'Lab Sheet' and 'Taxonomy' tabs.")

    # Read 'Lab Sheet' tab
    lab_sheet = pd.read_excel(xl, sheet_name='Lab Sheet', header=2)

    # Create mapping dict for Lab Sheet
    lab_mapping = dict(zip(lab_sheet['Sample ID'], lab_sheet['Process ID']))

    # Read 'Taxonomy' tab
    taxonomy = pd.read_excel(xl, sheet_name='Taxonomy', header=2)

    # List of taxonomy columns. These are ucfirst in the BOLD Excel file.
    taxonomy_columns = ['Phylum', 'Class', 'Order', 'Family', 'Subfamily', 'Tribe', 'Genus', 'Species', 'Subspecies']

    # Create mapping dict for Taxonomy
    taxonomy_mapping = {}
    for _, row in taxonomy.iterrows():
        sample_id = row['Sample ID']
        taxonomy_dict = {col.lower(): row[col] for col in taxonomy_columns}
        taxonomy_mapping[sample_id] = taxonomy_dict

    # Combine the two mappings
    result = {}
    for sample_id, process_id in lab_mapping.items():
        if sample_id in taxonomy_mapping:
            result[process_id] = taxonomy_mapping[sample_id]

    return result


def run_seqid(sequence):
    logging.info("Running sequence ID check")
    config = Config()
    method = config.get('idcheck')
    if method == 'boldigger2':
        return run_boldigger2(sequence)
    elif method == 'blast':
        return run_blast(sequence)
    else:
        logging.error(f"Invalid ID check method '{method}' specified in config file")
    pass


def run_blast(sequence):

    # Get config singleton
    config = Config()
    Entrez.email = config.get('email')
    target_level = config.get('level')

    # Run BLASTN, pparse result
    logging.info("Running BLASTN...")
    blast_result = NCBIWWW.qblast("blastn", "nt", sequence.seq)
    blast_records = NCBIXML.parse(blast_result)
    record = next(blast_records)
    top_hits = record.alignments[:10]  # Get top 10 hits

    # Fetch taxonomic lineages for top hits
    lineages = []
    for hit in top_hits:
        accession = hit.accession
        try:
            logging.info(f"Going to fetch nucleotide record for {accession}")
            summary = Entrez.esummary(db="nucleotide", id=accession)
            summary_record = Entrez.read(summary)[0]
            taxon_id = summary_record['TaxId']

            logging.info(f"Going to fetch taxonomy for {taxon_id}")
            handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
            taxonomy_record = Entrez.read(handle)[0]
            lineage = {rank['Rank']: rank['ScientificName'] for rank in taxonomy_record['LineageEx']}
            logging.debug(lineage)

            lineages.append(lineage)
            time.sleep(0.1)  # To avoid overloading NCBI servers
        except Exception as e:
            logging.error(f"Error fetching taxonomy for {accession}: {e}")

    # Apply majority rule consensus at the specified level
    taxa_at_level = list(set(lineage.get(target_level) for lineage in lineages if target_level in lineage))
    if taxa_at_level:
        logging.info(f'Taxa at {target_level} level: {taxa_at_level}')
        return taxa_at_level
    else:
        logging.warning(f"No taxa found at {target_level} level")
        return ["Unknown"]


def run_boldigger2(sequence):
    logging.info("Running boldigger2")

    # Create a temporary file for the input sequence
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as temp_input:
        SeqIO.write(sequence, temp_input, "fasta")
        temp_input_name = temp_input.name

    # Run boldigger2
    # TODO: fetch username and password from environment variables
    try:
        subprocess.run(['boldigger2', 'identify', temp_input_name], check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running boldigger2: {e}")
        raise

    # Process the resulting Excel file, which is named {temp_input_name}_identification_result.xlsx
    config = Config()
    basename = os.path.splitext(temp_input_name)[0]
    xlsx_filename = f"{basename}_identification_result.xlsx"
    df = pd.read_excel(xlsx_filename)
    query_rows = df.loc[df['ID'] == 'query']
    column_name = str(config.get('level')).lower()
    for _, row in query_rows.iterrows():
        logging.info(f'Found taxonomy: {row[column_name]}')
        return [row[column_name]]
