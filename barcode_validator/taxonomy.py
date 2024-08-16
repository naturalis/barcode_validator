import logging
import tempfile
import subprocess
import os
import time
import pandas as pd
import tarfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
from config import Config
from nbt.Phylo.BOLDXLSXIO import Parser as BOLDParser
from nbt.Phylo.NCBITaxdmp import Parser as NCBIParser


def read_bold_taxonomy(spreadsheet):
    logging.info("Reading BOLD taxonomy")
    return BOLDParser(spreadsheet).parse()


def read_ncbi_taxonomy(tarfile):
    logging.info("Reading NCBI taxonomy")
    tar = tarfile.open(tarfile, "r:gz")
    return NCBIParser(tar).parse()


def run_seqid(sequence, ncbi_tree):
    logging.info("Running sequence ID check")
    config = Config()
    method = config.get('idcheck')
    if method == 'boldigger2':
        return run_boldigger2(sequence)
    elif method == 'blast':
        return run_blast(sequence)
    elif method == 'localblast':
        return run_localblast(sequence, ncbi_tree)
    else:
        logging.error(f"Invalid ID check method '{method}' specified in config file")
    pass


def run_localblast(sequence, tree):
    logging.info("Running local BLASTN...")

    # Create a temporary file for the input sequence
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as temp_input:
        SeqIO.write(sequence, temp_input, "fasta")
        temp_input_name = temp_input.name
    blast_result = f"{temp_input_name}.tsv"  # output file name, as blast TSV (see -outfmt option)

    # Prepare local BLASTN run
    config = Config()
    os.environ['BLASTDB_LMDB_MAP_SIZE'] = config.get('BLASTDB_LMDB_MAP_SIZE')
    try:

        # Run BLASTN
        subprocess.run(['blastn',
                        '-db', config.get('blast_db'),
                        '-num_threads', config.get('num_threads'),
                        '-evalue', config.get('evalue'),
                        '-max_target_seqs', config.get('max_target_seqs'),
                        '-word_size', config.get('word_size'),
                        '-query', temp_input_name,
                        '-task', 'megablast',
                        '-outfmt', "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore staxids",
                        '-out', blast_result
                        ],
                       check=True)

        # Parse BLAST result
        distinct_taxids = set()
        with open(blast_result, 'r') as file:
            for line in file:
                columns = line.strip().split('\t')
                if columns:
                    taxid_field = columns[-1]
                    taxids = taxid_field.split(';')
                    distinct_taxids.update(taxid.strip() for taxid in taxids if taxid.strip())
        return collect_lineages(distinct_taxids, tree)

    # Handle exception
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running local BLASTN: {e}")
        raise


def collect_lineages(taxids, tree):
    lineages = []
    tips = []
    for tip in tree.get_terminals():
        taxid = tip.guids['taxon']
        if taxid in taxids:
            tips.append(tip)
    for tip in tips:
        lineage = []
        for node in tree.root.get_path(tip):
            if node.rank == str(config.get('level')).lower():
                lineage.append(node.name)
        lineages.append(lineage)
    return lineages


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
