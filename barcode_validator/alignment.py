import logging
import tempfile
import subprocess
from copy import deepcopy
from config import Config
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable


"""
alignment.py

This module contains the functions for technical QA on sequence results, including operations
to align/unalign the sequence against an HMM, calculate the length of the sequence inside the
barcode marker region, and calculate the number of stop codons, which is intended as an 
indication of off-target amplification of a pseudogene.
"""


def align_to_hmm(sequence):
    logging.info("Aligning sequence to HMM")
    config = Config()
    hmm_file = config.get('hmm_file')

    # Create a temporary file for the input sequence
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as temp_input:
        SeqIO.write(sequence, temp_input, "fasta")
        temp_input.close()
        infile = temp_input.name

        # Run hmmalign with the input file
        try:

            # Create a temporary file for the output alignment
            with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.sto') as temp_output:
                outfile = temp_output.name
                subprocess.run(['hmmalign', '--trim', '-o', outfile, '--outformat', 'Stockholm', hmm_file, infile])
                logging.debug(f"Going to parse outfile {outfile}")
                aligned_sequence = next(SeqIO.parse(outfile, 'stockholm'))

        except subprocess.CalledProcessError as e:
            logging.error(f"Error running hmmalign: {e}")
            raise

    logging.debug(f"Aligned sequence length: {len(aligned_sequence)}")

    return aligned_sequence


def marker_seqlength(sequence):
    logging.info("Calculating non-missing bases within marker region")

    # Clone the sequence, remove gaps and missing in the sequence
    cloned_seq = deepcopy(sequence)
    cloned_seq.letter_annotations = {}
    cloned_seq.seq = cloned_seq.seq.replace('-', '')
    cloned_seq.seq = cloned_seq.seq.replace('N', '')
    cloned_seq.seq = cloned_seq.seq.replace('n', '')
    logging.debug(f"Within marker sequence length: {len(cloned_seq)}")
    return len(cloned_seq)


def unalign_sequence(sequence):
    logging.info("Removing gaps from aligned sequence")
    unaligned_sequence = sequence.seq.replace('-', '')
    unaligned_sequence = unaligned_sequence.replace('~', '')
    sequence.seq = unaligned_sequence
    logging.debug(f"Unaligned sequence length: {len(unaligned_sequence)}")
    return sequence


def translate_sequence(dna_sequence):
    logging.info("Translating DNA sequence to amino acids")

    # Warn user
    logging.warning("PLEASE NOTE: we assume that the canonical 658bp COI-5P marker has an additional base at the start")
    logging.warning("PLEASE NOTE: this first base needs to be removed to arrive at a multiple of 3 for AA translation")
    logging.warning("PLEASE NOTE: here we remove this base so that the result is 657 bases")

    # Clone and phase sequence
    cloned_seq = deepcopy(dna_sequence)
    raw = cloned_seq.seq
    cloned_seq.letter_annotations = {}
    cloned_seq.seq = raw[1:]

    # Splice codons that have missing or gap characters
    seq_str = str(cloned_seq.seq).upper()
    valid_codons = []
    for i in range(0, len(seq_str), 3):
        codon = seq_str[i:i + 3]
        if len(codon) == 3 and all(base in 'ACGT' for base in codon):
            valid_codons.append(codon)
    cloned_seq.seq = ''.join(valid_codons)

    # Do the translation
    config = Config()
    table_idx = config.get('translation_table')
    ct = CodonTable.unambiguous_dna_by_id[table_idx]
    amino_acid_sequence = cloned_seq.translate(table=ct)
    logging.debug(f"Translated sequence: {amino_acid_sequence.seq}")
    return amino_acid_sequence


def parse_fasta(file_path):
    logging.info(f"Parsing FASTA file: {file_path}")
    with open(file_path, 'r') as file:
        record = next(SeqIO.parse(file, 'fasta'))
        process_id = record.id.split('_')[0]
    logging.debug(f"Parsed process ID: {process_id}")
    logging.debug(f"Sequence length: {len(record.seq)}")
    return process_id, record


def get_stop_codons(amino_acid_sequence):
    logging.info("Getting stop codon positions")
    return [i for i, char in enumerate(amino_acid_sequence.seq) if char == '*']
