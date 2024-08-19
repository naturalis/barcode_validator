import logging
import tempfile
import subprocess
import warnings
from copy import deepcopy
from barcode_validator.config import Config
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from Bio import BiopythonDeprecationWarning
warnings.filterwarnings("error", category=BiopythonDeprecationWarning)

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
    """
    Calculate the length of the sequence within the marker region,
    excluding gaps and ambiguous bases.
    :param sequence: A BioPython SeqRecord object
    :return: The length of the sequence within the marker region
    """
    logging.info("Calculating non-missing bases within marker region")

    # Operate on the cloned sequence string directly
    raw_seq = str(sequence.seq)
    raw_seq = raw_seq.replace('-', '')
    raw_seq = raw_seq.replace('N', '')
    raw_seq = raw_seq.replace('n', '')

    logging.debug(f"Within marker sequence length: {len(raw_seq)}")
    return len(raw_seq)


def num_ambiguous(sequence):
    logging.info("Calculating number of ambiguous bases in sequence")

    # Clone the sequence, count the number of ambiguous bases
    cloned_seq = deepcopy(sequence)
    cloned_seq.letter_annotations = {}
    return len([base for base in cloned_seq.seq if base not in 'acgtnACGTN-~?'])


#def unalign_sequence(sequence):
#    logging.info("Removing gaps from aligned sequence")
#    unaligned_sequence = sequence.seq.replace('-', '')
#    unaligned_sequence = unaligned_sequence.replace('~', '')
#    sequence.seq = unaligned_sequence
#    logging.debug(f"Unaligned sequence length: {len(unaligned_sequence)}")
#    return sequence


def unalign_sequence(sequence):
    """
    Remove gaps from an aligned sequence.

    :param sequence: A BioPython SeqRecord object
    :return: A new SeqRecord object with gaps removed
    """
    logging.info("Removing gaps from aligned sequence")
    if isinstance(sequence, SeqRecord):
        # Convert Seq to string, remove gaps, then convert back to Seq
        unaligned_sequence = str(sequence.seq).replace('-', '')
        return SeqRecord(
            Seq(unaligned_sequence),
            id=sequence.id,
            name=sequence.name,
            description=sequence.description
        )
    elif isinstance(sequence, Seq):
        # If it's just a Seq object, convert to string, remove gaps, then back to Seq
        return Seq(str(sequence).replace('-', ''))
    elif isinstance(sequence, str):
        # If it's a string, just remove the gaps
        return sequence.replace('-', '')
    else:
        raise TypeError(f"Unexpected type for sequence: {type(sequence)}")


def translate_sequence(dna_sequence):
    logging.info("Translating DNA sequence to amino acids")

    # Warn user
    logging.warning("PLEASE NOTE: we assume that the canonical 658bp COI-5P marker has an additional base at the start")
    logging.warning("PLEASE NOTE: this first base needs to be removed to arrive at a multiple of 3 for AA translation")
    logging.warning("PLEASE NOTE: here we remove this base so that the result is 657 bases")

    # Clone and phase sequence by starting from the second base (i.e. index 1 in 0-based indexing)
    cloned_seq = deepcopy(dna_sequence)
    raw = cloned_seq.seq
    cloned_seq.letter_annotations = {}
    cloned_seq.seq = raw[1:]

    # Skip codons that have missing, ambiguous or gap characters
    seq_str = str(cloned_seq.seq).upper()
    valid_codons = []
    for i in range(0, len(seq_str), 3):
        codon = seq_str[i:i + 3]
        if len(codon) == 3 and all(base in 'ACGT' for base in codon):
            valid_codons.append(codon)

    # Create a new Seq object from the valid codons
    valid_seq = Seq(''.join(valid_codons))

    # Do the translation
    config = Config()
    table_idx = config.get('translation_table')
    ct = CodonTable.unambiguous_dna_by_id[table_idx]
    amino_acid_sequence = valid_seq.translate(table=ct)

    # Create a new SeqRecord for the amino acid sequence
    amino_acid_record = SeqRecord(amino_acid_sequence,
                                  id=cloned_seq.id,
                                  name=cloned_seq.name,
                                  description="translated sequence")

    logging.debug(f"Translated sequence: {amino_acid_record.seq}")
    return amino_acid_record


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
