import logging
import tempfile
import subprocess
import warnings
import json
from copy import deepcopy
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


class SequenceHandler:

    @classmethod
    def align_to_hmm(cls, sequence, hmm_file):
        """
        Align a sequence to an HMM using hmmalign. The location of the HMM file is specified in the
        configuration file as 'hmm_file'. hmmalign is run with the '--trim' option to remove any
        leading or trailing ends of the input sequence that fall outside the HMM. The output is
        parsed as a Stockholm format file and returned as a SeqRecord object, i.e. it has the sequence
        annotations (which include per-residue posterior probabilities) preserved.
        :param sequence: A BioPython SeqRecord object
        :param hmm_file: Location of an HMM file
        :return: A BioPython SeqRecord object containing the aligned sequence
        """
        logging.info("Aligning sequence to HMM")
        if len(sequence.seq) == 0:
            logging.warning("Empty sequence provided for alignment")
            return None

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

    @classmethod
    def marker_seqlength(cls, sequence):
        """
        Calculate the length of the sequence within the marker region,
        excluding gaps, but including Ns and IUAPC ambiguity codes. Both '-' and '~' are considered gaps.
        :param sequence: A BioPython SeqRecord object
        :return: The length of the sequence within the marker region
        """
        logging.info("Calculating non-missing bases within marker region")

        # Operate on the cloned sequence string directly
        raw_seq = str(sequence.seq)
        raw_seq = raw_seq.replace('-', '')
        raw_seq = raw_seq.replace('~', '')

        logging.debug(f"Within marker sequence length: {len(raw_seq)}")
        return len(raw_seq)

    @classmethod
    def num_ambiguous(cls, sequence):
        """
        Calculate the number of ambiguous bases in a sequence. This is a count of all symbols that are not
        acgtACGT-~ or, in other words, all IUPAC single character ambiguity letters including n/N.

        :param sequence: A BioPython SeqRecord object
        :return: The number of ambiguous bases in the sequence
        """
        logging.info("Calculating number of ambiguous bases in sequence")
        return len([base for base in sequence.seq if base not in 'acgtACGT-~'])

    @classmethod
    def unalign_sequence(cls, sequence):
        """
        Remove gaps, i.e. '-' symbols, from an aligned sequence.
        Returns a new instance of the sequence with gaps removed.
        Note that, in the case of sequences aligned with hmmalign,
        this probably has undefined (i.e. bad) consequences for the per-residue
        posterior probability vector.

        :param sequence: A BioPython SeqRecord object
        :return: A new SeqRecord object with gaps removed
        """
        logging.info("Removing gaps from aligned sequence")
        if isinstance(sequence, SeqRecord):
            # Convert Seq to string, remove gaps, then convert back to Seq
            unaligned_sequence = str(sequence.seq).replace('-', '').replace('~', '')
            return SeqRecord(
                Seq(unaligned_sequence),
                id=sequence.id,
                name=sequence.name,
                description=sequence.description
            )
        elif isinstance(sequence, Seq):
            # If it's just a Seq object, convert to string, remove gaps, then back to Seq
            return Seq(str(sequence).replace('-', '').replace('~', ''))
        elif isinstance(sequence, str):
            # If it's a string, just remove the gaps
            return sequence.replace('-', '').replace('~', '')
        else:
            raise TypeError(f"Unexpected type for sequence: {type(sequence)}")

    @classmethod
    def translate_sequence(cls, dna_sequence, table_idx):
        """
        Translate a DNA sequence to amino acids using the translation table specified in the configuration file.
        The translation table is specified as an integer in the configuration file, which corresponds to the
        NCBI translation table ID. The translated sequence is returned as an amino acid SeqRecord object. Note that
        the first base of the DNA sequence is removed to ensure that the sequence length is a multiple of 3 (i.e.
        the canonical 658bp COI-5P marker is out of phase by one base).
        :param dna_sequence: A BioPython SeqRecord object containing a DNA sequence
        :param table_idx: Translation table index
        :return: A BioPython SeqRecord object containing the translated amino acid sequence
        """
        logging.info("Translating DNA sequence to amino acids")

        # Warn user
        logging.warning("NOTE: we assume that the canonical 658bp COI-5P marker has an additional base at the start")
        logging.warning("NOTE: this first base needs to be removed to arrive at a multiple of 3 for AA translation")
        logging.warning("NOTE: here we remove this base so that the result is 657 bases")

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
        ct = CodonTable.unambiguous_dna_by_id[table_idx]
        amino_acid_sequence = valid_seq.translate(table=ct)

        # Create a new SeqRecord for the amino acid sequence
        amino_acid_record = SeqRecord(amino_acid_sequence,
                                      id=cloned_seq.id,
                                      name=cloned_seq.name,
                                      description="translated sequence")

        logging.debug(f"Translated sequence: {amino_acid_record.seq}")
        return amino_acid_record

    @classmethod
    def parse_fasta(cls, file_path):
        """
        Parse a FASTA file and yield the process ID, sequence record, and optional JSON configuration for each entry.
        The process ID is the first part of the sequence ID, which is assumed to be separated by an underscore.
        Any JSON configuration data after the first '{' in the header is parsed and returned.

        :param file_path: Local path to the FASTA file
        :yield: A tuple containing the process ID, the sequence record, and the JSON configuration (or None) for each entry
        """
        logging.info(f"Parsing FASTA file: {file_path}")
        with open(file_path, 'r') as file:
            for record in SeqIO.parse(file, 'fasta'):
                process_id = record.id.split('_')[0]

                # Attempt to parse JSON from the description
                json_config = None
                json_start = record.description.find('{')
                if json_start != -1:
                    try:
                        json_str = record.description[json_start:]
                        json_config = json.loads(json_str)
                        # Remove the JSON part from the description
                        record.description = record.description[:json_start].strip()
                    except json.JSONDecodeError as e:
                        logging.warning(f"Failed to parse JSON for {process_id}: {e}")

                record.id = process_id
                logging.debug(f"Parsed process ID: {process_id}")
                logging.debug(f"Sequence length: {len(record.seq)}")
                logging.debug(f"JSON config: {json_config}")

                yield process_id, cls.unalign_sequence(record), json_config

    @classmethod
    def get_stop_codons(cls, amino_acid_sequence):
        """
        Get the positions of stop codons in an amino acid sequence, which are represented by the '*' character.
        :param amino_acid_sequence: A BioPython SeqRecord object containing an amino acid sequence
        :return: A list of integers representing the positions of stop codons in the sequence
        """
        logging.info("Getting stop codon positions")
        return [i for i, char in enumerate(amino_acid_sequence.seq) if char == '*']
