import importlib.resources
from pathlib import Path
import subprocess

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
from typing import Optional, List, Tuple
import tempfile
import re

from nbitk.Taxon import Taxon
from nbitk.config import Config
from .structural import StructuralValidator
from barcode_validator.dna_analysis_result import DNAAnalysisResult
from barcode_validator.constants import Marker


class ProteinCodingValidator(StructuralValidator):
    """
    Validator for protein-coding markers like COI-5P.

    This class extends StructuralValidator to collect protein-coding specific data,
    including sequence alignment, reading frame determination, and stop codon positions.
    These calculations form part of the structural validation process. This class
    uses nhmmer to search for gene regions by replacing gaps with Ns and aligning
    the complete sequence, then constructs the final sequence in HMM coordinate space.

    The workflow is:
    1. Remove tilde characters (preserve gaps - they represent missing gene regions)
    2. Replace gap characters with N characters to maintain sequence structure
    3. Run nhmmer on the single N-padded sequence
    4. Construct final sequence in HMM coordinate space using alignment coordinates
    5. Trim leading/trailing Ns and gaps while preserving internal Ns
    6. Analyze reading frames and stop codons

    In turn, this class exposes a single method for validation (`validate_marker_specific`).
    This method is invoked by the overall orchestrator as part of its composed validation
    process.

    Examples:
        >>> from nbitk.config import Config
        >>> config = Config()
        >>> config.load_config('/path/to/config.yaml')
        >>> validator = ProteinCodingValidator(config, '/path/to/hmm_profiles')
        >>> result = DNAAnalysisResult("foo")
        >>> validator.validate_marker_specific(SeqRecord(Seq('acgatgctacgag'),id="foo"), result, trans_table=5)

    :param config: Configuration object containing validation parameters
    """

    def __init__(self, config: Config):
        """Initialize the protein coding validator."""
        super().__init__(config)

    def validate_marker_specific(self, record: SeqRecord, result: DNAAnalysisResult) -> None:
        """
        Collect protein-coding specific measurements.
        Populates the result object with collected data.

        :param record: The DNA sequence record to validate
        :param result: Result object to populate with validation data
        """

        # Make sure that result.exp_taxon is assigned through tr.enrich_result() because this may not have been done
        # if we aren't doing taxonomic validation
        if result.exp_taxon is None:
            self.taxonomy_resolver.enrich_result(record, result)
            if result.error:
                return
        trans_table = int(self.get_translation_table(self.marker, result.exp_taxon))

        # Count original N characters in the sequence before any processing
        original_n_count = str(record.seq).upper().count('N')
        result.add_ancillary('ambig_original', str(original_n_count))
        self.logger.debug(f"Original N characters: {original_n_count}")

        # Extract and align gene region using N-padding nhmmer approach
        aligned_seq = self._align_sequence(record)
        if not aligned_seq:
            result.error = 'Gene region extraction or HMM alignment failed'
            return
        result.add_ancillary('nuc', str(aligned_seq.seq))

        # Recalculate marker length and ambiguities after HMM alignment
        result.ambiguities = self._calc_ambiguities(aligned_seq.seq)
        result.seq_length = self._calc_length(aligned_seq.seq)

        # Determine and store reading frame
        reading_frame, best_prot = self._determine_reading_frame(aligned_seq, trans_table)
        result.add_ancillary('reading_frame', str(reading_frame))
        self.logger.debug(f"Reading frame: {reading_frame}")

        # Find and store stop codon positions
        stop_positions = self._find_stop_codons(reading_frame, best_prot)
        result.stop_codons = stop_positions
        self.logger.debug(f"Stop codon positions: {stop_positions}")

    def _align_sequence(self, record: SeqRecord) -> Optional[SeqRecord]:
        """
        Process sequence using N-padding nhmmer approach.
        
        Process sequence by:
        1. Removing tilde characters (preserve gaps)
        2. Replacing gap characters with N characters
        3. Running nhmmer on the single N-padded sequence
        4. Constructing HMM space sequence from alignment results
        5. Trimming leading/trailing Ns and gaps
        
        :param record: The DNA sequence record to align
        :return: Sequence record in HMM coordinate space or None if alignment fails
        """
        # Handle hmm_profile_dir more robustly - check for type annotation issues
        hmm_file = None
        try:
            if hasattr(self, 'hmm_profile_dir') and self.hmm_profile_dir is not None:
                # Ensure we have a string, not a type annotation
                hmm_dir = str(self.hmm_profile_dir) if self.hmm_profile_dir else None
                if hmm_dir and hmm_dir != 'None':
                    hmm_file = Path(hmm_dir) / f"{self.marker.value}.hmm"
                    self.logger.debug(f"Using HMM file from profile dir: {hmm_file}")
        except (TypeError, AttributeError, ValueError) as e:
            self.logger.debug(f"Could not use hmm_profile_dir ({self.hmm_profile_dir}): {e}")
            hmm_file = None
        
        # Fallback to bundled HMM files
        if hmm_file is None or not hmm_file.exists():
            try:
                with importlib.resources.path("barcode_validator.hmm_files", f"{self.marker.value}.hmm") as path:
                    hmm_file = path
                    self.logger.debug(f"Using bundled HMM file: {hmm_file}")
            except Exception as e:
                self.logger.error(f"Could not locate bundled HMM file: {e}")
                return None
        # Verify HMM file exists
        if hmm_file is None or not hmm_file.exists():
            self.logger.error(f"HMM profile not found: {hmm_file}")
            return None

        try:
            # Step 1: Remove tilde characters (preserve gaps - they're biologically meaningful!)
            seq_str = str(record.seq).replace('~', '')
            self.logger.debug(f"After tilde removal: {seq_str[:100]}...")
            
            # Step 2: Replace all gap characters with N characters
            n_padded_seq = seq_str.replace('-', 'N')
            self.logger.debug(f"After gap-to-N replacement: {n_padded_seq[:100]}...")
            
            # Step 3: Run nhmmer on the single N-padded sequence
            nhmmer_result = self._run_nhmmer_on_sequence(n_padded_seq, record.id, hmm_file)
            if not nhmmer_result:
                self.logger.error("No significant nhmmer alignment found")
                return None
            
            # Step 4: Construct HMM space sequence using alignment coordinates
            hmm_sequence = self._construct_hmm_space_from_alignment(nhmmer_result, n_padded_seq)
            
            # Step 5: Trim leading/trailing Ns and gaps while preserving internal Ns
            trimmed_sequence = self._trim_sequence_ends(hmm_sequence)
            
            return SeqRecord(
                Seq(trimmed_sequence),
                id=record.id,
                description=record.description
            )

        except Exception as e:
            self.logger.error(f"Error in N-padding nhmmer alignment: {str(e)}")
            return None

    def _run_nhmmer_on_sequence(self, sequence: str, seq_id: str, hmm_file: Path) -> Optional[dict]:
        """
        Run nhmmer on a single sequence and return the best alignment result.
        
        :param sequence: The sequence to align
        :param seq_id: Sequence identifier
        :param hmm_file: Path to HMM file
        :return: Dictionary with alignment coordinates or None if no significant match
        """
        try:
            # Create FASTA file with the sequence
            with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as temp_input, \
                 tempfile.NamedTemporaryFile(mode='w+', suffix='.tbl', delete=False) as temp_tblout:
                
                temp_input.write(f">{seq_id}\n{sequence}\n")
                temp_input.flush()
                
                # Run nhmmer with separate tabular output file
                nhmmer_cmd = [
                    'nhmmer',
                    '--tblout', temp_tblout.name,  # Use separate file for tabular output
                    '--incE', '1e-3',
                    '--cpu', '1',
                    str(hmm_file),
                    temp_input.name
                ]
                
                self.logger.debug(f"Running nhmmer on sequence: {' '.join(nhmmer_cmd)}")
                result = subprocess.run(nhmmer_cmd, capture_output=True, text=True)
                
                if result.returncode != 0:
                    self.logger.error(f"nhmmer failed: {result.stderr}")
                    return None

                # LOG COMPLETE NHMMER OUTPUT FOR THIS SAMPLE
                self.logger.debug(f"=== COMPLETE nhmmer OUTPUT for {seq_id} ===")
                self.logger.debug(result.stdout)
                
                # Parse the clean tabular output file
                with open(temp_tblout.name, 'r') as f:
                    tabular_content = f.read()
                
                self.logger.debug(f"=== TABULAR OUTPUT for {seq_id} ===")
                self.logger.debug(tabular_content)
                
                # Parse tabular output to get best alignment
                alignment_result = self._parse_single_nhmmer_result(tabular_content, seq_id)
                
                return alignment_result

        except FileNotFoundError:
            self.logger.error("nhmmer not found in PATH")
            return None
        except Exception as e:
            self.logger.error(f"Error running nhmmer on sequence: {str(e)}")
            return None

    def _parse_single_nhmmer_result(self, tabular_content: str, seq_id: str) -> Optional[dict]:
        """
        Parse nhmmer tabular output for the best alignment result.
        
        :param tabular_content: Content of the .tbl file
        :param seq_id: Original sequence ID
        :return: Dictionary with alignment coordinates or None if no significant match
        """
        best_result = None
        best_evalue = float('inf')
        
        try:
            # Parse each line of tabular output
            for line in tabular_content.split('\n'):
                line = line.strip()
                
                # Skip comments and empty lines
                if not line or line.startswith('#'):
                    continue
                
                # Split fields
                fields = line.split()
                if len(fields) >= 15:  # Ensure we have all required fields
                    try:
                        target_name = fields[0]
                        hmm_from = int(fields[4])
                        hmm_to = int(fields[5])
                        seq_from = int(fields[6])
                        seq_to = int(fields[7])
                        evalue = float(fields[12])
                        score = float(fields[13])
                        bias = float(fields[14])
                        
                        self.logger.debug(f"Parsed alignment: {target_name} → HMM {hmm_from}-{hmm_to}, "
                                        f"Seq {seq_from}-{seq_to}, E-value: {evalue}, score: {score}")
                        
                        # Only include significant matches (should already be filtered by nhmmer --incE)
                        if evalue <= 1e-3 and target_name == seq_id:
                            # Keep the best (lowest E-value) result
                            if evalue < best_evalue:
                                best_evalue = evalue
                                best_result = {
                                    'target_name': target_name,
                                    'hmm_from': hmm_from,
                                    'hmm_to': hmm_to,
                                    'seq_from': seq_from,
                                    'seq_to': seq_to,
                                    'score': score,
                                    'bias': bias,
                                    'evalue': evalue
                                }
                                self.logger.debug(f"New best alignment: E-value={evalue}")
                        else:
                            self.logger.debug(f"Rejected alignment: E-value {evalue} > threshold 1e-3")
                    
                    except (ValueError, IndexError) as e:
                        self.logger.warning(f"Could not parse tabular line: {line} ({e})")
                        continue
                else:
                    self.logger.debug(f"Insufficient fields in tabular line: {line} (got {len(fields)}, need 15)")
            
        except Exception as e:
            self.logger.error(f"Error parsing nhmmer tabular output: {str(e)}")
        
        return best_result

    def _construct_hmm_space_from_alignment(self, alignment_result: dict, sequence: str) -> str:
        """
        Construct HMM coordinate space sequence from single alignment result.
        
        :param alignment_result: Dictionary with alignment coordinates
        :param sequence: The original N-padded sequence
        :return: Sequence string in HMM coordinate space
        """
        # Initialize 658-position HMM sequence with gaps
        hmm_length = 658  # COI-5P HMM length
        hmm_sequence = ['-'] * hmm_length
        
        self.logger.debug(f"Constructing HMM space sequence from alignment")
        
        # Extract alignment coordinates
        hmm_start = alignment_result['hmm_from'] - 1  # Convert to 0-based
        hmm_end = alignment_result['hmm_to']
        seq_start = alignment_result['seq_from'] - 1  # Convert to 0-based
        seq_end = alignment_result['seq_to']
        
        self.logger.debug(f"Placing sequence at HMM positions: {alignment_result['hmm_from']}-{alignment_result['hmm_to']} (1-based)")
        
        # Extract the aligned portion of the sequence
        aligned_sequence_portion = sequence[seq_start:seq_end]
        
        # Place sequence at HMM coordinates
        total_coverage = 0
        seq_pos = 0
        for hmm_pos in range(hmm_start, min(hmm_end, hmm_length)):
            if seq_pos < len(aligned_sequence_portion):
                hmm_sequence[hmm_pos] = aligned_sequence_portion[seq_pos]
                total_coverage += 1
                seq_pos += 1
        
        # Convert to string
        final_sequence = ''.join(hmm_sequence)
        
        # Log coverage statistics
        coverage_percentage = (total_coverage / hmm_length) * 100
        
        self.logger.debug(f"HMM space sequence construction complete:")
        self.logger.debug(f"  Total HMM length: {hmm_length} positions")
        self.logger.debug(f"  Covered positions: {total_coverage} ({coverage_percentage:.1f}%)")
        self.logger.debug(f"  Missing positions: {hmm_length - total_coverage} ({100-coverage_percentage:.1f}%)")
        self.logger.debug(f"  Final sequence length: {len(final_sequence)}")
        self.logger.debug(f"  Final sequence preview: {final_sequence[:100]}...")
        
        return final_sequence

    def _trim_sequence_ends(self, sequence: str) -> str:
        """
        Trim leading and trailing Ns and gaps while preserving internal Ns.
        
        :param sequence: Input sequence string
        :return: Trimmed sequence string
        """
        # Remove leading Ns and gaps
        start = 0
        while start < len(sequence) and sequence[start] in 'N-':
            start += 1
        
        # Remove trailing Ns and gaps
        end = len(sequence)
        while end > start and sequence[end-1] in 'N-':
            end -= 1
        
        trimmed = sequence[start:end]
        
        self.logger.debug(f"Sequence (N and gap) trimming:")
        self.logger.debug(f"  Original length: {len(sequence)}")
        self.logger.debug(f"  Trimmed length: {len(trimmed)}")
        self.logger.debug(f"  Removed from start: {start}")
        self.logger.debug(f"  Removed from end: {len(sequence) - end}")
        
        return trimmed

    @staticmethod
    def _get_complete_codons(seq, offset):
        """
        Extract complete codons from sequence starting at given offset.
        
        :param seq: Input sequence
        :param offset: Reading frame offset (0, 1, or 2)
        :return: String of complete codons (no gaps or Ns)
        """
        complete_codons = ""
        # Loop over the sequence in steps of 3, starting from the given offset.
        for i in range(offset, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            # Only add the codon if it is complete and contains no gap characters or Ns
            if '-' not in codon and 'N' not in codon:
                complete_codons += codon
        return complete_codons

    def _determine_reading_frame(self, aligned_seq: SeqRecord, trans_table: int) -> Tuple[int, SeqRecord]:
        """
        Determine the reading frame from HMM-aligned sequence.

        :param aligned_seq: HMM-aligned sequence
        :param trans_table: Translation table to use
        :return: Reading frame (0, 1, or 2) and best protein sequence
        """
        raw_seq = aligned_seq.seq
        best_frame = 0
        best_prot = ""
        min_stops = float('inf')

        self.logger.debug(f"Analysing reading frames")
        
        # Analyze all 3 reading frames
        frame_results = []
        for frame in range(3):
            coding_seq = Seq("".join(self._get_complete_codons(raw_seq, frame)))
            protein = coding_seq.translate(table=trans_table)
            stops = protein.count('*')
            
            frame_results.append({
                'frame': frame,
                'stops': stops,
                'protein_length': len(protein),
                'coding_length': len(coding_seq)
            })
            
            self.logger.debug(f"Frame {frame}: {stops} stop codons, {len(protein)} amino acids, {len(coding_seq)} coding nucleotides")

            if stops < min_stops:
                min_stops = stops
                best_frame = frame
                best_prot = protein

        # Log detailed results for best frame
        self.logger.debug(f"=== READING FRAME ANALYSIS RESULTS ===")
        self.logger.debug(f"Best reading frame: {best_frame}")
        self.logger.debug(f"Stop codons in best frame: {min_stops}")
        self.logger.debug(f"Protein length: {len(best_prot)} amino acids")
        
        # LOG COMPLETE AMINO ACID SEQUENCE FOR BEST FRAME
        self.logger.debug(f"=== TRANSLATED AMINO ACID SEQUENCE (Frame {best_frame}) ===")
        self.logger.debug(f"Translation table: {trans_table}")
        self.logger.debug(f"Amino acid sequence: {str(best_prot)}")
        
        return best_frame, best_prot

    def _find_stop_codons(self, reading_frame, protein) -> List[int]:
        """
        Find positions of stop codons in the sequence.

        :param reading_frame: The reading frame of the protein sequence (0, 1, or 2)
        :param protein: The protein sequence record
        :return: List of stop codon positions in nucleotide coordinates
        """
        stop_positions = []
        
        self.logger.debug(f"Searching for stop codons in protein sequence (length: {len(protein)})")
        
        # Find all stop codons (excluding the natural terminal stop)
        for i, aa in enumerate(protein[:-1]):  # Exclude last amino acid (natural terminal)
            if aa == '*':
                # Convert protein position to nucleotide position in the reading frame
                nuc_pos = (i * 3) + reading_frame
                stop_positions.append(nuc_pos)
                self.logger.debug(f"Stop codon found: protein position {i} → nucleotide position {nuc_pos}")

        self.logger.debug(f"=== STOP CODON ANALYSIS ===")
        self.logger.debug(f"Reading frame: {reading_frame}")
        self.logger.debug(f"Total stop codons found: {len(stop_positions)}")
        self.logger.debug(f"Stop codon positions (nucleotide coordinates): {stop_positions}")
        self.logger.debug(f"=== END STOP CODON ANALYSIS ===")
        
        return stop_positions

    def get_translation_table(self, marker: Marker, taxon: Taxon) -> int:
        """
        Determine the appropriate translation table based on marker and taxonomy.

        :param marker: The genetic Marker (enum) being analyzed
        :param taxon: The NCBI taxon object representing at the highest a family in the taxonomy.
        :return: Translation table index (int) for use with Biopython
        """
        taxonomy_dict = self.taxonomy_resolver.get_lineage_dict(taxon)

        if marker in [Marker.MATK, Marker.RBCL]:
            if taxonomy_dict.get('family') == 'Balanophoraceae':
                return 32
            return 11

        elif marker == Marker.COI_5P:
            phylum = taxonomy_dict.get('phylum')
            tax_class = taxonomy_dict.get('class')
            family = taxonomy_dict.get('family')

            if phylum == 'Chordata':
                if tax_class == 'Ascidiacea':
                    return 13

                # This list was gradually arrived at through the discovery of mismatches between BOLD and NCBI
                # taxonomies at class level. Desparately needs a generic solution.
                elif tax_class in ['Actinopteri', 'Actinopterygii', 'Myxini', 'Elasmobranchii', 'Petromyzonti', 'Amphibia', 'Mammalia', 'Aves', 'Reptilia']:
                    return 2

            elif phylum == 'Hemichordata':
                if family == 'Cephalodiscidae':
                    return 33
                elif family == 'Rhabdopleuridae':
                    return 24

            elif phylum in ['Echinodermata', 'Platyhelminthes']:
                return 9

            # Default invertebrate mitochondrial code for other invertebrates
            if phylum != 'Chordata':
                return 5

        # Fall back to standard code if no specific rules match
        return 1

    @staticmethod
    def requires_resolver() -> bool:
        return True

    @staticmethod
    def requires_marker() -> bool:
        return True

    @staticmethod
    def requires_hmmalign() -> bool:
        return False  # We no longer use hmmalign (replaced with nhmmer)

    @staticmethod
    def requires_nhmmer() -> bool:
        """nhmmer is required for gene region detection."""
        return True