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
    uses a fragment-based nhmmer approach to search for gene regions within longer 
    barcode sequences, then constructs the final sequence in HMM coordinate space.

    The workflow is:
    1. Remove tilde characters (preserve gaps - they represent missing gene regions)
    2. Split sequence on gaps to create biologically meaningful fragments
    3. Run nhmmer on multi-FASTA of fragments to find HMM matches
    4. Construct final sequence in HMM coordinate space
    5. Analyze reading frames and stop codons

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

        # Extract and align gene region using fragment-based nhmmer approach
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
        Process sequence using fragment-based nhmmer approach.
        
        Process sequence by:
        1. Removing tilde characters (preserve gaps)
        2. Splitting on gaps to create fragments
        3. Running nhmmer on multi-FASTA of fragments
        4. Constructing HMM space sequence from results
        
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
            
            # Step 2: Split sequence on gaps to create fragments
            fragments = self._split_sequence_on_gaps(seq_str)
            if not fragments:
                self.logger.error("No fragments found after gap splitting")
                return None
                
            self.logger.debug(f"Created {len(fragments)} fragments from gap splitting")
            
            # Step 3: Run nhmmer on multi-FASTA of fragments
            fragment_matches = self._search_fragments_with_nhmmer(fragments, record.id, hmm_file)
            if not fragment_matches:
                self.logger.error("No fragments matched HMM profile")
                return None
            
            # Step 4: Construct HMM space sequence
            hmm_sequence = self._construct_hmm_space_from_fragments(fragment_matches)
            
            return SeqRecord(
                Seq(hmm_sequence),
                id=record.id,
                description=record.description
            )

        except Exception as e:
            self.logger.error(f"Error in fragment-based nhmmer alignment: {str(e)}")
            return None

    def _split_sequence_on_gaps(self, seq_str: str) -> List[Tuple[str, int, int]]:
        """
        Split sequence on gap characters and return fragments with their positions.
        
        :param seq_str: Input sequence string
        :return: List of tuples (fragment_sequence, start_position, end_position)
        """
        fragments_info = []
        current_fragment = ""
        start_pos = 0
        
        for i, char in enumerate(seq_str):
            if char == '-':
                if current_fragment:
                    # End current fragment
                    fragments_info.append((current_fragment, start_pos, i))
                    current_fragment = ""
                # Find start of next fragment (skip consecutive gaps)
                start_pos = i + 1
            else:
                if not current_fragment:
                    start_pos = i
                current_fragment += char
        
        # Add final fragment if exists
        if current_fragment:
            fragments_info.append((current_fragment, start_pos, len(seq_str)))
        
        # Log fragment information
        for i, (fragment, start_pos, end_pos) in enumerate(fragments_info):
            self.logger.debug(f"Fragment {i+1}: positions {start_pos}-{end_pos}, length {len(fragment)}")
            self.logger.debug(f"Fragment {i+1} sequence: {fragment[:50]}...")
        
        return fragments_info

    def _search_fragments_with_nhmmer(self, fragments: List[Tuple[str, int, int]], seq_id: str, hmm_file: Path) -> List[dict]:
        """
        Run nhmmer on multi-FASTA of fragments to find HMM matches.
        
        :param fragments: List of (fragment_sequence, start_pos, end_pos) tuples
        :param seq_id: Original sequence identifier
        :param hmm_file: Path to HMM file
        :return: List of fragment matches with HMM coordinates
        """
        try:
            # Create multi-FASTA file with all fragments
            with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as temp_input, \
                 tempfile.NamedTemporaryFile(mode='w+', suffix='.tbl', delete=False) as temp_tblout:
                
                for i, (fragment_seq, start_pos, end_pos) in enumerate(fragments):
                    fragment_id = f"{seq_id}_fragment_{i+1}"
                    temp_input.write(f">{fragment_id}\n{fragment_seq}\n")
                
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
                
                self.logger.debug(f"Running nhmmer on {len(fragments)} fragments: {' '.join(nhmmer_cmd)}")
                result = subprocess.run(nhmmer_cmd, capture_output=True, text=True)
                
                if result.returncode != 0:
                    self.logger.error(f"nhmmer failed: {result.stderr}")
                    return []

                # LOG COMPLETE NHMMER OUTPUT FOR THIS SAMPLE
                self.logger.debug(f"=== COMPLETE nhmmer OUTPUT for {seq_id} ===")
                self.logger.debug(result.stdout)
                self.logger.debug(f"=== END nhmmer OUTPUT for {seq_id} ===")
                
                # Parse the clean tabular output file
                with open(temp_tblout.name, 'r') as f:
                    tabular_content = f.read()
                
                self.logger.debug(f"=== TABULAR OUTPUT for {seq_id} ===")
                self.logger.debug(tabular_content)
                self.logger.debug(f"=== END TABULAR OUTPUT for {seq_id} ===")
                
                # Parse tabular output to get fragment matches
                fragment_matches = self._parse_nhmmer_tabular_output(tabular_content, fragments, seq_id)
                
                self.logger.debug(f"Found {len(fragment_matches)} fragment matches:")
                for match in fragment_matches:
                    self.logger.debug(f"  {match['fragment_id']}: HMM {match['hmm_from']}-{match['hmm_to']}, E-value: {match['evalue']}")
                
                return fragment_matches

        except FileNotFoundError:
            self.logger.error("nhmmer not found in PATH")
            return []
        except Exception as e:
            self.logger.error(f"Error running nhmmer on fragments: {str(e)}")
            return []

    def _parse_nhmmer_tabular_output(self, tabular_content: str, fragments: List[Tuple[str, int, int]], seq_id: str) -> List[dict]:
        """
        Parse nhmmer tabular output (.tbl file) for fragment matches.
        
        :param tabular_content: Content of the .tbl file
        :param fragments: Original fragment list for reference
        :param seq_id: Original sequence ID
        :return: List of fragment matches with HMM coordinates
        """
        fragment_matches = []
        
        try:
            # Create fragment lookup for easy reference
            fragment_lookup = {}
            for i, (fragment_seq, orig_start, orig_end) in enumerate(fragments):
                fragment_id = f"{seq_id}_fragment_{i+1}"
                fragment_lookup[fragment_id] = {
                    'index': i,
                    'seq': fragment_seq,
                    'orig_start': orig_start,
                    'orig_end': orig_end
                }
            
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
                        evalue = float(fields[12])
                        score = float(fields[13])
                        bias = float(fields[14])
                        
                        self.logger.debug(f"Parsed tabular line: {target_name} → HMM {hmm_from}-{hmm_to}, E-value: {evalue}, score: {score}")
                        
                        # Look up fragment info
                        if target_name in fragment_lookup:
                            frag_info = fragment_lookup[target_name]
                            
                            fragment_match = {
                                'fragment_id': target_name,
                                'fragment_index': frag_info['index'],
                                'fragment_seq': frag_info['seq'],
                                'original_start': frag_info['orig_start'],
                                'original_end': frag_info['orig_end'],
                                'score': score,
                                'bias': bias,
                                'evalue': evalue,
                                'hmm_from': hmm_from,
                                'hmm_to': hmm_to
                            }
                            
                            # Only include significant matches (should already be filtered by nhmmer --incE)
                            if evalue <= 1e-3:
                                fragment_matches.append(fragment_match)
                                self.logger.debug(f"Added significant match: {target_name} → HMM {hmm_from}-{hmm_to}, E-value={evalue}")
                            else:
                                self.logger.debug(f"Rejected {target_name}: E-value {evalue} > threshold 1e-3")
                        else:
                            self.logger.warning(f"Unknown target name in tabular output: {target_name}")
                    
                    except (ValueError, IndexError) as e:
                        self.logger.warning(f"Could not parse tabular line: {line} ({e})")
                        continue
                else:
                    self.logger.debug(f"Insufficient fields in tabular line: {line} (got {len(fields)}, need 15)")
            
            # Sort by HMM position for logical ordering
            fragment_matches.sort(key=lambda x: x['hmm_from'])
            
        except Exception as e:
            self.logger.error(f"Error parsing nhmmer tabular output: {str(e)}")
        
        return fragment_matches


    def _construct_hmm_space_from_fragments(self, fragment_matches: List[dict]) -> str:
        """
        Construct HMM coordinate space sequence from fragment matches.
        
        :param fragment_matches: List of fragment matches with HMM coordinates
        :return: Sequence string in HMM coordinate space
        """
        # Initialize 658-position HMM sequence with gaps
        hmm_length = 658  # COI-5P HMM length
        hmm_sequence = ['-'] * hmm_length
        
        self.logger.debug(f"Constructing HMM space sequence from {len(fragment_matches)} fragment matches")
        
        total_coverage = 0
        
        # Place each fragment at its HMM coordinates
        for match in fragment_matches:
            fragment_seq = match['fragment_seq']
            hmm_start = match['hmm_from'] - 1  # Convert to 0-based (1-based HMM → 0-based Python)
            hmm_end = match['hmm_to']
            
            # Log the coordinate conversion for clarity
            self.logger.debug(f"Placing {match['fragment_id']} ({len(fragment_seq)}bp):")
            self.logger.debug(f"  nhmmer HMM positions: {match['hmm_from']}-{match['hmm_to']} (1-based)")
            self.logger.debug(f"  Python array indices: {hmm_start}-{hmm_end-1} (0-based)")
            
            # Place fragment sequence at HMM positions
            seq_pos = 0
            for hmm_pos in range(hmm_start, min(hmm_end, hmm_length)):
                if seq_pos < len(fragment_seq):
                    if hmm_sequence[hmm_pos] == '-':  # Only place if position is empty
                        hmm_sequence[hmm_pos] = fragment_seq[seq_pos]
                        total_coverage += 1
                    else:
                        self.logger.warning(f"Overlap detected at HMM position {hmm_pos+1} (1-based)")
                    seq_pos += 1
            
            self.logger.debug(f"Successfully placed {seq_pos} nucleotides from {match['fragment_id']}")
        
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

    @staticmethod
    def _get_complete_codons(seq, offset):
        """
        Extract complete codons from sequence starting at given offset.
        
        :param seq: Input sequence
        :param offset: Reading frame offset (0, 1, or 2)
        :return: String of complete codons (no gaps)
        """
        complete_codons = ""
        # Loop over the sequence in steps of 3, starting from the given offset.
        for i in range(offset, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            # Only add the codon if it is complete and contains no gap characters.
            if '-' not in codon:
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

        self.logger.debug(f"Analyzing reading frames for sequence length: {len(raw_seq)}")
        
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
        self.logger.debug(f"=== END AMINO ACID SEQUENCE ===")
        
        # Log comparison of all frames
        self.logger.debug(f"=== FRAME COMPARISON ===")
        for result in frame_results:
            self.logger.debug(f"Frame {result['frame']}: {result['stops']} stops, {result['protein_length']} aa, {result['coding_length']} coding bp")
        self.logger.debug(f"=== END FRAME COMPARISON ===")
        
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
                elif tax_class in ['Actinopteri', 'Amphibia', 'Mammalia', 'Aves', 'Reptilia']:
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
        return False  # We no longer use hmmalign as primary method

    @staticmethod
    def requires_nhmmer() -> bool:
        """nhmmer is required for gene region detection."""
        return True