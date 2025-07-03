import requests
import time
import json
from typing import Optional, Dict, Any, Union, Set
from pathlib import Path
import re
from io import StringIO

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from nbitk.Taxon import Taxon
from nbitk.config import Config
from barcode_validator.constants import TaxonomicRank

from barcode_validator.idservices.idservice import IDService


class BOLD(IDService):
    """BOLD Systems implementation of the IDService interface."""

    def __init__(self, config: Config, base_url: str = "https://id.boldsystems.org"):
        """
        Initialize the BOLD ID service.

        :param config: Configuration object
        :param base_url: Base URL for the BOLD ID service
        """
        super().__init__(config)
        self.base_url = base_url.rstrip('/')
        self.session = requests.Session()

        # Set headers to mimic browser behavior
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (compatible; BOLD-ID-Python-Client/1.0)',
            'Accept': 'application/json, text/plain, */*',
            'Accept-Language': 'en-US,en;q=0.9',
        })

        # Default search parameters
        self.default_database = "all.tax-derep"  # Animal Library (Public + Private)
        self.default_min_identity = 0.94
        self.default_min_overlap = 100
        self.default_max_hits = 25
        self.default_order = 3
        self.poll_interval = 5
        self.max_wait_time = 300

    def identify_record(self, record: SeqRecord, level: TaxonomicRank = TaxonomicRank.FAMILY, extent: Taxon = None) -> Set[Taxon]:
        """
        Identify the taxonomic classification of a sequence record using BOLD ID service.

        :param record: A Bio.SeqRecord object containing the sequence to identify
        :param level: The taxonomic rank at which to return results (default: FAMILY)
        :param extent: A Taxon object representing the extent of the search (currently unused)
        :return: Set of Taxon objects representing possible taxonomic classifications at the specified level
        """
        self.logger.info(f"Identifying sequence record: {record.id}")

        try:
            # Convert SeqRecord to FASTA format
            fasta_content = self._seqrecord_to_fasta(record)

            # Submit to BOLD ID service
            results = self._identify_sequences(fasta_content)

            # Parse results and extract taxa at the specified level
            taxa_set = self._parse_results_to_taxa(results, level)

            self.logger.info(f"Found {len(taxa_set)} taxa at {level.value} level for record {record.id}")
            return taxa_set

        except Exception as e:
            self.logger.error(f"Error identifying record {record.id}: {str(e)}")
            return set()

    def _seqrecord_to_fasta(self, record: SeqRecord) -> str:
        """
        Convert a SeqRecord to FASTA format string.

        :param record: SeqRecord to convert
        :return: FASTA formatted string
        """
        output = StringIO()
        SeqIO.write(record, output, "fasta")
        return output.getvalue()

    def _identify_sequences(self, fasta_content: str) -> Dict[str, Any]:
        """
        Submit sequences to BOLD ID service and return results.

        :param fasta_content: FASTA formatted sequence string
        :return: Dictionary containing identification results
        """
        # Validate FASTA format
        if not self._is_valid_fasta(fasta_content):
            raise ValueError("Invalid FASTA format")

        # Submit sequences
        sub_id = self._submit_sequences(fasta_content)
        self.logger.info(f"BOLD submission ID: {sub_id}")

        # Poll for results
        return self._wait_for_results(sub_id)

    def _submit_sequences(self, fasta_content: str) -> str:
        """
        Submit sequences to BOLD ID service.

        :param fasta_content: FASTA formatted sequence string
        :return: Submission ID for tracking the request
        """

        # Prepare form data
        files = {
            'file': ('sequences.fas', fasta_content, 'text/plain')
        }

        params = {
            'db': self.default_database,
            'mi': str(self.default_min_identity),
            'mo': str(self.default_min_overlap),
            'maxh': str(self.default_max_hits),
            'order': str(self.default_order)
        }

        # Submit to submission endpoint
        url = f"{self.base_url}/submission"
        response = self.session.post(url, files=files, params=params)

        if response.status_code != 200:
            try:
                error_detail = response.json().get('detail', 'Unknown error')
            except:
                error_detail = f"HTTP {response.status_code}: {response.text}"
            raise requests.RequestException(f"BOLD submission failed: {error_detail}")

        result = response.json()
        return result['sub_id']

    def _wait_for_results(self, sub_id: str) -> Dict[str, Any]:
        """
        Poll for results until completion or timeout.

        :param sub_id: Submission ID to track
        :return: Dictionary containing the analysis results
        """

        start_time = time.time()

        while time.time() - start_time < self.max_wait_time:
            # Check status
            status_url = f"{self.base_url}/status/{sub_id}"
            response = self.session.get(status_url)

            if response.status_code == 200:
                status_data = response.json()

                if status_data.get('status') == 'completed':
                    # Get results
                    results_url = f"{self.base_url}/results/{sub_id}"
                    results_response = self.session.get(results_url)

                    if results_response.status_code == 200:
                        return results_response.json()
                    else:
                        raise requests.RequestException(f"Failed to retrieve BOLD results: {results_response.text}")

                elif status_data.get('status') == 'failed':
                    error_msg = status_data.get('error', 'Unknown error')
                    raise requests.RequestException(f"BOLD analysis failed: {error_msg}")

                # Still processing, wait and try again
                self.logger.debug(f"BOLD status: {status_data.get('status', 'processing')}...")
                time.sleep(self.poll_interval)

            else:
                self.logger.warning(f"BOLD status check failed (HTTP {response.status_code}), retrying...")
                time.sleep(self.poll_interval)

        raise TimeoutError(f"BOLD results not available within {self.max_wait_time} seconds")

    def _parse_results_to_taxa(self, results: Dict[str, Any], level: TaxonomicRank) -> Set[Taxon]:
        """
        Parse BOLD ID results and extract taxa at the specified taxonomic level.

        :param results: Results dictionary from BOLD ID service
        :param level: Taxonomic rank to extract
        :return: Set of Taxon objects at the specified level
        """
        taxa_set = set()

        # The exact structure of BOLD results may vary, so this is a best-effort implementation
        # Based on typical BOLD API responses, results usually contain matches with taxonomic info

        try:
            # Handle different possible result structures
            matches = results.get('matches', [])
            if not matches and 'results' in results:
                matches = results['results']

            for match in matches:
                # Extract taxonomic information from the match
                taxonomy = self._extract_taxonomy_from_match(match)

                if taxonomy:
                    # Get the taxon at the specified level
                    taxon_at_level = self._get_taxon_at_level(taxonomy, level)
                    if taxon_at_level:
                        taxa_set.add(taxon_at_level)

        except Exception as e:
            self.logger.warning(f"Error parsing BOLD results: {str(e)}")

        return taxa_set

    def _extract_taxonomy_from_match(self, match: Dict[str, Any]) -> Optional[Dict[str, str]]:
        """
        Extract taxonomic information from a BOLD match result.

        :param match: Single match result from BOLD
        :return: Dictionary of taxonomic classifications or None if not found
        """
        taxonomy = {}

        # Common fields in BOLD results
        taxonomic_fields = [
            'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species',
            'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'
        ]

        for field in taxonomic_fields:
            if field in match and match[field]:
                taxonomy[field.lower()] = match[field]

        # Also check for nested taxonomy objects
        if 'taxonomy' in match:
            tax_obj = match['taxonomy']
            for field in taxonomic_fields:
                if field in tax_obj and tax_obj[field]:
                    taxonomy[field.lower()] = tax_obj[field]

        return taxonomy if taxonomy else None

    def _get_taxon_at_level(self, taxonomy: Dict[str, str], level: TaxonomicRank) -> Optional[Taxon]:
        """
        Extract a Taxon object at the specified taxonomic level.

        :param taxonomy: Dictionary of taxonomic classifications
        :param level: Desired taxonomic rank
        :return: Taxon object at the specified level, or None if not found
        """
        level_name = level.value.lower()

        if level_name in taxonomy:
            taxon_name = taxonomy[level_name]

            if self.taxonomy_resolver:
                # Use the taxonomy resolver if available
                try:
                    return self.taxonomy_resolver.resolve_name(taxon_name, level)
                except Exception as e:
                    self.logger.warning(f"Could not resolve taxon '{taxon_name}' at level '{level_name}': {str(e)}")

            # Create a basic Taxon object without resolver
            # Note: This may need adjustment based on your Taxon class implementation
            return Taxon(name=taxon_name, taxonomic_rank=level.value)

        return None

    def _is_valid_fasta(self, content: str) -> bool:
        """
        Check if content is valid FASTA format.

        :param content: String content to validate
        :return: True if valid FASTA format, False otherwise
        """
        content = content.strip()
        if not content:
            return False

        lines = content.split('\n')
        has_header = False
        has_sequence = False

        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                has_header = True
            elif has_header and re.match(r'^[ATCGNRYSWKMBDHV\-]+$', line.upper()):
                has_sequence = True

        return has_header and has_sequence

    def set_search_parameters(self, database: str = None, min_identity: float = None,
                              min_overlap: int = None, max_hits: int = None):
        """
        Set search parameters for BOLD ID service.

        :param database: Database to search against
        :param min_identity: Minimum similarity threshold (0.70-1.00)
        :param min_overlap: Minimum overlap length (1-1000)
        :param max_hits: Maximum hits per sequence (1-1000)
        """
        if database is not None:
            self.default_database = database
        if min_identity is not None:
            self.default_min_identity = min_identity
        if min_overlap is not None:
            self.default_min_overlap = min_overlap
        if max_hits is not None:
            self.default_max_hits = max_hits

    def get_available_databases(self) -> Dict[str, str]:
        """
        Get information about available BOLD databases.

        :return: Dictionary mapping database codes to descriptions
        """
        return {
            "public.tax-derep": "Animal Library (Public) - Non-redundant COI sequences",
            "species": "Animal Species-Level Library (Public + Private) - COI sequences with species designation",
            "all.tax-derep": "Animal Library (Public + Private) - Non-redundant COI sequences",
            "DS-CANREF22": "Validated Canadian Arthropod Library - Validated Canadian arthropod records",
            "public.plants": "Plant Library (Public) - Non-redundant plant rbcL, matK, and ITS sequences",
            "public.fungi": "Fungi Library (Public) - Non-redundant fungal ITS and 18S sequences",
            "all.animal-alt": "Animal Secondary Markers (Public) - Non-redundant 18S and 12S sequences",
            "DS-IUCNPUB": "Validated Animal Red List Library - Validated animal Red List species"
        }

    @staticmethod
    def requires_resolver():
        """
        Return True if this service benefits from a taxonomy resolver.

        :return: True, as BOLD service benefits from taxonomy resolution
        """
        return True

    @staticmethod
    def requires_blastn():
        """
        Return False as BOLD ID service doesn't require local BLAST.

        :return: False, as no local BLAST is needed
        """
        return False