import json
import time
import datetime
from typing import Set, Dict, List, Optional
from Bio.SeqRecord import SeqRecord
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from json.decoder import JSONDecodeError

from nbitk.Taxon import Taxon
from nbitk.config import Config
from barcode_validator.constants import TaxonomicRank
from .idservice import IDService


class BOLD(IDService):
    """
    BOLD identification service implementation.

    This service uses the BOLD Systems identification engine to identify
    taxonomic information from sequence records. This is a self-contained
    implementation that doesn't require the id_engine module.
    """

    def __init__(self, config: Config):
        super().__init__(config)

        # Get BOLD-specific configuration with defaults
        self.database = config.get('bold_database', 1)  # Default to public.bin-tax-derep
        self.operating_mode = config.get('bold_operating_mode', 1)  # Default to 94% similarity
        self.timeout = config.get('bold_timeout', 300)  # 5 minutes default timeout

        # Build URL and parameters from config
        self.base_url, self.params = self._build_url_params()

        self.logger.info(
            f"Initialized BOLD service with database {self.database}, operating mode {self.operating_mode}")

    def _build_url_params(self) -> tuple:
        """
        Build the base URL and parameters for BOLD requests.

        :return: A tuple containing the base URL and a dictionary of parameters
        """
        # Database mapping
        idx_to_database = {
            1: "public.bin-tax-derep",
            2: "species",
            3: "all.bin-tax-derep",
            4: "DS-CANREF22",
            5: "public.plants",
            6: "public.fungi",
            7: "all.animal-alt",
            8: "DS-IUCNPUB",
        }

        # Operating mode mapping
        idx_to_operating_mode = {
            1: {"mi": 0.94, "maxh": 25},
            2: {"mi": 0.9, "maxh": 50},
            3: {"mi": 0.75, "maxh": 100},
        }

        if self.database not in idx_to_database:
            raise ValueError(f"Invalid database: {self.database}. Must be 1-8.")

        if self.operating_mode not in idx_to_operating_mode:
            raise ValueError(f"Invalid operating mode: {self.operating_mode}. Must be 1-3.")

        params = {
            "db": idx_to_database[self.database],
            "mi": idx_to_operating_mode[self.operating_mode]["mi"],
            "mo": 100,
            "maxh": idx_to_operating_mode[self.operating_mode]["maxh"],
            "order": 3,
        }

        base_url = "https://id.boldsystems.org/submission?db={}&mi={}&mo={}&maxh={}&order={}".format(
            params["db"], params["mi"], params["mo"], params["maxh"], params["order"]
        )

        return base_url, params

    def _submit_sequence(self, record: SeqRecord) -> str:
        """
        Submit a sequence to BOLD and return the submission ID.

        :param record: A Bio.SeqRecord object containing the sequence to submit
        :return: The submission ID returned by BOLD
        """
        # Create session with retry strategy
        session = requests.Session()
        session.headers.update({
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/98.0.4758.82 Safari/537.36"
        })
        retry_strategy = Retry(total=10, backoff_factor=1)
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("https://", adapter)

        # Format sequence data for submission
        data = f">{record.id}\n{record.seq}\n"
        files = {"file": ("submitted.fas", data, "text/plain")}

        try:
            while True:
                try:
                    # Submit the POST request
                    self.logger.debug(f"Submitting sequence {record.id} to BOLD")
                    response = session.post(self.base_url, params=self.params, files=files)
                    response.raise_for_status()
                    result = json.loads(response.text)
                    break
                except (JSONDecodeError, requests.RequestException) as e:
                    self.logger.warning(f"Request failed: {e}, retrying in 60 seconds")
                    time.sleep(60)

            sub_id = result['sub_id']
            self.logger.debug(f"Received submission ID: {sub_id}")
            return sub_id
        finally:
            session.close()

    def _wait_for_and_get_results(self, sub_id: str, record_id: str) -> List[Dict]:
        """
        Wait for BOLD processing to complete and return parsed results.

        :param sub_id: The submission ID returned by BOLD
        :param record_id: The ID of the sequence record to match results against
        :return: A list of dictionaries containing the parsed results
        """
        results_url = f"https://id.boldsystems.org/submission/results/{sub_id}"

        # Create session for polling results
        session = requests.Session()
        retry_strategy = Retry(total=5, backoff_factor=1)
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("https://", adapter)

        start_time = time.time()

        try:
            while time.time() - start_time < self.timeout:
                try:
                    self.logger.debug(f"Checking BOLD results for submission {sub_id}")
                    response = session.get(results_url)
                    response.raise_for_status()

                    # Try to parse JSON response
                    data = response.json()

                    # Check if processing is complete
                    # The response should be a dict with 'results', 'seqid', and 'sequence' keys
                    if (data and isinstance(data, dict) and
                            'results' in data and 'seqid' in data):

                        # Check if this is our sequence
                        if data.get("seqid") == record_id:
                            return self._parse_bold_result(data)
                        else:
                            self.logger.warning(
                                f"Results available but sequence ID mismatch: expected {record_id}, got {data.get('seqid')}")
                            return []

                    # Results not ready yet, wait and retry
                    self.logger.debug(f"Results not ready yet for submission {sub_id}, waiting...")
                    time.sleep(10)

                except requests.RequestException as e:
                    self.logger.debug(f"Request error while checking results: {e}, retrying...")
                    time.sleep(10)
                except json.JSONDecodeError:
                    # Results might not be ready if we can't parse JSON
                    self.logger.debug("Results not ready (invalid JSON), waiting...")
                    time.sleep(10)

            raise TimeoutError(f"BOLD processing timed out after {self.timeout} seconds")

        finally:
            session.close()

    def _parse_bold_result(self, result_data: Dict) -> List[Dict]:
        """
        Parse a BOLD result response into our standard format.

        :param result_data: The raw result data from BOLD
        :return: A list of dictionaries containing parsed results
        """
        results = []
        bold_results = result_data.get("results", {})

        if bold_results:
            for key, result_info in bold_results.items():
                # Parse the key format: process_id|primer|bin_uri|taxid|etc
                key_parts = key.split("|")
                process_id = key_parts[0] if len(key_parts) > 0 else ""
                primer = key_parts[1] if len(key_parts) > 1 else ""
                bin_uri = key_parts[2] if len(key_parts) > 2 else ""
                taxid = key_parts[3] if len(key_parts) > 3 else ""

                # Extract alignment metrics
                pident = result_info.get("pident", 0.0)
                bitscore = result_info.get("bitscore", 0.0)
                evalue = result_info.get("evalue", 1.0)

                # Extract taxonomy information
                taxonomy = result_info.get("taxonomy", {})
                result_dict = {
                    "phylum": taxonomy.get("phylum"),
                    "class": taxonomy.get("class"),
                    "order": taxonomy.get("order"),
                    "family": taxonomy.get("family"),
                    "subfamily": taxonomy.get("subfamily"),  # Include subfamily if present
                    "genus": taxonomy.get("genus"),
                    "species": taxonomy.get("species"),
                    "pct_identity": pident,
                    "bitscore": bitscore,
                    "evalue": evalue,
                    "process_id": process_id,
                    "primer": primer,
                    "bin_uri": bin_uri,
                    "taxid": taxid,
                    "taxid_count": taxonomy.get("taxid_count", "")
                }
                results.append(result_dict)
        else:
            # No matches found
            self.logger.debug(f"No matches found for sequence {result_data.get('seqid', 'unknown')}")

        return results

    def _build_taxonomic_trees(self, results: List[Dict]) -> List[Taxon]:
        """
        Build taxonomic trees from BOLD results, merging taxa with the same name.

        :param results: A list of dictionaries containing BOLD results
        :return: A list of root Taxon objects representing the taxonomic trees
        """
        # Dictionary to store unique taxa by (name, rank) pairs
        taxa_registry = {}

        # Define taxonomic levels in order from highest to lowest
        levels = [
            ("phylum", TaxonomicRank.PHYLUM),
            ("class", TaxonomicRank.CLASS),
            ("order", TaxonomicRank.ORDER),
            ("family", TaxonomicRank.FAMILY),
            ("subfamily", TaxonomicRank.FAMILY),  # Treat subfamily as family level for now
            ("genus", TaxonomicRank.GENUS),
            ("species", TaxonomicRank.SPECIES)
        ]

        # Track root taxa (highest level taxa with no parents)
        root_taxa = set()

        for result in results:
            previous_taxon = None

            # Build taxonomy chain from highest to lowest level
            for level_name, rank in levels:
                taxon_name = result.get(level_name)
                if not taxon_name or not taxon_name.strip():
                    continue

                taxon_name = taxon_name.strip()
                taxon_key = (taxon_name, rank.value)

                # Get or create taxon
                if taxon_key not in taxa_registry:
                    confidence = result.get("pct_identity", 0.0) / 100.0
                    taxon = Taxon(
                        name=taxon_name,
                        taxonomic_rank=rank.value,
                        confidence=confidence
                    )
                    taxa_registry[taxon_key] = taxon

                    # If this is the first (highest) level taxon, it's a root
                    if previous_taxon is None:
                        root_taxa.add(taxon)
                else:
                    taxon = taxa_registry[taxon_key]
                    # Update confidence if this result has higher confidence
                    current_confidence = result.get("pct_identity", 0.0) / 100.0
                    if current_confidence > taxon.confidence:
                        taxon.confidence = current_confidence

                # Link to parent if exists
                if previous_taxon is not None:
                    # Check if this taxon is already a child of the previous taxon
                    if taxon not in previous_taxon.clades:
                        previous_taxon.clades.append(taxon)

                previous_taxon = taxon

        return list(root_taxa)

    def _extract_taxa_at_level(self, trees: List[Taxon], level: TaxonomicRank) -> Set[Taxon]:
        """
        Extract all taxa at the specified taxonomic level from trees.

        :param trees: A list of root Taxon objects representing taxonomic trees
        :param level: The taxonomic rank at which to extract taxa
        :return: A set of Taxon objects at the specified level
        """
        taxa_at_level = set()

        def traverse_tree(taxon: Taxon, target_level: str):
            if taxon.taxonomic_rank == target_level:
                taxa_at_level.add(taxon)

            # Continue traversing children
            for child in taxon.clades:
                traverse_tree(child, target_level)

        for tree in trees:
            traverse_tree(tree, level.value)

        return taxa_at_level

    def identify_record(self, record: SeqRecord, level: TaxonomicRank = TaxonomicRank.FAMILY, extent: Taxon = None) -> \
    Set[Taxon]:
        """
        Identify the taxonomic classification of a sequence record using BOLD.

        :param record: A Bio.SeqRecord object containing the sequence to identify
        :param level: The taxonomic rank at which to return results (default: 'family')
        :param extent: Ignored for BOLD service
        :return: A set of Taxon objects representing the possible taxonomic classifications
        """
        try:
            self.logger.info(f"Identifying sequence {record.id} using BOLD")

            # Submit sequence to BOLD and get submission ID
            sub_id = self._submit_sequence(record)

            # Wait for processing and get results
            results = self._wait_for_and_get_results(sub_id, record.id)

            if not results:
                self.logger.warning(f"No BOLD results found for sequence {record.id}")
                return set()

            # Build taxonomic trees
            trees = self._build_taxonomic_trees(results)

            # Extract taxa at requested level
            taxa_at_level = self._extract_taxa_at_level(trees, level)

            self.logger.info(f"Found {len(taxa_at_level)} taxa at {level.value} level for sequence {record.id}")
            return taxa_at_level

        except Exception as e:
            self.logger.error(f"Error identifying sequence {record.id}: {str(e)}")
            return set()

    @staticmethod
    def requires_resolver():
        """BOLD service does not require a taxonomy resolver."""
        return False

    @staticmethod
    def requires_blastn():
        """BOLD service does not require BLAST."""
        return False