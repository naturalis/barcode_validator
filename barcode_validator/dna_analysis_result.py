from nbitk.Taxon import Taxon
from nbitk.config import Config
from barcode_validator.constants import TaxonomicRank, Marker, ValidationMode
from barcode_validator.criteria import MarkerCriteria
from typing import List, Optional, Tuple
import yaml
import csv

"""
The DNAAnalysisResult and DNAAnalysisResultSet classes are used to store and manipulate the results of barcode 
validation analyses. The DNAAnalysisResult class represents an analysis result for a single sequence, while the 
DNAAnalysisResultSet class represents a set of results, typically for a multifasta file. The aim of this design
is to be able to represent the results in various formats, with tabular data currently being the primary focus.
(The implementation could be expanded to serialize the results in other formats, such as RO-crate or JSON-LD.)

Under basic circumstances, the output has the following columns:
- sequence_id: An identifier for the sequence that is unique within the dataset (first word of the FASTA header)
- ambig_basecount: The number of ambiguous bases in the sequence within the barcode region
- ambig_full_basecount: The number of ambiguous bases in the full sequence
- dataset: The dataset name (e.g. the multifasta file name)
- error: An error message, if any
- identification: The expected taxon identification at level `identification_rank`
- identification_method: The method used for taxon identification (e.g. BLAST)
- identification_rank: The taxonomic level at which the identification was made (e.g. family)
- nuc_basecount: The number of nucleotides in the sequence within the barcode region
- nuc_full_basecount: The number of nucleotides in the full sequence
- obs_taxon: The distinct taxa observed via `identification_method` at level `identification_rank` 
- species: The expected species name, if provided
- stop_codons: The number of stop codons in the sequence

Furthermore, a validation procedure may be accompanied by ancillary data, which can be added to the result object.
There are two distinct scenarios where ancillary data may be added:
1. A CSV file is provided with additional analytics for each sequence, e.g. as produced by a skimming assembly pipeline.
   In this case, the CSV is joined with the results (by way of the sequence_id), and any additional columns from the
   CSV are added to the result object as ancillary data.
2. A YAML file is provided with additional metadata at the level of the dataset, e.g. the settings that were used for
   the skimming pipeline that produced the sequences. In this case, the YAML file is read and the data is added to the
   result object as ancillary data (i.e. every result object in the set will have the same ancillary data).
In order to ensure that the columns are consistent across all result objects in the set, the columns are tracked in a
global set called `columns`. This set is updated whenever a new column is added to any result object.
"""

# Module-level columns set
columns = set()

# Initial columns that should always be present
base_columns = {
    'sequence_id',
    'ambig_basecount',
    'ambig_full_basecount',
    'dataset',
    'error',
    'identification',
    'identification_method',
    'identification_rank',
    'nuc_basecount',
    'nuc_full_basecount',
    'obs_taxon',
    'species',
    'stop_codons',
    'sequence',
    'group_id'
}

def reset_columns():
    """Reset the global columns set to initial state"""
    global columns # noqa: F824
    columns.clear()  # Using clear() instead of reassignment

def initialize_columns():
    """Initialize the columns set with base columns"""
    reset_columns()
    global columns # noqa: F824
    columns.update(base_columns)

# Initialize columns at module import
initialize_columns()

class DNAAnalysisResult:

    def __init__(self, sequence_id: str, dataset: str = None, config: Config = None, group_id: str = None, criteria: MarkerCriteria = None):
        """
        Initialize a DNAAnalysisResult object.
        :param sequence_id: The sequence identifier
        :param dataset: The dataset name (e.g. the multifasta file name)
        """

        # Make a default config object for COI-5P
        if config is None:
            config = Config()
            config.config_data = {
                'seq_length': 500,
                'stop_codons': [],
                'ambiguities': 0
            }
            config.initialized = True
        self.config = config
        self.sequence_id: str = sequence_id
        self.criteria: MarkerCriteria = criteria
        self.data: dict = {'sequence_id': sequence_id, 'ambig_basecount': None,
                           'ambig_full_basecount': None, 'dataset': dataset, 'error': None,
                           'identification': None, 'identification_method': 'BLAST',
                           'identification_rank': None, 'nuc_basecount': None,
                           'nuc_full_basecount': None, 'obs_taxon': [],
                           'species': None, 'stop_codons': [], 'sequence': None,
                            'group_id': group_id, 'ancillary': {}}

    @property
    def marker_criteria(self) -> Optional[MarkerCriteria]:
        """
        Getter for the marker criteria.
        :return: A MarkerCriteria object
        """
        return self.criteria

    @marker_criteria.setter
    def marker_criteria(self, criteria: MarkerCriteria) -> None:
        """
        Setter for the marker criteria.
        :param criteria: A MarkerCriteria object
        :return:
        """
        if not isinstance(criteria, MarkerCriteria):
            raise ValueError("criteria must be a MarkerCriteria object")
        self.criteria = criteria

    @property
    def group_id(self) -> Optional[str]:
        """
        Getter for the group identifier.
        :return: A string representing the group identifier
        """
        return self.data['group_id']

    @group_id.setter
    def group_id(self, group_id: str) -> None:
        """
        Setter for the group identifier.
        :param group_id: A string representing the group identifier
        :return:
        """
        self.data['group_id'] = group_id

    @property
    def ancillary(self) -> dict:
        """
        Getter for the ancillary data.
        :return: A dictionary representing the ancillary data
        """
        return self.data['ancillary']

    @ancillary.setter
    def ancillary(self, data: dict) -> None:
        """
        Setter for the ancillary data.
        :param data: A dictionary representing the ancillary data
        :return:
        """
        if not isinstance(data, dict):
            raise ValueError("Ancillary data must be a dictionary")
        columns.update(data.keys())
        self.data['ancillary'].update(data)

    def add_ancillary(self, key: str, value: str) -> None:
        """
        Add an ancillary data item.
        :param key: A string representing the key
        :param value: A string representing the value
        :return:
        """
        columns.update([key])
        self.data['ancillary'][key] = value

    @property
    def error(self) -> Optional[str]:
        """
        Getter for the error message.
        :return: A string representing the error message
        """
        return self.data['error']

    @error.setter
    def error(self, error: str) -> None:
        """
        Setter for the error message.
        :param error: A string representing the error message
        :return:
        """
        self.data['error'] = error

    @property
    def dataset(self) -> Optional[str]:
        """
        Getter for the dataset.
        :return: A string representing the dataset (e.g. FASTA file name)
        """
        return self.data['dataset']

    @dataset.setter
    def dataset(self, dataset: str) -> None:
        """
        Setter for the dataset.
        :param dataset: A string representing the dataset (e.g. FASTA file name)
        :return:
        """
        self.data['dataset'] = dataset

    @property
    def level(self) -> Optional[TaxonomicRank]:
        """
        Getter for the taxonomic level.
        :return: A string representing the taxonomic level
        """
        return self.data['identification_rank']

    @level.setter
    def level(self, level: TaxonomicRank) -> None:
        """
        Setter for the taxonomic level.
        :param level: A TaxanomicRank enum
        :return:
        """
        self.data['identification_rank'] = level

    @property
    def seq_length(self) -> Optional[int]:
        """
        Getter for the sequence length within the marker region.
        :return: an integer representing the sequence length
        """
        return self.data['nuc_basecount']

    @seq_length.setter
    def seq_length(self, value: int) -> None:
        """
        Setter for the sequence length within the marker region.
        :param value: an integer representing the sequence length
        :return:
        """
        if not isinstance(value, int) or value < 0:
            raise ValueError("seq_length must be a positive integer")
        self.data['nuc_basecount'] = value

    @property
    def full_length(self) -> Optional[int]:
        """
        Getter for the full sequence length.
        :return: an integer representing the sequence length
        """
        return self.data['nuc_full_basecount']

    @full_length.setter
    def full_length(self, value: int) -> None:
        """
        Setter for the full sequence length.
        :param value: an integer representing the sequence length
        :return:
        """
        if not isinstance(value, int) or value < 0:
            raise ValueError("full_length must be a positive integer")
        self.data['nuc_full_basecount'] = value

    @property
    def obs_taxon(self) -> List[Taxon]:
        """
        Getter for the observed taxon.
        :return: A list of strings representing the observed taxon
        """
        return self.data['obs_taxon']

    @obs_taxon.setter
    def obs_taxon(self, taxa) -> None:
        """
        Setter for the observed taxon.
        :param taxa: A list or set of Taxon objects representing the observed taxon
        :return:
        """
        # Handle both list and set inputs
        if isinstance(taxa, set):
            taxa = list(taxa)

        if not isinstance(taxa, list) or not all(isinstance(item, Taxon) for item in taxa):
            raise ValueError("obs_taxon must be a list or set of Taxon objects")
        self.data['obs_taxon'] = taxa

    def add_obs_taxon(self, taxon: Taxon) -> None:
        """
        Add an observed taxon to the list.
        :param taxon: A string representing the observed taxon
        :return:
        """
        if not isinstance(taxon, Taxon):
            raise ValueError("Taxon must be a Taxon object")
        if taxon not in self.data['obs_taxon']:
            self.data['obs_taxon'].append(taxon)

    @property
    def exp_taxon(self) -> Optional[Taxon]:
        """
        Getter for the expected taxon.
        :return: A Taxon object representing the expected taxon
        """
        return self.data['identification']

    @exp_taxon.setter
    def exp_taxon(self, taxon: Taxon) -> None:
        """
        Setter for the expected taxon.
        :param taxon: A Taxon object representing the expected taxon
        :return:
        """
        if not isinstance(taxon, Taxon):
            raise ValueError("exp_taxon must be a Taxon object")
        self.data['identification'] = taxon

    @property
    def species(self) -> Optional[Taxon]:
        """
        Getter for the species name.
        :return: A Taxon object representing the species name
        """
        return self.data['species']

    @species.setter
    def species(self, species: Taxon) -> None:
        """
        Setter for the species name.
        :param species: A Taxon object representing the species name
        :return:
        """
        if not isinstance(species, Taxon):
            raise ValueError("species must be a Taxon object")
        self.data['species'] = species

    @property
    def stop_codons(self) -> List[int]:
        """
        Getter for the stop codons.
        :return: A list of integers representing the stop codon positions
        """
        return self.data['stop_codons']

    @stop_codons.setter
    def stop_codons(self, codon_positions: List[int]) -> None:
        """
        Setter for the stop codons.
        :param codon_positions: A list of integers representing the stop codon positions
        :return:
        """
        if not isinstance(codon_positions, list) or not all(isinstance(x, int) and x >= 0 for x in codon_positions):
            raise ValueError("stop_codons must be a list of non-negative integers")
        self.data['stop_codons'] = codon_positions

    @property
    def ambiguities(self) -> Optional[int]:
        """
        Getter for the number of ambiguities within the marker region.
        :return: An integer representing the number of ambiguities
        """
        return self.data['ambig_basecount']

    @ambiguities.setter
    def ambiguities(self, n_ambiguities: int) -> None:
        """
        Setter for the number of ambiguities within the marker region.
        :param n_ambiguities: An integer representing the number of ambiguities
        :return:
        """
        if not isinstance(n_ambiguities, int) or n_ambiguities < 0:
            raise ValueError("ambiguities must be a non-negative integer")
        self.data['ambig_basecount'] = n_ambiguities

    @property
    def full_ambiguities(self) -> Optional[int]:
        """
        Getter for the total number of ambiguities.
        :return: An integer representing the number of ambiguities
        """
        return self.data['ambig_full_basecount']

    @full_ambiguities.setter
    def full_ambiguities(self, n_ambiguities: int) -> None:
        """
        Setter for the total number of ambiguities.
        :param n_ambiguities: An integer representing the number of ambiguities
        :return:
        """
        if not isinstance(n_ambiguities, int) or n_ambiguities < 0:
            raise ValueError("ambiguities must be a non-negative integer")
        self.data['ambig_full_basecount'] = n_ambiguities

    def add_stop_codon(self, position: int) -> None:
        """
        Add a stop codon position to the list.
        :param position: An integer representing the stop codon position
        :return:
        """
        if not isinstance(position, int) or position < 0:
            raise ValueError("Stop codon position must be a non-negative integer")
        self.data['stop_codons'].append(position)

    def check_length(self) -> bool:
        """
        Check if the sequence length meets the minimum requirement.
        :return: A boolean indicating whether the sequence length is valid
        """
        if self.marker_criteria is None:
            raise ValueError("Marker criteria is not set")
        return self.seq_length >= int(self.marker_criteria.min_length) if self.seq_length is not None else False

    def check_taxonomy(self) -> bool:
        """
        Check if expected taxon is in the observed taxon list.
        :return: A boolean indicating whether the taxonomy check passed
        """
        return self.exp_taxon in [taxon for taxon in self.obs_taxon] if self.obs_taxon and self.exp_taxon else False

    def check_pseudogene(self) -> bool:
        """
        Check if the sequence contains stop codons, i.e. if the list of stop codon locations is empty.
        :return: A boolean indicating whether the sequence is a pseudogene
        """
        if self.marker_criteria is None:
            raise ValueError("Marker criteria is not set")
        return len(self.stop_codons) <= self.marker_criteria.max_stop_codons

    def check_ambiguities(self) -> bool:
        """
        Check if the sequence contains ambiguities, i.e. if the number of ambiguities is zero.
        :return: A boolean indicating whether the sequence contains ambiguities
        """
        if self.marker_criteria is None:
            raise ValueError("Marker criteria is not set")
        return self.ambiguities <= self.marker_criteria.max_ambiguities

    def check_seq_quality(self) -> bool:
        """
        Check if the sequence passes the quality checks for ambiguities and early stop codons
        :return: A boolean indicating whether the sequence passes all checks
        """
        return self.check_pseudogene() and self.check_ambiguities() and self.check_length()

    def passes_all_checks(self) -> bool:
        """
        Check if the sequence passes all quality checks.
        :return: A boolean indicating whether the sequence passes all checks
        """
        return self.check_taxonomy() and self.check_seq_quality()

    def get_values(self) -> list:
        """
        Get the values of the result object.
        :return: A list of values representing the result object
        """
        values = []
        for key in self.result_fields():
            if key == 'identification':
                exp_taxon_name = self.exp_taxon.name if self.exp_taxon else None
                values.append(exp_taxon_name)
            elif key == 'species':
                species_name = self.species.name if self.species else None
                values.append(species_name)
            elif key == 'obs_taxon':
                obs = [taxon.name for taxon in self.obs_taxon]
                values.append(",".join(obs))
            elif key == 'stop_codons':
                values.append(len(self.stop_codons))
            elif key in self.data['ancillary']:
                anc = self.data.get('ancillary')[key]
                values.append(str(anc))
            else:
                if key in self.data:
                    values.append(self.data[key])
                else:
                    values.append(None)
        return values

    def __str__(self) -> str:
        """
        String representation of the result object.
        :return: A tab-separated string representing the result object
        """
        return '\t'.join(map(str, self.get_values()))

    @classmethod
    def result_fields(cls) -> List[str]:
        """
        Returns a tab-separated string containing the result fields.
        :return:
        """
        return sorted(item for item in columns if item is not None)


class DNAAnalysisResultSet:
    def __init__(self, results: List[DNAAnalysisResult], config: Config = None):
        # Reset and initialize columns for this new result set
        initialize_columns()

        # Set the configuration object
        if config is None:
            config = Config()
            config.config_data = { 'group_id_separator': '_' }
            config.initialized = True
        self.config = config

        self.results = results

        # Update columns based on all results in the set
        for result in results:
            # Add any ancillary columns from existing results
            if result.data['ancillary']:
                columns.update(result.data['ancillary'].keys())


    def __str__(self) -> str:
        """
        String representation of the result set.
        :return: A tab-separated string representing the result set
        """
        header = '\t'.join(DNAAnalysisResult.result_fields())
        contents = '\n'.join([str(result) for result in self.results])
        return header + "\n" + contents

    def to_string(self, output_format: str = 'tsv') -> str:
        """
        Convert the result set to a string.
        :param output_format: The output format, tsv or fasta (default: tsv)
        :return: A string representing the result set
        """
        if output_format.lower() == 'tsv':
            return str(self)
        elif output_format.lower() == 'fasta':
            return "\n".join([f">{result.sequence_id}\n{result.data['ancillary']['nuc']}" for result in self.results])
        else:
            raise ValueError(f"Output format '{output_format}' not supported")

    def add_yaml_file(self, file: str):
        """
        Join the YAML file to the results.
        :param file: YAML file
        :return:
        """
        # Read the YAML file and join it with the results
        with open(file, 'r') as yamlfile:
            yaml_data = yaml.safe_load(yamlfile)
            # Update columns first with all possible keys from YAML
            columns.update(yaml_data.keys())
            # Then add data to each result
            for result in self.results:
                for key, value in yaml_data.items():
                    result.add_ancillary(key, value)

    def add_csv_file(self, file: str):
        """
        Join the CSV file to the results.
        :param file: CSV file
        :return:
        """
        # The CSV file has process IDs, not sequence IDs, so we need to map them to the results
        process_id_to_result = {}
        for result in self.results:
            seqid = result.sequence_id
            process_id = seqid.split(self.config.get('group_id_separator'))[0]
            process_id_to_result[process_id] = result

        with open(file, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            # Update columns first with all possible fields from CSV
            columns.update(field for field in reader.fieldnames if field != 'ID')

            # Reset file pointer to start
            csvfile.seek(0)
            reader = csv.DictReader(csvfile)

            # Then add data to each result
            for row in reader:
                process_id = row['ID']
                if process_id in process_id_to_result:
                    result = process_id_to_result[process_id]
                    for key, value in row.items():
                        if key != 'ID':  # Avoid duplicating the process_id
                            result.add_ancillary(key, value)

    def triage(self, mode: ValidationMode, aggregate: bool = False) -> 'DNAAnalysisResultSet':
        """
        Perform triage on the result set.
        :return: A new DNAAnalysisResultSet object containing the triaged results
        """
        if mode == ValidationMode.STRUCTURAL:
            triaged_results = [result for result in self.results if result.check_seq_quality()]
        elif mode == ValidationMode.TAXONOMIC:
            triaged_results = [result for result in self.results if result.check_taxonomy()]
        else:
            triaged_results = [result for result in self.results if result.passes_all_checks()]

        if aggregate:

            # Pick the longest sequence for each group_id
            aggregated = {}
            for item in triaged_results:
                group_id = item.group_id
                if group_id not in aggregated or item.seq_length > aggregated[group_id].seq_length:
                    aggregated[group_id] = item
            filtered_results = list(aggregated.values())

            # Specify that this is a BOLD submission
            for item in filtered_results:
                item.add_ancillary('BOLD_submission', item.group_id)

            return DNAAnalysisResultSet(filtered_results, self.config)
        else:
            return DNAAnalysisResultSet(triaged_results, self.config)