from nbitk.Taxon import Taxon
from typing import List, Optional, Tuple
import yaml
import csv

columns = set()

levels = [
    'kingdom',
    'phylum',
    'class',
    'order',
    'family',
    'subfamily',
    'tribe',
    'genus',
    'species',
    'subspecies'
]

class DNAAnalysisResult:

    def __init__(self, process_id: str):
        self.process_id: str = process_id
        self.data: dict = {
            'process_id': process_id,
            'nuc_basecount': None,
            'nuc_full_basecount': None,
            'obs_taxon': [],
            'identification': None,
            'species': None,
            'stop_codons': [],
            'ambig_basecount': None,
            'ambig_full_basecount': None,
            'identification_rank': None,
            'error': None,
            'identification_method': 'BLAST',
        }
        columns.update(self.data.keys())
        columns.discard(None)
        self.data['ancillary'] = {}

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
    def level(self) -> Optional[str]:
        """
        Getter for the taxonomic level.
        :return: A string representing the taxonomic level
        """
        return self.data['identification_rank']

    @level.setter
    def level(self, level: str) -> None:
        """
        Setter for the taxonomic level.
        :param level: A string representing the taxonomic level
        :return:
        """
        if not isinstance(level, str) or level.lower() not in levels:
            raise ValueError(f"level must be a string from {levels}")
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
    def obs_taxon(self, taxa: List[Taxon]) -> None:
        """
        Setter for the observed taxon.
        :param taxa: A list of strings representing the observed taxon
        :return:
        """
        if not isinstance(taxa, list) or not all(isinstance(item, Taxon) for item in taxa):
            raise ValueError("obs_taxon must be a list of Taxon objects")
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
        return self.seq_length >= 500 if self.seq_length is not None else False

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
        return len(self.stop_codons) == 0

    def check_ambiguities(self) -> bool:
        """
        Check if the sequence contains ambiguities, i.e. if the number of ambiguities is zero.
        :return: A boolean indicating whether the sequence contains ambiguities
        """
        return self.ambiguities == 0

    def check_seq_quality(self) -> bool:
        """
        Check if the sequence passes the quality checks for ambiguities and early stop codons
        :return: A boolean indicating whether the sequence passes all checks
        """
        return self.check_pseudogene() and self.check_ambiguities()

    def passes_all_checks(self) -> bool:
        """
        Check if the sequence passes all quality checks.
        :return: A boolean indicating whether the sequence passes all checks
        """
        return self.check_length() and self.check_taxonomy() and self.check_seq_quality()

    def calculate_ranks(self, verbosity: int = 2) -> Tuple[int, int, str]:
        """
        Calculate barcode_rank and full_rank, and generate messages based on verbosity.

        :param verbosity: 1=errors only, 2=errors+warnings, 3=errors+warnings+info
        :return: Tuple of (barcode_rank, full_rank, messages)
        """
        barcode_rank = 8
        full_rank = 6
        messages = []

        # Calculate barcode_rank
        if self.seq_length is not None and self.ambiguities is not None:
            if self.seq_length >= 650 and self.ambiguities == 0:
                barcode_rank = 1
                if verbosity >= 3:
                    messages.append("\U0001F947\U0001F31F BIN compliant, perfect")
            elif self.seq_length >= 500 and self.ambiguities == 0:
                barcode_rank = 2
                if verbosity >= 3:
                    messages.append("\U0001F947 BIN compliant")
            elif self.seq_length >= 650 and 1 <= self.ambiguities <= 6:
                barcode_rank = 3
                if verbosity >= 3:
                    messages.append("\U0001F948 BIN compliant")
                if verbosity >= 2:
                    messages.append("\u2753 Marker may be chimeric")
            elif self.seq_length >= 500 and 1 <= self.ambiguities <= 6:
                barcode_rank = 4
                if verbosity >= 3:
                    messages.append("\U0001F949 BIN compliant")
                if verbosity >= 2:
                    messages.append("\u2753 Marker may be chimeric")
            elif 400 <= self.seq_length < 500 and self.ambiguities == 0:
                barcode_rank = 5
                if verbosity >= 3:
                    messages.append("\u26A0 Useful marker sequence")
                if verbosity >= 2:
                    messages.append("\u2757 Not BIN compliant")
            elif 300 <= self.seq_length < 400 and self.ambiguities == 0:
                barcode_rank = 6
                if verbosity >= 3:
                    messages.append("\u26A0 Useful marker sequence")
                if verbosity >= 2:
                    messages.append("\u203C Not BIN compliant")
            elif (self.seq_length < 300 and self.ambiguities == 0) or (
                    self.seq_length < 500 and 1 <= self.ambiguities <= 6):
                barcode_rank = 7
                if verbosity >= 3:
                    messages.append("\u26A0 Useful marker sequence")
                if verbosity >= 2:
                    messages.append("\u203C Not BIN compliant")

        if barcode_rank == 8 and verbosity >= 1:
            messages.append("\u26D4 Unacceptable marker sequence")

        # Calculate full_rank
        if self.full_length is not None and self.full_ambiguities is not None:
            if self.full_length >= 1500 and self.full_ambiguities == 0:
                full_rank = 1
                if verbosity >= 3:
                    messages.append("\U0001F947 Excellent full length, no ambiguities")
            elif self.full_length >= 1000 and self.full_ambiguities == 0:
                full_rank = 2
                if verbosity >= 3:
                    messages.append("\U0001F948 Good full length, no ambiguities")
            elif self.full_length >= 1500 and 1 <= self.full_ambiguities < 15:
                full_rank = 3
                if verbosity >= 3:
                    messages.append("\U0001F947 Excellent full length")
                if verbosity >= 2:
                    messages.append("\u2757 Some ambiguities in full sequence")
            elif self.full_length >= 1000 and 1 <= self.full_ambiguities < 15:
                full_rank = 4
                if verbosity >= 3:
                    messages.append("\U0001F948 Good full length")
                if verbosity >= 2:
                    messages.append("\u2757 Some ambiguities in full sequence")
            elif self.full_length < 1000 and 0 <= self.full_ambiguities < 15:
                full_rank = 5
                if verbosity >= 3:
                    messages.append("\U0001F949 Acceptable full length")
                if self.full_ambiguities > 0 and verbosity >= 2:
                    messages.append("\u2757 Some ambiguities in full sequence")

        if full_rank == 6 and verbosity >= 1:
            messages.append("\u26D4 Unacceptable full sequence")

        return barcode_rank, full_rank, "\n".join(messages)

    def get_values(self) -> list:
        """
        String representation of the result object.
        :return: A list of values representing the result object
        """
        values = []
        for key in self.result_fields(self.level):
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
            elif key in levels:
                exp_taxon_name = self.exp_taxon.name if self.exp_taxon else None
                values.append(exp_taxon_name)
            else:
                values.append(self.data[key])
        return values

    def __str__(self) -> str:
        """
        String representation of the result object.
        :return: A tab-separated string representing the result object
        """
        return '\t'.join(map(str, self.get_values()))

    @classmethod
    def result_fields(cls, level: str = 'family') -> List[str]:
        """
        Returns a tab-separated string containing the result fields.
        :return:
        """
        columns.update([level])
        return sorted(item for item in columns if item is not None and item != 'exp_taxon')


class DNAAnalysisResultSet:

    def __init__(self, results: List[DNAAnalysisResult], level: str = 'family'):
        self.results = results
        self.level = level

    def __str__(self) -> str:
        """
        String representation of the result set.
        :return: A tab-separated string representing the result set
        """
        header = '\t'.join(DNAAnalysisResult.result_fields(self.level))
        contents = '\n'.join([str(result) for result in self.results])
        return header + "\n" + contents

    def add_yaml_file(self, file: str):
        """
        Join the YAML file to the results.
        :param file: YAML file
        :return:
        """

        # Read the YAML file and join it with the results
        with open(file, 'r') as yamlfile:
            yaml_data = yaml.safe_load(yamlfile)
            for result in self.results:
                for key, value in yaml_data.items():
                    result.ancillary[key] = value

    def add_csv_file(self, file: str):
        """
        Join the CSV file to the results.
        :param file: CSV file
        :return:
        """

        # Read the CSV file and join it with the results
        process_id_to_result = {result.process_id: result for result in self.results}
        with open(file, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                process_id = row['Process ID']
                if process_id in process_id_to_result:
                    result = process_id_to_result[process_id]

                    # Add all fields from CSV to the ancillary dictionary
                    for key, value in row.items():
                        if key != 'Process ID':  # Avoid duplicating the process_id
                            result.ancillary[key] = value