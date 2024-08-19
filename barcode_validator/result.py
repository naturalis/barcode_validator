from barcode_validator.config import Config
import logging


def result_fields():
    """
    Returns a tab-separated string containing the result fields.
    :return:
    """
    return "Process ID\tSequence Length\tObserved Taxon\tExpected Taxon\tSpecies\tStop Codons\tAmbiguities\tPasses " \
           "All Checks"


class DNAAnalysisResult:
    def __init__(self, process_id):
        self.process_id = process_id
        self._seq_length = None
        self._obs_taxon = []  # Changed to an empty list
        self._exp_taxon = None
        self._species = None
        self._stop_codons = []
        self._ambiguities = None

    @property
    def seq_length(self):
        """
        Getter for the sequence length.
        :return: an integer representing the sequence length
        """
        return self._seq_length

    @seq_length.setter
    def seq_length(self, value):
        """
        Setter for the sequence length.
        :param value: an integer representing the sequence length
        :return:
        """
        if not isinstance(value, int) or value <= 0:
            raise ValueError("seq_length must be a positive integer")
        self._seq_length = value

    @property
    def obs_taxon(self):
        """
        Getter for the observed taxon.
        TODO: Change this to a list of Taxon objects
        :return: A list of strings representing the observed taxon
        """
        return self._obs_taxon

    @obs_taxon.setter
    def obs_taxon(self, value):
        """
        Setter for the observed taxon.
        TODO: Change this to a list of Taxon objects
        :param value: A list of strings representing the observed taxon
        :return:
        """
        if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
            logging.error(value)
            raise ValueError("obs_taxon must be a list of strings")
        self._obs_taxon = value

    def add_obs_taxon(self, taxon):
        """
        Add an observed taxon to the list.
        TODO: Change this to a Taxon object
        :param taxon: A string representing the observed taxon
        :return:
        """
        if not isinstance(taxon, str):
            raise ValueError("Taxon must be a string")
        if taxon not in self._obs_taxon:
            self._obs_taxon.append(taxon)

    @property
    def exp_taxon(self):
        """
        Getter for the expected taxon.
        :return: A string representing the expected taxon
        """
        return self._exp_taxon

    @exp_taxon.setter
    def exp_taxon(self, value):
        """
        Setter for the expected taxon.
        TODO: Change this to a Taxon object
        :param value: A string representing the expected taxon
        :return:
        """
        if not isinstance(value, str):
            raise ValueError("exp_taxon must be a string")
        self._exp_taxon = value

    @property
    def species(self):
        """
        Getter for the species name.
        TODO: Change this to a Taxon object
        :return: A string representing the species name
        """
        return self._species

    @species.setter
    def species(self, value):
        """
        Setter for the species name.
        TODO: Change this to a Taxon object
        :param value: A string representing the species name
        :return:
        """
        if not isinstance(value, str):
            raise ValueError("species must be a string")
        self._species = value

    @property
    def stop_codons(self):
        """
        Getter for the stop codons.
        :return: A list of integers representing the stop codon positions
        """
        return self._stop_codons

    @stop_codons.setter
    def stop_codons(self, value):
        """
        Setter for the stop codons.
        :param value: A list of integers representing the stop codon positions
        :return:
        """
        if not isinstance(value, list) or not all(isinstance(x, int) and x >= 0 for x in value):
            raise ValueError("stop_codons must be a list of non-negative integers")
        self._stop_codons = value

    @property
    def ambiguities(self):
        """
        Getter for the number of ambiguities.
        :return: An integer representing the number of ambiguities
        """
        return self._ambiguities

    @ambiguities.setter
    def ambiguities(self, value):
        """
        Setter for the number of ambiguities.
        :param value: An integer representing the number of ambiguities
        :return:
        """
        if not isinstance(value, int) or value < 0:
            raise ValueError("ambiguities must be a non-negative integer")
        self._ambiguities = value

    def add_stop_codon(self, position):
        """
        Add a stop codon position to the list.
        :param position: An integer representing the stop codon position
        :return:
        """
        if not isinstance(position, int) or position < 0:
            raise ValueError("Stop codon position must be a non-negative integer")
        self._stop_codons.append(position)

    def check_length(self):
        """
        Check if the sequence length meets the minimum requirement.
        :return: A boolean indicating whether the sequence length is valid
        """
        config = Config()
        min_length = config.get('min_seq_length')
        return self.seq_length >= min_length if self.seq_length is not None else False

    def check_taxonomy(self):
        """
        Check if expected taxon is in the observed taxon list.
        TODO: Change this to compare Taxon objects
        :return: A boolean indicating whether the taxonomy check passed
        """
        return self.exp_taxon in [taxon for taxon in self.obs_taxon] if self.obs_taxon and self.exp_taxon else False

    def check_pseudogene(self):
        """
        Check if the sequence contains stop codons, i.e. if the list of stop codon locations is empty.
        :return: A boolean indicating whether the sequence is a pseudogene
        """
        return len(self.stop_codons) == 0

    def passes_all_checks(self):
        """
        Check if the sequence passes all quality checks.
        :return: A boolean indicating whether the sequence passes all checks
        """
        return self.check_length() and self.check_taxonomy() and self.check_pseudogene()

    def __str__(self):
        """
        String representation of the result object.
        :return: A tab-separated string containing the result fields
        """
        results = [self.process_id, self.seq_length, ', '.join(self.obs_taxon), self.exp_taxon, self.species,
                   self.stop_codons, self.ambiguities, self.passes_all_checks()]
        return '\t'.join(map(str, results))

