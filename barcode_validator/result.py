from config import Config


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
        return self._seq_length

    @seq_length.setter
    def seq_length(self, value):
        if not isinstance(value, int) or value <= 0:
            raise ValueError("seq_length must be a positive integer")
        self._seq_length = value

    @property
    def obs_taxon(self):
        return self._obs_taxon

    @obs_taxon.setter
    def obs_taxon(self, value):
        if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
            logging.error(value)
            raise ValueError("obs_taxon must be a list of strings")
        self._obs_taxon = value

    def add_obs_taxon(self, taxon):
        if not isinstance(taxon, str):
            raise ValueError("Taxon must be a string")
        if taxon not in self._obs_taxon:
            self._obs_taxon.append(taxon)

    @property
    def exp_taxon(self):
        return self._exp_taxon

    @exp_taxon.setter
    def exp_taxon(self, value):
        if not isinstance(value, str):
            raise ValueError("exp_taxon must be a string")
        self._exp_taxon = value

    @property
    def species(self):
        return self._species

    @species.setter
    def species(self, value):
        if not isinstance(value, str):
            raise ValueError("species must be a string")
        self._species = value

    @property
    def stop_codons(self):
        return self._stop_codons

    @stop_codons.setter
    def stop_codons(self, value):
        if not isinstance(value, list) or not all(isinstance(x, int) and x >= 0 for x in value):
            raise ValueError("stop_codons must be a list of non-negative integers")
        self._stop_codons = value

    @property
    def ambiguities(self):
        return self._ambiguities

    @ambiguities.setter
    def ambiguities(self, value):
        if not isinstance(value, int) or value < 0:
            raise ValueError("ambiguities must be a non-negative integer")
        self._ambiguities = value

    def add_stop_codon(self, position):
        if not isinstance(position, int) or position < 0:
            raise ValueError("Stop codon position must be a non-negative integer")
        self._stop_codons.append(position)

    def check_length(self):
        config = Config()
        min_length = config.get('min_seq_length')
        return self.seq_length >= min_length if self.seq_length is not None else False

    def check_taxonomy(self):
        return self.exp_taxon in [taxon for taxon in self.obs_taxon] if self.obs_taxon and self.exp_taxon else False

    def check_pseudogene(self):
        return len(self.stop_codons) == 0

    def passes_all_checks(self):
        return self.check_length() and self.check_taxonomy() and self.check_pseudogene()

    def __str__(self):
        results = [self.process_id, self.seq_length, ', '.join(self.obs_taxon), self.exp_taxon, self.species,
                   self.stop_codons, self.ambiguities, self.passes_all_checks()]
        return '\t'.join(map(str, results))