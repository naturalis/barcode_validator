from config import Config


class DNAAnalysisResult:
    def __init__(self, process_id):
        self.process_id = process_id
        self._seq_length = None
        self._obs_taxon = None
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
        if not isinstance(value, str):
            raise ValueError("obs_taxon must be a string")
        self._obs_taxon = value

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
        return self.obs_taxon.lower() == self.exp_taxon.lower() if self.obs_taxon and self.exp_taxon else False

    def check_pseudogene(self):
        return len(self.stop_codons) == 0

    def passes_all_checks(self):
        return self.check_length() and self.check_taxonomy() and self.check_pseudogene()

    def __str__(self):
        return f"DNAAnalysisResult for {self.process_id}:\n" \
               f" Sequence Length: {self.seq_length}\n" \
               f" Observed Taxon: {self.obs_taxon}\n" \
               f" Expected Taxon: {self.exp_taxon}\n" \
               f" Species: {self.species}\n" \
               f" Stop Codons: {self.stop_codons}\n" \
               f" Ambiguities: {self.ambiguities}\n" \
               f" Passes All Checks: {self.passes_all_checks()}"

