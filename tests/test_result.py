import pytest
from nbitk.Taxon import Taxon
from nbitk.config import Config

from barcode_validator.constants import TaxonomicRank
from barcode_validator.dna_analysis_result import DNAAnalysisResult

def test_init():

    # New object must have a seq ID and a data set ID and will instantiate a config object in the constructor
    res = DNAAnalysisResult('seq1', 'set1')
    assert res.dataset == 'set1'
    assert res.sequence_id == 'seq1'
    assert isinstance(res.config, Config)

def test_taxonomy():
    res = DNAAnalysisResult('seq1', 'set1')

    # set 'species', i.e. the expected identification at the lowest known level
    assert res.species is None
    res.species = Taxon('Homo sapiens', taxonomic_rank = TaxonomicRank.SPECIES.value)

    # set the expected taxon, i.e. the expected identification at the validation level
    assert res.exp_taxon is None
    res.exp_taxon = Taxon('Hominidae', taxonomic_rank= TaxonomicRank.FAMILY.value)

    # set the observed taxon list, i.e. whatever came out of the ID service
    assert len(res.obs_taxon) == 0
    res.obs_taxon = [Taxon('Hominidae', taxonomic_rank= TaxonomicRank.FAMILY.value)]

    # check if valid
    assert res.check_taxonomy()

def test_length():
    config = Config()
    config.config_data = { 'seq_length': 500 }
    config.initialized = True
    res = DNAAnalysisResult('seq1', 'set1', config)

    # test sequence length
    assert res.seq_length is None
    res.seq_length = 500
    assert res.seq_length == 500
    assert res.config.get('seq_length') == 500
    assert res.check_length()

    # test full length, which may be applicable for sequences with defined marker regions
    assert res.full_length is None
    res.full_length = 501
    assert res.full_length == 501

def test_ambiguities():
    config = Config()
    config.config_data = { 'ambiguities': 6 }
    config.initialized = True
    res = DNAAnalysisResult('seq1', 'set1', config)

    # test ambiguities
    assert res.ambiguities is None
    res.ambiguities = 6
    assert res.ambiguities == 6
    assert res.config.get('ambiguities') == 6
    assert res.check_ambiguities()

    # test full ambiguities, which may be applicable for sequences with defined marker regions
    assert res.full_ambiguities is None
    res.full_ambiguities = 7
    assert res.full_ambiguities == 7

def test_stop_codons():
    config = Config()
    config.config_data = { 'stop_codons': [] }
    config.initialized = True
    res = DNAAnalysisResult('seq1', 'set1', config)

    # test stop codons
    assert len(res.stop_codons) == 0
    res.stop_codons = [ 323 ]
    assert res.stop_codons == [ 323 ]
    assert res.config.get('stop_codons') == []
    assert not res.check_pseudogene()

def test_all():

    # Create a result object that should pass all checks
    config = Config()
    config.config_data = {
        'stop_codons': [],
        'ambiguities': 6,
        'seq_length': 500
    }
    config.initialized = True
    res = DNAAnalysisResult('seq1', 'set1', config)
    res.seq_length = 501
    res.ambiguities = 5
    res.stop_codons = []
    res.species = Taxon('Homo sapiens', taxonomic_rank = TaxonomicRank.SPECIES.value)
    res.exp_taxon = Taxon('Hominidae', taxonomic_rank=TaxonomicRank.FAMILY.value)
    res.obs_taxon = [ Taxon('Hominidae', taxonomic_rank=TaxonomicRank.FAMILY.value) ]
    assert res.passes_all_checks()

if __name__ == '__main__':
    pytest.main()