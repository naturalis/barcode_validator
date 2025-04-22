import pytest
from barcode_validator.criteria import (
    MarkerCriteria, COI5PCriteria, MATKCriteria, RBCLCriteria,
    ITSCriteria, ITS2Criteria, MarkerCriteriaFactory
)
from barcode_validator.constants import Marker
from nbitk.config import Config


class TestMarkerCriteria:
    """Test suite for the MarkerCriteria base class."""

    def test_default_initialization(self):
        """Test that default initialization sets expected values."""
        criteria = MarkerCriteria()
        assert criteria.min_length == 500
        assert criteria.max_ambiguities == 0
        assert criteria.max_stop_codons == 0
        assert criteria.marker_type is None

    def test_config_initialization(self):
        """Test initialization with configuration data."""
        config = Config()
        config.config_data = {
            'min_length': 600,
            'max_ambiguities': 10,
            'max_stop_codons': 2
        }
        config.initialized = True
        criteria = MarkerCriteria(config)
        assert criteria.min_length == 600
        assert criteria.max_ambiguities == 10
        assert criteria.max_stop_codons == 2

    def test_update_from_config(self):
        """Test updating from configuration data."""
        criteria = MarkerCriteria()
        config = Config()
        config.config_data = {
            'min_length': 600,
            'max_ambiguities': 10,
            'max_stop_codons': 2
        }
        config.initialized = True
        criteria.update_from_config(config)
        assert criteria.min_length == 600
        assert criteria.max_ambiguities == 10
        assert criteria.max_stop_codons == 2

    def test_to_dict(self):
        """Test conversion to dictionary."""
        criteria = MarkerCriteria()
        criteria.marker_type = Marker.COI_5P
        expected = {
            'min_length': 500,
            'max_ambiguities': 0,
            'max_stop_codons': 0,
            'marker_type': 'COI-5P'
        }
        assert criteria.to_dict() == expected

    def test_string_representation(self):
        """Test string representation."""
        criteria = MarkerCriteria()
        criteria.marker_type = Marker.COI_5P
        assert "COI-5P Criteria: min_length=500" in str(criteria)
        assert "max_ambiguities=0" in str(criteria)
        assert "max_stop_codons=0" in str(criteria)


class TestCOI5PCriteria:
    """Test suite for the COI5PCriteria class."""

    def test_default_initialization(self):
        """Test that default initialization sets expected values for COI-5P."""
        criteria = COI5PCriteria()
        assert criteria.min_length == 500
        assert criteria.max_ambiguities == 6
        assert criteria.max_stop_codons == 0
        assert criteria.marker_type == Marker.COI_5P

    def test_config_override(self):
        """Test that configuration data overrides defaults."""
        config = Config()
        config.config_data = {
            'min_length': 600,
            'max_ambiguities': 10,
            'max_stop_codons': 2
        }
        config.initialized = True
        criteria = COI5PCriteria(config)
        assert criteria.min_length == 600
        assert criteria.max_ambiguities == 10
        assert criteria.max_stop_codons == 2
        assert criteria.marker_type == Marker.COI_5P


class TestMATKCriteria:
    """Test suite for the MATKCriteria class."""

    def test_default_initialization(self):
        """Test that default initialization sets expected values for matK."""
        criteria = MATKCriteria()
        assert criteria.min_length == 700
        assert criteria.max_ambiguities == 20
        assert criteria.max_stop_codons == 0
        assert criteria.marker_type == Marker.MATK


class TestRBCLCriteria:
    """Test suite for the RBCLCriteria class."""

    def test_default_initialization(self):
        """Test that default initialization sets expected values for rbcL."""
        criteria = RBCLCriteria()
        assert criteria.min_length == 550
        assert criteria.max_ambiguities == 15
        assert criteria.max_stop_codons == 0
        assert criteria.marker_type == Marker.RBCL


class TestITSCriteria:
    """Test suite for the ITSCriteria class."""

    def test_default_initialization(self):
        """Test that default initialization sets expected values for ITS."""
        criteria = ITSCriteria()
        assert criteria.min_length == 500
        assert criteria.max_ambiguities == 30
        assert criteria.max_stop_codons is None  # Not relevant for non-coding
        assert criteria.marker_type == Marker.ITS


class TestITS2Criteria:
    """Test suite for the ITS2Criteria class."""

    def test_default_initialization(self):
        """Test that default initialization sets expected values for ITS2."""
        criteria = ITS2Criteria()
        assert criteria.min_length == 350
        assert criteria.max_ambiguities == 20
        assert criteria.max_stop_codons is None  # Not relevant for non-coding
        assert criteria.marker_type == Marker.ITS2


class TestMarkerCriteriaFactory:
    """Test suite for the MarkerCriteriaFactory class."""

    def test_get_criteria(self):
        """Test getting criteria by marker type."""
        criteria = MarkerCriteriaFactory.get_criteria(Marker.COI_5P)
        assert isinstance(criteria, COI5PCriteria)

        criteria = MarkerCriteriaFactory.get_criteria(Marker.MATK)
        assert isinstance(criteria, MATKCriteria)

        criteria = MarkerCriteriaFactory.get_criteria(Marker.RBCL)
        assert isinstance(criteria, RBCLCriteria)

        criteria = MarkerCriteriaFactory.get_criteria(Marker.ITS)
        assert isinstance(criteria, ITSCriteria)

        criteria = MarkerCriteriaFactory.get_criteria(Marker.ITS2)
        assert isinstance(criteria, ITS2Criteria)

    def test_from_string(self):
        """Test getting criteria by marker string."""
        criteria = MarkerCriteriaFactory.from_string('COI-5P')
        assert isinstance(criteria, COI5PCriteria)

        criteria = MarkerCriteriaFactory.from_string('matK')
        assert isinstance(criteria, MATKCriteria)

        criteria = MarkerCriteriaFactory.from_string('rbcL')
        assert isinstance(criteria, RBCLCriteria)

        criteria = MarkerCriteriaFactory.from_string('ITS')
        assert isinstance(criteria, ITSCriteria)

        criteria = MarkerCriteriaFactory.from_string('ITS2')
        assert isinstance(criteria, ITS2Criteria)

    def test_unknown_marker(self):
        """Test that unknown marker types raise ValueError."""
        with pytest.raises(ValueError):
            MarkerCriteriaFactory.from_string('unknown_marker')

    def test_get_available_markers(self):
        """Test getting available markers."""
        markers = MarkerCriteriaFactory.get_available_markers()
        assert 'COI-5P' in markers
        assert 'matK' in markers
        assert 'rbcL' in markers
        assert 'ITS' in markers
        assert 'ITS2' in markers
        assert len(markers) == 5  # Update if more markers are added