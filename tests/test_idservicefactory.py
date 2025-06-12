import pytest
from barcode_validator.idservices.factory import IDServiceFactory
from barcode_validator.constants import RefDB
from nbitk.config import Config

@pytest.fixture
def config():
    """Fixture to provide a Config object for testing."""
    config = Config()
    config.config_data = { "reference_taxonomy": "ncbi", "log_level": "DEBUG" }
    config.initialized = True
    return config

def test_create_idservice_ncbi(config):
    """Test creating an NCBI ID service."""
    db = RefDB(config.get('reference_taxonomy'))
    id_service = IDServiceFactory.create_idservice(config, db)
    assert id_service is not None
    assert id_service.__class__.__name__ == "NCBI"