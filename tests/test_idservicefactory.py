import pytest
from barcode_validator.idservices.factory import IDServiceFactory
from barcode_validator.constants import RefDB
from barcode_validator.config.schema_config import SchemaConfig

@pytest.fixture
def config():
    """Fixture to provide a Config object for testing."""
    config = SchemaConfig()
    config.set('log_level', 'DEBUG')
    config.set('reflib_taxonomy_type', 'ncbi')
    return config

def test_create_idservice_ncbi(config):
    """Test creating an NCBI ID service."""
    db = RefDB(config.get('reflib_taxonomy_type'))
    id_service = IDServiceFactory.create_idservice(config, db)
    assert id_service is not None
    assert id_service.__class__.__name__ == "NCBI"

if __name__ == '__main__':
    pytest.main(['-v'])