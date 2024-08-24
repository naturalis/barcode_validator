import pytest
import json
from barcode_validator.provenance import ROCrateGenerator


@pytest.fixture
def crate():
    return ROCrateGenerator('ro-crate-meta.json', 'common_meta.json')


def test_initialization(crate):
    assert crate.context == "https://w3id.org/ro/crate/1.1/context"
    assert len(crate.graph) == 2
    assert crate.graph[0]["@id"] == "ro-crate-meta.json"
    assert crate.graph[1]["@id"] == "./"
    assert crate.graph[1]["mainEntity"]["@id"] == "common_meta.json"


def test_add_common_metadata(crate):
    crate.add_common_metadata("common_meta.json", "Test Metadata", "Test Description")
    metadata = next(item for item in crate.graph if item["@id"] == "common_meta.json")
    assert metadata["name"] == "Test Metadata"
    assert metadata["description"] == "Test Description"
    assert metadata["encodingFormat"] == "application/json"


def test_add_file(crate):
    crate.add_file("test.fa", "Test File", "Test FASTA File", "application/x-fasta")
    file_metadata = next(item for item in crate.graph if item["@id"] == "test.fa")
    assert file_metadata["name"] == "Test File"
    assert file_metadata["description"] == "Test FASTA File"
    assert file_metadata["encodingFormat"] == "application/x-fasta"


def test_add_dna_analysis(crate):
    result = {
        "@type": "DNAAnalysisResult",
        "processId": "TEST-01",
        "sequenceLength": 100,
        "passesAllChecks": True
    }
    crate.add_dna_analysis("#test_analysis", "test.fa", result)
    analysis = next(item for item in crate.graph if item["@id"] == "#test_analysis")
    assert analysis["object"]["@id"] == "test.fa"
    assert analysis["result"] == result


def test_add_software(crate):
    crate.add_software("#test_software", "TestSoft", "1.0", "Test Software Description")
    software = next(item for item in crate.graph if item["@id"] == "#test_software")
    assert software["name"] == "TestSoft"
    assert software["version"] == "1.0"
    assert software["description"] == "Test Software Description"


def test_generate_json(crate):
    json_output = crate.generate_json()
    assert "@context" in json_output
    assert "@graph" in json_output
    assert len(json_output["@graph"]) == 2  # Initial two items


def test_save_to_file(crate, tmp_path):
    file_path = tmp_path / "test_crate.json"
    crate.save_to_file(str(file_path))
    assert file_path.exists()
    with open(file_path, 'r') as f:
        saved_data = json.load(f)
    assert "@context" in saved_data
    assert "@graph" in saved_data


def test_full_crate_generation():
    crate = ROCrateGenerator('ro-crate-meta.json', 'common_meta.json')
    crate.add_common_metadata("common_meta.json", "Common Metadata", "Test Description")
    crate.add_file("test.fa", "Test FASTA", "Test FASTA File", "application/x-fasta")
    crate.add_dna_analysis("#analysis1", "test.fa", {
        "@type": "DNAAnalysisResult",
        "processId": "TEST-01",
        "sequenceLength": 100,
        "passesAllChecks": True
    })
    crate.add_software("#software", "TestSoft", "1.0", "Test Software")

    json_output = crate.generate_json()
    assert len(json_output["@graph"]) == 6  # 2 initial + 4 added items


def test_duplicate_entries(crate):
    crate.add_file("test.fa", "Test File 1", "Description 1", "application/x-fasta")
    crate.add_file("test.fa", "Test File 2", "Description 2", "application/x-fasta")
    file_entries = [item for item in crate.graph if item["@id"] == "test.fa"]
    assert len(file_entries) == 2  # Allows duplicates, might want to handle this differently
