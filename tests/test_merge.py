import pytest
import sys
import os

from barcode_validator.result import DNAAnalysisResult, DNAAnalysisResultSet

@pytest.fixture
def result_set():
    results = [
        DNAAnalysisResult("BGENL001-23"),
        DNAAnalysisResult("BGENL002-23"),
        DNAAnalysisResult("BGENL003-23"),
        DNAAnalysisResult("BGENL004-23"),
        DNAAnalysisResult("BGENL005-23")
    ]
    return DNAAnalysisResultSet(results)

@pytest.fixture
def data_dir():
    # Construct the path to the data directory
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'examples', 'rich_set')

def test_add_csv_file(result_set, data_dir):
    csv_file = os.path.join(data_dir, 'mge_fastp_r13s100_nocontam.csv')
    result_set.add_csv_file(csv_file)

    # Check if the CSV data was correctly added to each result
    for result in result_set.results:
        assert 'cov_avg' in result.data['ancillary']
        cov_avg = result.data['ancillary']['cov_avg']
        assert cov_avg != ''  # Ensure it's not empty
        try:
            float_cov_avg = float(cov_avg)
            assert str(int(float_cov_avg)).isdigit()  # Check if the integer part is a digit string
        except ValueError:
            pytest.fail(f"cov_avg '{cov_avg}' is not a valid float string")

    # Check a specific value (you may need to adjust this based on your actual data)
    assert result_set.results[0].data['ancillary']['cov_avg'] != ''

def test_add_yaml_file(result_set, data_dir):
    yaml_file = os.path.join(data_dir, 'mge_fastp_r13s100_nocontam.yaml')
    result_set.add_yaml_file(yaml_file)

    # Check if the YAML data was correctly added to each result
    expected_keys = [
        'samples_file', 'protein_reference_file', 'output_dir', 'genes', 'r', 's', 'run_name'
    ]

    for result in result_set.results:
        for key in expected_keys:
            assert key in result.data['ancillary'], f"Key '{key}' not found in ancillary data"

    # Check specific values
    first_result = result_set.results[0].data['ancillary']
    assert first_result['samples_file'] == "/gpfs/nhmfsa/bulk/share/data/mbl/share/workspaces/groups/genomics-collections/BGE/2024-02-01/DataDelivery_2024-02-01_18-24-39_snpseq00629/files/WK-3860/read_paths.csv"
    assert first_result['protein_reference_file'] == "/gpfs/nhmfsa/bulk/share/data/mbl/share/scratch/MGE/protein_references/benchmarking_data_570_refs-contam_refs_final14/benchmarking_data_570_taxonomy_gene_fetch_sum_out.csv"
    assert first_result['output_dir'] == "/gpfs/nhmfsa/bulk/share/data/mbl/share/scratch/MGE/cox1/benchmarking/MGE-fastp_pipeline_alt_params"
    assert first_result['genes'] == ['cox1']
    assert first_result['r'] == [1.3]
    assert first_result['s'] == [100]
    assert first_result['run_name'] == "mge_fastp_r13s100_nocontam"


def test_result_fields(result_set, data_dir):
    # Add CSV file
    csv_file = os.path.join(data_dir, 'mge_fastp_r13s100_nocontam.csv')
    result_set.add_csv_file(csv_file)

    # Add YAML file
    yaml_file = os.path.join(data_dir, 'mge_fastp_r13s100_nocontam.yaml')
    result_set.add_yaml_file(yaml_file)

    # Check result fields
    fields = result_set.results[0].result_fields()
    assert 'cov_avg' in fields
    assert 'samples_file' in fields
    assert 'protein_reference_file' in fields
    assert 'output_dir' in fields
    assert 'genes' in fields
    assert 'r' in fields
    assert 's' in fields
    assert 'run_name' in fields

    # Check values for a single result
    values = result_set.results[0].get_values()
    assert any(isinstance(v, str) and v != '' and float(v) > 0 for v in values)  # cov_avg from CSV
    assert any(
        v == "/gpfs/nhmfsa/bulk/share/data/mbl/share/workspaces/groups/genomics-collections/BGE/2024-02-01/DataDelivery_2024-02-01_18-24-39_snpseq00629/files/WK-3860/read_paths.csv"
        for v in values)  # samples_file from YAML

    # Check string representation of the entire result set
    result_set_str = str(result_set)

    # Split the string into lines
    lines = result_set_str.split('\n')

    # Check header
    header = lines[0].split('\t')
    assert 'cov_avg' in header
    assert 'samples_file' in header
    assert 'genes' in header
    assert 'r' in header
    assert 's' in header
    assert 'run_name' in header

    # Check content (first data line)
    first_data_line = lines[1].split('\t')
    assert "/gpfs/nhmfsa/bulk/share/data/mbl/share/workspaces/groups/genomics-collections/BGE/2024-02-01/DataDelivery_2024-02-01_18-24-39_snpseq00629/files/WK-3860/read_paths.csv" in first_data_line
    assert "['cox1']" in first_data_line
    assert "[1.3]" in first_data_line
    assert "[100]" in first_data_line
    assert "mge_fastp_r13s100_nocontam" in first_data_line

if __name__ == '__main__':
    pytest.main()