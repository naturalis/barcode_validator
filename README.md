# DNA Barcode Validator

A Python-based toolkit for validating DNA barcode sequences through structural and taxonomic validation. This tool 
helps ensure sequence quality and taxonomic accuracy for submissions to the Barcode of Life Data System (BOLD)
and to Naturalis's Core Sequence Cloud.

## Features

- **Structural validation** of DNA barcodes:
  - Sequence length requirements
  - Ambiguous base detection
  - Stop codon analysis for protein-coding markers
  - HMM-based alignment for codon phase detection

- **Taxonomic validation**:
  - Validation against ID service reference databases (BOLD, Galaxy BLAST)
  - Flexible taxonomy mapping (NSR, NCBI or BOLD)
  - Support for multiple taxonomic ranks

- **Triaging and filtering**:
  - Automatic selection of best valid sequence per specimen or sample group
  - Support for assembly attempt grouping
  - Customizable triage criteria

- **Input/Output**:
  - Support for FASTA and tabular input formats
  - BOLD Excel spreadsheet integration
  - Detailed validation reports in TSV format
  - Filtered FASTA output for valid sequences
  - Integration with Galaxy workflow platform

## Installation

### Using pip

Install the barcode validator from PyPI using pip:

```bash
pip install barcode-validator
```

Note: Additional dependencies (BLAST and HMMER) may need to be installed separately depending on your use case.

### Using bioconda

The recommended way to install a complete environment with all dependencies is using 
[bioconda](https://bioconda.github.io/):

```bash
conda create -n barcode-validator
conda activate barcode-validator
conda install -c bioconda barcode-validator blast hmmer
```

This will install the barcode validator along with BLAST and HMMER, which are required for taxonomic and structural 
validation respectively.

## Usage

### Command Line Interface

The barcode validator can be run as a Python module:

```bash
python -m barcode_validator [options]
```

Below are detailed examples for common use cases, particularly for the **BGE (Bioscan Genomics Europe)** and 
**ARISE** projects.

#### BGE Use Case: Structural and Taxonomic Validation with Assembly Triage

The BGE use case involves validating sequences from multiple assembly attempts per specimen, selecting the best valid 
sequence per specimen, and performing taxonomic validation using the Galaxy BLAST web service.

**Input requirements:**
- FASTA file with sequences where IDs are formatted as `processID_assemblyAttemptID`
- CSV file with assembly metrics (optional but recommended)
- BOLD Excel spreadsheet with 'Lab Sheet' and 'Taxonomy' tabs
- Galaxy API credentials for taxonomic validation

**Example: Two-stage validation**

First, perform structural validation with early triage:

```bash
# Set up input and output files
INPUT_FASTA=data/sequences.fasta
INPUT_CSV=data/metrics.csv
BOLD_EXCEL=data/bold_spreadsheet.xlsx
STRUCTVAL_FASTA=data/structval_out.fasta
STRUCTVAL_TSV=data/structval_out.tsv

# Run structural validation
python -m barcode_validator \
  --input-file $INPUT_FASTA \
  --csv-file $INPUT_CSV \
  --mode structural \
  --marker COI-5P \
  --input-resolver format=bold \
  --input-resolver file=$BOLD_EXCEL \
  --output-fasta $STRUCTVAL_FASTA \
  --output-tsv $STRUCTVAL_TSV \
  --triage-config group_id_separator=_ \
  --triage-config group_by_sample=true \
  --log-level INFO 2> structval.log
```

Then, perform taxonomic validation on the triaged results:

```bash
# Set Galaxy credentials
export GALAXY_API_KEY=your_galaxy_api_key
export GALAXY_DOMAIN=galaxy.naturalis.nl

# Set output files
TAXVAL_FASTA=data/taxonval_out.fasta
TAXVAL_TSV=data/taxonval_out.tsv

# Run taxonomic validation
python -m barcode_validator \
  --input-file $STRUCTVAL_FASTA \
  --csv-file $INPUT_CSV \
  --mode taxonomic \
  --marker COI-5P \
  --input-resolver format=bold \
  --input-resolver file=$BOLD_EXCEL \
  --output-fasta $TAXVAL_FASTA \
  --output-tsv $TAXVAL_TSV \
  --taxon-validation method=galaxy \
  --taxon-validation rank=family \
  --taxon-validation min_identity=0.8 \
  --taxon-validation max_target_seqs=100 \
  --log-level INFO 2> taxonval.log
```

**Example: Combined validation in one step**

For more thorough validation where all structurally valid sequences are checked taxonomically before triaging:

```bash
python -m barcode_validator \
  --input-file $INPUT_FASTA \
  --csv-file $INPUT_CSV \
  --mode both \
  --marker COI-5P \
  --input-resolver format=bold \
  --input-resolver file=$BOLD_EXCEL \
  --output-fasta data/validated_out.fasta \
  --output-tsv data/validated_out.tsv \
  --taxon-validation method=galaxy \
  --taxon-validation rank=family \
  --taxon-validation min_identity=0.8 \
  --taxon-validation max_target_seqs=100 \
  --triage-config group_id_separator=_ \
  --triage-config group_by_sample=true \
  --log-level INFO 2> validation.log
```

#### ARISE Use Case: Simple Validation Without Assembly Grouping

The ARISE use case involves validating individual sequences (one per specimen) with both structural and taxonomic 
validation using the BOLD web service.

**Input requirements:**
- FASTA file with sequences where the first word of the definition line is the process ID
- BOLD Excel spreadsheet with 'Lab Sheet' and 'Taxonomy' tabs

**Example:**

```bash
python -m barcode_validator \
  --input-file input_sequences.fasta \
  --mode both \
  --marker COI-5P \
  --input-resolver format=bold \
  --input-resolver file=bold_spreadsheet.xlsx \
  --output-fasta validated_sequences.fasta \
  --output-tsv validation_report.tsv \
  --taxon-validation method=bold \
  --taxon-validation rank=family \
  --taxon-validation min_identity=0.8 \
  --taxon-validation max_target_seqs=100 \
  --triage-config group_by_sample=false \
  --log-level ERROR 2> validation.log
```

This produces a FASTA file with valid sequences and a TSV file with detailed validation results.

#### Common Options

- `--input-file`: Input FASTA file with sequences to validate
- `--csv-file`: Optional CSV file with additional metrics
- `--mode`: Validation mode (`structural`, `taxonomic`, or `both`)
- `--marker`: Marker gene (e.g., `COI-5P`)
- `--input-resolver format=bold`: Use BOLD spreadsheet for taxonomy
- `--input-resolver file=<path>`: Path to BOLD Excel spreadsheet
- `--output-fasta`: Output FASTA file with valid sequences
- `--output-tsv`: Output TSV file with validation results
- `--taxon-validation method=<bold|galaxy>`: Taxonomic validation service
- `--taxon-validation rank=<rank>`: Taxonomic rank to validate at
- `--taxon-validation min_identity=<float>`: Minimum identity threshold (0-1)
- `--taxon-validation max_target_seqs=<int>`: Maximum number of BLAST hits to consider
- `--triage-config group_id_separator=<char>`: Separator for parsing group IDs
- `--triage-config group_by_sample=<true|false>`: Enable/disable assembly attempt grouping
- `--log-level`: Logging level (`ERROR`, `INFO`, `DEBUG`)

### Galaxy Integration

The structural validation functionality is also available through the Galaxy platform:

- **Galaxy Toolshed**: The barcode validator structural validation tool is available in the Galaxy Toolshed, 
  enabling easy installation into any Galaxy instance.
- **Web-based interface**: Users can upload sequence files, configure validation parameters through the GUI, 
  run validations, and download results.
- **Workflow integration**: The tool can be incorporated into Galaxy workflows for automated processing pipelines.

To use the tool in Galaxy:
1. Install the tool from the Galaxy Toolshed (search for "barcode validator")
2. Upload your sequence files to your Galaxy history
3. Configure validation parameters through the GUI
4. Run the validation
5. View results and download validation reports and filtered sequences

For taxonomic validation through Galaxy, the tool can connect to Galaxy's BLAST web service using API credentials.

## Contributing

We welcome contributions! Please see:
- [Contributing Guidelines](CONTRIBUTING.md)
- [Code of Conduct](CODE_OF_CONDUCT.md)

When contributing, please ensure:
- Code follows PEP 8 style guidelines
- All functions and classes include docstrings
- New features include unit tests
- All tests pass before submitting a pull request

## Testing

Run the test suite to verify your installation:

```bash
# Run all tests
pytest

# Run with coverage report
pytest --cov=barcode_validator

# Run specific test suites
pytest tests/bge/
pytest tests/arise/
```

The test suite includes comprehensive examples for BGE and ARISE use cases.

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE.md](LICENSE.md) file for details.

## Citation

If you use this software in your research, please cite:

[Citation information to be added]

## Contact

For questions, issues, or contributions:
- GitHub Issues: https://github.com/naturalis/barcode_validator/issues
- Email: rutger.vos@naturalis.nl

## Acknowledgments

This tool was developed to support the Bioscan Genomics Europe (BGE) and ARISE projects, as well as general DNA 
barcoding initiatives at Naturalis Biodiversity Center.