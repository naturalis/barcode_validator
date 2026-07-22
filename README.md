[![status](https://joss.theoj.org/papers/62a4e1902c74fc953677a993fb2c1854/status.svg)](https://joss.theoj.org/papers/62a4e1902c74fc953677a993fb2c1854)

# DNA Barcode Validator

A Python-based toolkit for validating DNA barcode sequences through structural and taxonomic validation. This tool 
helps ensure sequence quality and taxonomic accuracy for submissions to the Barcode of Life Data System (BOLD)
and to Naturalis's DNA domain within the BioCloud.

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

Below are detailed examples for common use cases, particularly for the **BGE (Biodiversity Genomics Europe)** and 
**ARISE** projects.

#### BGE Use Case: Structural and Taxonomic Validation with Assembly Triage

The BGE use case involves validating sequences from multiple genome skimming assembly attempts per specimen, selecting 
the best valid sequence per specimen, and performing taxonomic validation using the Galaxy BLAST web service.

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
validation using the BOLD web service. This case is focused on fresh specimen sequencing, either by ONT or Sanger,
which lacks the brute forcing that BGE requires.

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
- `--log-level`: Logging level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`)

### Galaxy Integration

Galaxy features in two distinct ways, which have different requirements. They are easily confused, so they are
described separately here.

#### Barcode Validator as a Galaxy Tool

Structural validation is wrapped as a Galaxy tool and published on the
[Galaxy Toolshed](https://toolshed.g2.bx.psu.edu/) as `rvosa/barcode_validator_structural`, and as the suite
`rvosa/suite_barcode_validator`. This lets non-technical users validate sequences from their Galaxy history
through a web form, and lets the validation be embedded in Galaxy workflows.

Installing a tool from the Toolshed requires **administrator privileges** on the target Galaxy instance; ordinary
users cannot do this themselves. The tool is not preinstalled on the large public servers, so searching for
"barcode validator" in the tool panel of usegalaxy.org or usegalaxy.eu will return nothing. To use it, either ask
the administrator of your institutional instance to install it, or run your own instance.

Once installed:
1. Upload your sequence files to your Galaxy history
2. Configure validation parameters through the GUI
3. Run the validation
4. View results and download validation reports and filtered sequences

#### Galaxy as a Taxonomic Validation Backend

This is a separate feature. `--taxon-validation method=galaxy` does *not* run this package on a Galaxy server.
Instead, it submits your sequences to a different tool, "Identify reads with blastn and find taxonomy", on a
remote Galaxy instance, and parses the taxonomic lineages that come back.

> **Only `galaxy.naturalis.nl` is supported at present.** That tool, and the curated reference databases it
> searches (e.g. "BOLD species only no duplicates"), are deployed on the Naturalis instance only. An API key from
> another server will not work: the run fails with
> `Tool 'Identify reads with blastn and find taxonomy' not found`.

The instance and the credentials are read from the environment:

```bash
export GALAXY_DOMAIN=galaxy.naturalis.nl   # currently the only instance hosting the required tool
export GALAXY_API_KEY=your_galaxy_api_key  # User -> Preferences -> Manage API Key
```

If you do not have an account on galaxy.naturalis.nl, use one of the backends that needs no Galaxy credentials:

- `--taxon-validation method=bold` or `method=bolddistilled`: the BOLD identification service, which needs only an
  internet connection. This is the backend used in the ARISE example above.
- `--taxon-validation method=blast` combined with `--local-blast db=/path/to/db`: a local BLAST+ database that you
  supply yourself.

Accounts on galaxy.naturalis.nl are available by default to Naturalis staff. Furthermore, access can be requested in 
specific cases for scientists connected to institutions that are members of SURF, the Dutch cooperative for academic
computing. If you qualify, you can apply for access via bioinformatics@naturalis.nl

## Architecture

For detailed information about the software architecture, including class hierarchies and process flow diagrams, 
see [Architecture Documentation](ARCHITECTURE.md).

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

If you use this software in your research, please cite this repository. An application note (to JOSS) that 
describes the software is in preparation.

## Contact

For questions, issues, or contributions:
- GitHub Issues: https://github.com/naturalis/barcode_validator/issues
- Email: rutger.vos@naturalis.nl

## Acknowledgments

This tool was developed to support the Biodiversity Genomics Europe (BGE) and ARISE projects, as well as general DNA 
barcoding initiatives at Naturalis Biodiversity Center.
