# DNA Barcode Validator

A Python-based toolkit for validating DNA barcode sequences through structural and taxonomic validation. This tool 
helps ensure sequence quality and taxonomic accuracy for submissions to the Barcode of Life Data System (BOLD)
and to Naturalis's Core Sequence Cloud.

## Features

- Structural validation of DNA barcodes:
  - Sequence length requirements
  - Ambiguous base detection
  - Stop codon analysis for protein-coding markers
  - HMM-based alignment for codon phase detection

- Taxonomic validation:
  - BLAST-based validation against reference databases
  - Flexible taxonomy mapping (NSR or BOLD)
  - Integration with NCBI taxonomy
  - Support for multiple taxonomic ranks

- Input/Output:
  - Support for FASTA and tabular input formats
  - Detailed validation reports
  - Filtered FASTA output for valid sequences
  - Integration with Galaxy workflow platform

## Installation

### Using Conda

```bash
# Create and activate conda environment
conda env create -f environment.yml
conda activate barcode-validator

# Installation of Python dependencies invoked by conda
# pip install -r requirements.txt
```

### Manual Installation

```bash
# Install required command line tools
sudo apt-get install hmmer ncbi-blast+

# Install Python package
pip install .
```

## Usage

### Command Line Interface

```bash
# Using NSR taxonomy, we do this when we validate CSC dumps
python barcode_validator \
  --input_file sequences.tsv \
  --exp_taxonomy nsr.zip \
  --exp_taxonomy_type nsr \
  --config config.yml \
  --output_file results.tsv

# Using BOLD taxonomy, this is the typical process for BGE where we prepare BOLD uploads
barcode-validator \
  --input_file sequences.fasta \
  --exp_taxonomy bold.xlsx \
  --exp_taxonomy_type bold \
  --config config.yml \
  --output_file results.tsv \
  --emit_valid_fasta --output_fasta valid.fasta
```

### Galaxy Integration

The tool is available as a Galaxy tool wrapper, enabling web-based usage through the Galaxy platform. Users can:
1. Upload sequence files to their Galaxy history
2. Configure validation parameters through the GUI
3. Run validations and view results within Galaxy
4. Download validation reports and filtered sequences

## Configuration

The tool uses YAML configuration files for flexible setup:

```yaml
# Example config.yml
marker: COI-5P
validation_rank: family
taxonomic_backbone: bold
blast_db: /path/to/blast/db
hmm_profile_dir: /path/to/hmm/profiles
log_level: INFO
```

See `config/config.yml` for a complete configuration template with documentation.

## Contributing

We welcome contributions! Please see:
- [Contributing Guidelines](CONTRIBUTING.md)
- [Code of Conduct](CODE_OF_CONDUCT.md)

## Testing

```bash
# Run test suite
pytest

# Run with coverage
pytest --cov=barcode_validator
```

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE.md](LICENSE.md) file for details.

## Citation

If you use this software in your research, please cite:

[Citation information to be added]

## Contact

[Contact information to be added]