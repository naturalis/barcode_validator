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
  - Validation against ID service reference databases
  - Flexible taxonomy mapping (NSR, NCBI or BOLD)
  - Support for multiple taxonomic ranks

- Input/Output:
  - Support for FASTA and tabular input formats
  - Detailed validation reports
  - Filtered FASTA output for valid sequences
  - Integration with Galaxy workflow platform

## Installation

The easiest way to install a contained environment for the barcode validator is using  
[bioconda](https://pypi.org/project/barcode-validator/):

```bash
conda create -n barcode-validator
conda activate barcode-validator
conda install -c bioconda barcode-validator blast hmmer
```

## Usage

### Command Line Interface

#### Case 1: Validating FASTA against BOLD taxonomy and perform triage for BGE submissions 

```bash

python barcode_validator \
  --input-file <input_fasta> \
  --csv-file <input_csv> \
  --yaml-file <input_yaml> \
  --mode both \
  --marker COI-5P \
  --input-resolver format=bold \
  --input-resolver file=<bold_excel> \
  --output-fasta <output_fasta> \
  --output-tsv <output_tsv> \
  --taxon-validation method=bold \
  --taxon-validation rank=family \
  --taxon-validation min_identity=0.8 \
  --taxon-validation max_target_seqs=100 \
  --triage-config group_id_separator=_ \
  --triage-config group_by_sample=true \
  --log-level ERROR 2> <log_file>
```

### Galaxy Integration

The tool is available as a Galaxy tool wrapper, enabling web-based usage through the Galaxy platform. Users can:
1. Upload sequence files to their Galaxy history
2. Configure validation parameters through the GUI
3. Run validations and view results within Galaxy
4. Download validation reports and filtered sequences

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