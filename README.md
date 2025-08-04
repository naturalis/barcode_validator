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

The easiest way to install a contained environment for the barcode validator is using  
[bioconda](https://pypi.org/project/barcode-validator/):

```bash
conda create -n barcode-validator
conda activate barcode-validator
conda install -c bioconda barcode-validator blast hmmer
```

## Usage

### Command Line Interface

#### Case 1: Validating a FASTA file against BOLD taxonomy to generate a tabular report

```bash
python barcode_validator \
  --input_file data/BGE00196_MGE-BGE_r1_1.3_1.5_s50_100.fasta \
  --exp_taxonomy examples/bold.xlsx \
  --exp_taxonomy_type bold \
  --output_format tsv \
  --mode both \
  --log_level WARNING \
  > results.tsv \
  2> results.log
```

- `--input_file`: Path to the input FASTA file containing sequences to validate. The first word in the header line 
  should be the BOLD process ID, followed by an underscore '_', and then a suffix that makes the sequence unique.
  (The underscore separator can be changed in the configuration file under `group_id_separator`).
- `--exp_taxonomy`: Path to the 'expected taxonomy' file, i.e. what the sequences are expected to be. In this case,
    this is a BOLD spreadsheet in Excel format.
- `--exp_taxonomy_type`: Type of expected taxonomy, either `nsr` (Nederlands Soortenregister) or `bold`. In this case,
   we are using `bold` to validate against the BOLD database.
- `--output_format`: Format of the output report. Options are `tsv` (tab-separated values) or `fasta` (filtered FASTA).
   In this case, we are generating a tabular (tsv) report.
- `--mode`: Validation mode. Options are `structural`, `taxonomic`, or `both`. In this case, we are validating both 
  structural and taxonomic aspects of the sequences.
- `--log_level`: Set the logging level. Options are `DEBUG`, `INFO`, `WARNING`, `ERROR`, or `CRITICAL`. In this case,
  we are setting it to `WARNING` to only log warnings and errors.
- `> results.tsv`: Redirects the output to a file named `results.tsv`.
- `2> results.log`: Redirects warnings and error messages to a file named `results.log`.

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