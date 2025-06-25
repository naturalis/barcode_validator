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

## Configuration

When setting up the BLAST environment, the following environment variables should be set correctly:

- `BLASTDB`: Path to the BLAST database directory. This must be the *directory* within which the BLAST databases are 
             stored, not the database files themselves. The directory must contain (many) files starting with `nt`. 
             Furthermore, the directory must contain the files `taxdb.btd`, `taxdb.bti`, and `taxonomy4blast.sqlite3`.
             (The `nt.*` files are the indexed sequences, the other files help BLAST running taxonomically constrained 
             queries. All can be fetched into the correct folder using the command `update_blastdb.pl --decompress nt`
             from the NCBI BLAST+ package.)
- `BLASTDB_LMDB_MAP_SIZE`: Optionally, set the size of the LMDB map for the BLAST database. This is useful for large 
             databases and can be set to a value like `1000G` (1 TB RAM) to ensure sufficient RAM for the initial
             map of the BLAST database. At Naturalis, we discovered that this mostly functions as a threshold: if you
             set it too low, BLAST will fail to start. Empirically, this is around 512G. Higher values above the 
             threshold have no effect on performance, they are simply a means to discover you don't have enough RAM 
             available for the BLAST database.


## Usage

### Command Line Interface

#### Case 1: Validating a FASTA file against BOLD taxonomy to generate a tabular report

```bash
python barcode_validator \
  --input_file data/BGE00196_MGE-BGE_r1_1.3_1.5_s50_100.fasta \
  --exp_taxonomy examples/bold.xlsx \
  --exp_taxonomy_type bold \
  --config config/config.yml \
  --output_format tsv \
  --log_level DEBUG > results.tsv
```

- `--input_file`: Path to the input FASTA file containing sequences to validate. The first word in the header line 
  should be the BOLD process ID, followed by an underscore '_', and then a suffix that makes the sequence unique.
- `--exp_taxonomy`: Path to the 'expected taxonomy' file, i.e. what the sequences are expected to be. In this case,
    this is a BOLD spreadsheet in Excel format.
- `--exp_taxonomy_type`: Type of expected taxonomy, either `nsr` (Nederlands Soortenregister) or `bold`. In this case,
   we are using `bold` to validate against the BOLD database.
- `--config`: Path to the configuration file. Almost certainly, you will want to update the `config/config.yml` file
   to specify the BLAST database name, configuration of the BLAST search, and other parameters, and the location of the
   NCBI taxonomy database (as *.tar.gz).
- `--output_format`: Format of the output report. Options are `tsv` (tab-separated values) or `fasta` (filtered FASTA).
   In this case, we are generating a tabular (tsv) report.
- `--log_level`: Set the logging level. Options are `DEBUG`, `INFO`, `WARNING`, `ERROR`, or `CRITICAL`. In this case,
  we are setting it to `DEBUG` for detailed output.
- `> results.tsv`: Redirects the output to a file named `results.tsv`.

Note: the config file has a parameter `blast_db`. This should be set to the name of the BLAST database you want to use. 
The *name* of the database is the path to the 'file stem' of the database, without the `.nhr`, `.nin`, etc. So it is
*not* the name of the directory, but that of the indexed sequence files without the extensions.

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