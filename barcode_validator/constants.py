from enum import Enum

class RefDB(Enum):
    NCBI = "ncbi"
    BOLD = "bold"

class TaxonomicRank(Enum):
    KINGDOM = "kingdom"
    PHYLUM = "phylum"
    CLASS = "class"
    ORDER = "order"
    FAMILY = "family"
    GENUS = "genus"
    SPECIES = "species"

class TaxonomicBackbone(Enum):
    NSR = "nsr"    # NSR taxonomy from DwC-A
    BOLD = "bold"  # BOLD taxonomy from Excel
    NCBI = "ncbi"  # NCBI taxonomy from taxdump

class Marker(Enum):
    COI_5P = "COI-5P"
    MATK = "matK"
    RBCL = "rbcL"
    ITS = "ITS"
    ITS2 = "ITS2"

class ValidationMode(Enum):
    STRUCTURAL = "structural"
    TAXONOMIC = "taxonomic"
    BOTH = "both"