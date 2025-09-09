from enum import Enum

# This file maintains constants that are passed around through the application.
# The general flow is that the variables that users provided are anchored on the
# schema.yaml, and prepared and validated as command line arguments by schema_config.py.
# There, cli.py picks them up from user input on the command line and passes them on
# to orchestrator.py. The orchestrator then translates the 'dirty' input from users
# into the constants that are defined here. Beyond that, other classes should therefore
# NOT by dealing with dirty strings or config variables (ok, beyond logging verbosity).

class RefDB(Enum):
    BLAST = "blast"
    BOLD = "bold"
    GALAXY = "galaxy"
    BOLDDISTILLED = "bolddistilled"

class TaxonomicRank(Enum):
    KINGDOM = "kingdom"
    PHYLUM = "phylum"
    CLASS = "class"
    ORDER = "order"
    FAMILY = "family"
    SUBFAMILY = "subfamily"
    GENUS = "genus"
    SPECIES = "species"
    SUBSPECIES = "subspecies"

def index_for_rank(rank: TaxonomicRank) -> int:
    """Return the index for a given taxonomic rank."""
    rank_order = [
        TaxonomicRank.KINGDOM,
        TaxonomicRank.PHYLUM,
        TaxonomicRank.CLASS,
        TaxonomicRank.ORDER,
        TaxonomicRank.FAMILY,
        TaxonomicRank.SUBFAMILY,
        TaxonomicRank.GENUS,
        TaxonomicRank.SPECIES,
        TaxonomicRank.SUBSPECIES
    ]
    return rank_order.index(rank)

def rank_by_index(index: int) -> TaxonomicRank:
    """Return the taxonomic rank for a given index."""
    rank_order = [
        TaxonomicRank.KINGDOM,
        TaxonomicRank.PHYLUM,
        TaxonomicRank.CLASS,
        TaxonomicRank.ORDER,
        TaxonomicRank.FAMILY,
        TaxonomicRank.SUBFAMILY,
        TaxonomicRank.GENUS,
        TaxonomicRank.SPECIES,
        TaxonomicRank.SUBSPECIES
    ]
    if 0 <= index < len(rank_order):
        return rank_order[index]
    else:
        raise ValueError(f"Index {index} is out of range for taxonomic ranks.")

class TaxonomicBackbone(Enum):
    DWC = "dwc"    # NSR taxonomy from DwC-A
    BOLD = "bold"  # BOLD taxonomy from Excel
    NCBI = "ncbi"  # NCBI taxonomy from taxdump
    BOLDDISTILLED = "bolddistilled"  # BOLD distilled taxonomy from TSV

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