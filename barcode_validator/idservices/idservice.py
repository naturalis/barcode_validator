from typing import Set
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.BaseTree import Tree
from nbitk.Taxon import Taxon
from barcode_validator.taxonomy_resolver import TaxonomicRank


class IDService():
    """
    Abstract base class for identification services.
    
    This class defines the interface that all identification services must implement.
    Identification services are used to identify taxonomic information from sequence
    records, typically by comparing them against reference databases.
    """
    
    def identify_record(self, record: SeqRecord, level: TaxonomicRank = TaxonomicRank.FAMILY, extent: Taxon = None) -> Set[Taxon]:
        """
        Identify the taxonomic classification of a sequence record. This method returns the set of
        higher taxon objects (at the specified taxonomic rank) that subtend the specific matches
        obtained by the underlying identification service. By default, this is a set of families.
        
        Optionally, the extent of the search can be constrained by providing a Taxon object, in
        which case the search will be limited to that higher taxon. Typically, this will be 
        something like the taxonomic class within which the record is expected to be found. This
        feature is provided as a means to limit the search space. 
        
        For example, when running identifications via BLAST, the `extent` is the taxon (with NCBI 
        taxon ID) that is passed to the `-taxids` parameter. Necessarily, the `extent` must be a 
        taxon at a higher level than the `level` at which the identification is being made.
        Note that the BLAST runner uses Metazoa (taxon:33208) as the default extent.
        
        :param record: A Bio.SeqRecord object containing the sequence to identify
        :param level: The taxonomic rank at which to return results (default: 'family')
        :param extent: A Taxon object representing the extent of the search (default: None)
        :return: A list of Taxon objects representing the possible taxonomic classifications
        """
        pass
