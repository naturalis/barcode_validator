import json
from nbitk.Taxon import Taxon
from typing import Dict, List, Any
from barcode_validator.result import DNAAnalysisResult


class DNAAnalysisEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, DNAAnalysisResult):
            return {
                "@type": "DNAAnalysisResult",
                "processId": obj.process_id,
                "markerLength": obj.seq_length,
                "fullLength": obj.full_length,
                "observedTaxon": [str(taxon.name) for taxon in obj.obs_taxon],
                "expectedTaxon": str(obj.exp_taxon.name),
                "species": str(obj.species.name),
                "stopCodons": obj.stop_codons,
                "markerAmbiguities": obj.ambiguities,
                "fullAmbiguities": obj.full_ambiguities,
                "passesAllChecks": obj.passes_all_checks(),
                "barcodeRank": obj.calculate_ranks()[0],
                "fullRank": obj.calculate_ranks()[1]
            }
        return super().default(obj)


class ROCrateGenerator:
    def __init__(self, location: str, meta_location: str):
        self.context = "https://w3id.org/ro/crate/1.1/context"
        self.graph = [
            {
                "@type": "CreativeWork",
                "@id": location,
                "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"},
                "about": {"@id": "./"}
            },
            {
                "@id": "./",
                "@type": "Dataset",
                "additionalType": "https://schema.org/Dataset",
                "mainEntity": {"@id": meta_location}
            }
        ]
        self.encoder = DNAAnalysisEncoder()

    def add_common_metadata(self, location: str, name: str, description: str):
        common_metadata = {
            "@id": location,
            "@type": "CreativeWork",
            "name": name,
            "description": description,
            "encodingFormat": "application/json"
        }
        self.graph.append(common_metadata)

    def add_file(self, file_id: str, name: str, description: str, encoding_format: str):
        file_metadata = {
            "@id": file_id,
            "@type": "File",
            "name": name,
            "description": description,
            "encodingFormat": encoding_format
        }
        self.graph.append(file_metadata)

    def add_dna_analysis(self, idref: str, object_id: str, result: DNAAnalysisResult):
        dna_analysis = {
            "@id": idref,
            "@type": "CreativeWork",
            "name": "DNA Sequence Analysis",
            "description": "Analysis results of the DNA sequence from the FASTA file",
            "object": {"@id": object_id},
            "result": json.loads(self.encoder.encode(result))  # Use the encoder here
        }
        self.graph.append(dna_analysis)

    def add_multiple_dna_analyses(self, analyses: List[Dict[str, Any]]):
        for analysis in analyses:
            self.add_dna_analysis(analysis['idref'], analysis['object_id'], analysis['result'])

    def add_software(self, idref: str, name: str, version: str, description: str):
        software = {
            "@id": idref,
            "@type": "SoftwareApplication",
            "name": name,
            "version": version,
            "description": description
        }
        self.graph.append(software)

    def generate_json(self) -> Dict[str, Any]:
        return {
            "@context": self.context,
            "@graph": self.graph
        }

    def save_to_file(self, filename: str):
        with open(filename, 'w') as f:
            json.dump(self.generate_json(), f, indent=2, cls=DNAAnalysisEncoder)


# Example usage
if __name__ == "__main__":
    crate = ROCrateGenerator('ro-crate-meta.json', 'common_meta.json')

    crate.add_common_metadata(
        "common_meta.json",
        "Common Metadata",
        "Common metadata for the DNA analysis project"
    )

    crate.add_file(
        "mge.fa",
        "Sample FASTA File",
        "FASTA file containing the DNA sequence that was analyzed",
        "application/x-fasta"
    )

    result = DNAAnalysisResult("BGENL191-23")
    result.seq_length = 658
    result.obs_taxon = [Taxon(name="Trichoptera", taxonomic_rank="family")]
    result.exp_taxon = Taxon(name="Trichoptera", taxonomic_rank="family")
    result.species = Taxon(name="Halesus tessellatus", taxonomic_rank="species")
    result.stop_codons = []
    result.ambiguities = 0
    result.full_length = 1509
    result.full_ambiguities = 4

    crate.add_dna_analysis(
        "#analysis1",
        "mge.fa",
        result
    )

    crate.add_software(
        "#software",
        "barcode_validator",
        "1.0",
        "Validates DNA sequences against a reference database"
    )

    crate.save_to_file("ro-crate-metadata.json")