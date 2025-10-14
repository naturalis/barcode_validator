# Architecture Documentation

This document provides visual representations of the barcode_validator architecture through Mermaid diagrams. The diagrams are organized into modular components suitable for manuscript publication.

## Class Inheritance Hierarchies

### Figure 1: Validator Class Hierarchy

```mermaid
graph TB
    AbstractValidator["AbstractValidator"]
    StructuralValidator["StructuralValidator"]
    TaxonomicValidator["TaxonomicValidator"]
    ProteinCodingValidator["ProteinCodingValidator"]
    NonCodingValidator["NonCodingValidator"]
    
    AbstractValidator --> StructuralValidator
    AbstractValidator --> TaxonomicValidator
    StructuralValidator --> ProteinCodingValidator
    StructuralValidator --> NonCodingValidator
    
    style AbstractValidator fill:#e1f5ff
```

The validator hierarchy implements a Strategy pattern where `AbstractValidator` defines the interface for sequence validation. `StructuralValidator` validates sequence quality (length, ambiguities, marker-specific features), while `TaxonomicValidator` performs reverse taxonomy via identification services. Concrete validators (`ProteinCodingValidator`, `NonCodingValidator`) implement marker-specific validation logic.

### Figure 2: IDService Class Hierarchy

```mermaid
graph TB
    IDService["IDService"]
    BLAST["BLAST"]
    BOLD["BOLD"]
    GalaxyBLAST["GalaxyBLAST"]
    
    IDService --> BLAST
    IDService --> BOLD
    IDService --> GalaxyBLAST
    
    style IDService fill:#e1f5ff
```

The IDService hierarchy provides pluggable identification backends. `BLAST` uses local NCBI BLAST+ searches, `BOLD` queries the BOLD Systems identification engine, and `GalaxyBLAST` uses the Galaxy web service for distributed BLAST searches.

### Figure 3: TaxonResolver Class Hierarchy

```mermaid
graph TB
    TaxonResolver["TaxonResolver"]
    BoldResolver["BoldResolver"]
    NCBIResolver["NCBIResolver"]
    NSRResolver["NSRResolver"]
    
    TaxonResolver --> BoldResolver
    TaxonResolver --> NCBIResolver
    TaxonResolver --> NSRResolver
    
    style TaxonResolver fill:#e1f5ff
```

The TaxonResolver hierarchy handles taxonomic name resolution across different taxonomic backbones. `BoldResolver` parses BOLD process IDs, `NCBIResolver` works with NCBI taxonomy IDs, and `NSRResolver` handles Netherlands Species Register identifiers.

## Validation Process Flow

### Figure 4: Initialization and Setup Phase

```mermaid
flowchart TD
    Start([User runs CLI]) --> CLI[BarcodeValidatorCLI]
    CLI --> ParseArgs[Parse arguments]
    ParseArgs --> CreateOrch[Create ValidationOrchestrator]
    CreateOrch --> InitMode{Check mode}
    
    InitMode -->|structural or both| InitStruct[Initialize StructuralValidator]
    InitMode -->|taxonomic or both| InitTaxon[Initialize TaxonomicValidator]
    
    InitStruct --> SV[StructuralValidator ready]
    InitTaxon --> TV[TaxonomicValidator + IDService ready]
    
    SV --> ParseInput[Parse input file]
    TV --> ParseInput
    ParseInput --> CreateRS[Create DNAAnalysisResultSet]
    CreateRS --> Ready[Ready for validation]
    
    style Start fill:#fff3e1
    style Ready fill:#e1ffe1
```

The initialization phase parses command-line arguments, creates the orchestrator, and initializes validators based on the selected mode (structural, taxonomic, or both). The orchestrator uses factory patterns to instantiate appropriate validator implementations and supporting services (IDService, TaxonResolver, BLASTN, HMM profiles).

### Figure 5: Structural Validation Phase

```mermaid
flowchart TD
    Start([For each sequence]) --> CheckMarker[Check marker type]
    CheckMarker --> ValidateLength[Validate length]
    ValidateLength --> ValidateAmbig[Count ambiguous bases]
    ValidateAmbig --> MarkerSpecific{Marker type}
    
    MarkerSpecific -->|Protein-coding| CheckStopCodons[Check stop codons<br/>via translation]
    MarkerSpecific -->|Non-coding| CalcGC[Calculate GC content]
    
    CheckStopCodons --> RecordResult[Record result]
    CalcGC --> RecordResult
    RecordResult --> End([Next sequence])
    
    style Start fill:#fff3e1
    style End fill:#e1ffe1
```

Structural validation examines sequence quality without external database queries. Each sequence is checked for correct length, ambiguous base content, and marker-specific features. Protein-coding sequences are translated to detect premature stop codons, while non-coding sequences have their GC content calculated.

### Figure 6: Taxonomic Validation Phase

```mermaid
flowchart TD
    Start([Structurally valid sequences]) --> EnrichTaxa[Enrich with expected taxonomy]
    EnrichTaxa --> Batch[Batch sequences by constraint]
    Batch --> IDService{IDService type}
    
    IDService -->|BLAST| RunBLAST[Run BLAST search]
    IDService -->|BOLD| CallBOLD[Query BOLD API]
    IDService -->|Galaxy| CallGalaxy[Call Galaxy service]
    
    RunBLAST --> CollectTaxa[Collect observed taxa]
    CallBOLD --> CollectTaxa
    CallGalaxy --> CollectTaxa
    
    CollectTaxa --> Compare[Compare observed vs expected]
    Compare --> RecordResult[Record taxonomic result]
    RecordResult --> Output[Write results to TSV/FASTA]
    
    style Start fill:#fff3e1
    style Output fill:#e1ffe1
```

Taxonomic validation uses reverse taxonomy to verify specimen identifications. Sequences are enriched with expected taxonomy from the input resolver, then batched by taxonomic constraint for efficient identification. The selected IDService performs searches and returns observed taxa, which are compared against expected taxa to validate the original identification.

