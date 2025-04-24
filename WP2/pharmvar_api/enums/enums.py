from enum import Enum

class Function(str, Enum):
    DECREASED = "decreased function"
    NOT_ASSIGNED = "function not assigned"
    INCREASED = "increased function"
    NORMAL = "normal function"  
    POSSIBLY_DECREASED = "possibly decreased"
    SEVERELY_DECREASED = "severely decreased"
    UNCERTAIN = "uncertain function"
    UNKNOWN = "unknown function"

class EvidenceLevel(str, Enum):
    DEFINITIVE = "Definitive"
    LIMITED = "Limited"
    MODERATE = "Moderate"

class ReferenceBase(str, Enum):
    A = "A"
    C = "C"
    G = "G"
    T = "T"

class ReferenceCollection(str, Enum):
    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"
    REFSEQ_GENE = "RefSeqGene"
    REFSEQ_TRANSCRIPT = "RefSeqTranscript"

class ReferenceLocationType(str, Enum):
    ATG_START = "ATG Start"
    SEQUENCE_START = "Sequence Start"

class FilterError(Enum):
    INVALID_ATTRS = "Invalid attribute(s): {}"
    TYPE_MISMATCH = "Type mismatch for {}: expected {}, got {}"