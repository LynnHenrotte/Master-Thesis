from enum import Enum

class VariantEndpoint(Enum):
    ALL = "variants"
    ALLELE = "variants/allele/{identifier}"
    GENE = "variants/gene/{symbol}"
    RSID = "variants/rsid/{rsId}"
    RSID_IMPACT = "variants/rsid/{rsId}/impact"
    RSID_FREQUENCY = "variants/rsid/{rsId}/variant-frequency"
    SPDI = "variants/spdi/{spdi}"   
    SPDI_IMPACT = "variants/spdi/{spdi}/impact"
    SPDI_FREQUENCY = "variants/spdi/{spdi}/variant-frequency"


class AlleleEndPoint(Enum):
    ALL = "alleles"
    ACTIVE = "alleles/list"
    IDENTIFIER = "alleles/{identifier}" # returns a single allele by identifier (PV_ID or allele name)
    ALLELE_NAME = "alleles/{identifier}/allele-name"
    EVIDENCE_LEVEL = "alleles/{identifier}/evidence_level"
    function = "alleles/{identifier}/function"
    PV_ID = "alleles/{identifier}/pvid"
    REFERENCES = "alleles/{identifier}/references"
    VARIANTS = "alleles/{identifier}/variants"


