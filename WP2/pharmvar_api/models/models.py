from typing import List, Dict, Optional, Any
from enums.enums import Function, EvidenceLevel, ReferenceBase, ReferenceCollection, ReferenceLocationType, FilterError
from exceptions.exceptions import ValidationError
from .base_collection import BaseCollection

class Result:
    def __init__(self, status_code: int, message: str = "", data: List[Dict] = None):
        """
        Result returned from low-level RestAdapter
        :param status_code: Standard HTTP Status code
        :param message: Human readable result
        :param data: Python List of Dictionaries (or maybe just a single Dictionary on error)
        """
        self.status_code = int(status_code)
        self.message = str(message)
        self.data = data if data else []


class Variant:
    def __init__(self, hgvs: str = None,
                impact: str = None,
                position: str = None,
                referenceCollections: List[str] = None,
                referenceLocation: str = None,
                referenceSequence: str = None,
                rsId: str = None,
                url: str = None,
                variantFrequency: List[Dict] = None,
                variantId: str = None):
        """
        Represents a variant
        :param hgvs: The HGVS syntax for the variant
        :param impact: The amino acid impact for the variant
        :param position: The Position String for the variant as it is stored in the PharmVar database
        :param reference_collections: List of reference collections
        :param reference_location: The reference location type defining the start position for counting
        :param reference_sequence: The reference sequence for the variant
        :param rs_id: The RS ID Number for the variant
        :param url: The PharmVar URL for the variant
        :param variant_frequency: List of variant frequencies
        :param variant_id: The transient variant ID used for identifying corresponding positions across multiple reference sequences
        """
        self.hgvs = hgvs
        self.impact = impact
        self.position = position
        self.reference_collections = referenceCollections
        self.reference_location = referenceLocation
        self.reference_sequence = referenceSequence
        self.rs_id = rsId
        self.url = url
        self.variant_frequency = variantFrequency
        self.variant_id = variantId

    def __repr__(self) -> str:
        return f"Variant(hgvs={self.hgvs}, impact={self.impact}, position={self.position}, reference_collections={self.reference_collections}, reference_location={self.reference_location}, reference_sequence={self.reference_sequence}, rs_id={self.rs_id}, url={self.url}, variant_frequency={self.variant_frequency}, variant_id={self.variant_id})"

class VariantGroup:
    """
    A representation of a PharmVar Variant Group

    Args:
    hgvsGene (str): The HGVS value for the RefSeq Gene
    hgvsGrch37 (str): The HGVS value for GRCh37
    hgvsGrch38 (str): The HGVS value for GRCh38
    hgvsTranscript (str): The HGVS value for the RefSeq Transcript
    impact (str): The amino acid impact for the variant
    label (str): The group label
    rsIds (List[str]): List of RS IDs
    variantIds (List[str]): List of variant IDs
    """
    def __init__(self,
                hgvsGene: str = None,
                hgvsGrch37: str = None,
                hgvsGrch38: str = None,
                hgvsTranscript: str = None,
                impact: str = None,
                label: str = None,
                rsIds: List[str] = None,
                variantIds: List[str] = None):

        self.hgvs_gene = hgvsGene
        self.hgvs_grch37 = hgvsGrch37
        self.hgvs_grch38 = hgvsGrch38
        self.hgvs_transcript = hgvsTranscript
        self.impact = impact
        self.label = label
        self.rs_ids = rsIds
        self.variant_ids = variantIds

    def __repr__(self) -> str:
        return (f"VariantGroup(label={self.label}, "
                f"impact={self.impact}, "
                f"hgvs_gene={self.hgvs_gene})")

class VariantGroupCollection(BaseCollection):
    def __init__(self, data: List[Dict] = None):
        """
        Represents a list of variant groups
        :param variant_groups: List of VariantGroup objects
        """
        if not data:
            data = []
        self.variant_groups = [VariantGroup(**group) for group in data] if data else []

    def __len__(self) -> int:
        return len(self.variant_groups)

    def __getitem__(self, index: int) -> VariantGroup:
        return self.variant_groups[index]

    def __iter__(self):
        return iter(self.variant_groups)

    def __repr__(self) -> str:
        return f"VariantGroupList(variant_groups={self.variant_groups})"

    def _get_items(self):
        return self.variant_groups
    

class VariantCollection(BaseCollection):
    def __init__(self, data: List[Dict] = None):
        """
        Represents a list of variants
        :param variants: List of Variant objects
        """
        # unpack dictionary if variant data is in dictionary format, otherwise it should be a Variant object
        self.variants = [Variant(**variant) if isinstance(variant, dict) else variant for variant in data] if data else []

    ###todo: create method that gets all the variants with an impact.

    def __len__(self) -> int:
        return len(self.variants)

    def __getitem__(self, index: int) -> Variant:
        return self.variants[index]

    def __iter__(self):
        return iter(self.variants)

    def __repr__(self) -> str:
        return f"VariantList(variants={self.variants})"
    
    def _get_items(self):
        return self.variants

class Allele:
    """
        A representation of a PharmVar allele
        
        Args:
            activeInd (bool): Indicates if the allele is the most recent, active version
            alleleName (str): The name of the allele
            alleleType (str): The allele type ('Core' or 'Sub')
            coreAllele (str): The allele name of the associated core allele for a sub-allele
            description (str): Description for complex alleles
            evidenceLevel (str): The allele evidence level (D=Definitive, M=Moderate, L=Limited)
            function (str): The CPIC Clinical Function for the allele
            geneSymbol (str): The gene symbol for the associated gene
            hgvs (str): The definition of the allele in HGVS notation
            legacyLabel (str): The legacy name for the allele
            pvId (str): The PharmVar ID
            references (List[Dict]): The allele references
            url (str): The PharmVar URL for the allele
            variantGroups (List[Dict]): List of variant groups associated with the allele
            variants (List[Dict]): List of variants
    """
    def __init__(self, 
                 activeInd: bool = False,
                 alleleName: str = None,
                 alleleType: str = None,
                 coreAllele: str = None,
                 description: str = None,
                 evidenceLevel: str = None,
                 function: str = None,
                 geneSymbol: str = None,
                 hgvs: str = None,
                 legacyLabel: str = None,
                 pvId: str = None,
                 references: List[Dict] = None,
                 url: str = None,
                 variantGroups: List[Dict] = None,
                 variants: List[Dict] = None):
        
        self.active_ind = activeInd
        self.allele_name = alleleName
        self.allele_type = alleleType
        self.core_allele = coreAllele
        self.description = description
        self.evidence_level = evidenceLevel
        self.function = function
        self.gene_symbol = geneSymbol
        self.hgvs = hgvs
        self.legacy_label = legacyLabel
        self.pv_id = pvId
        self.references = references
        self.url = url
        self.variant_groups = VariantGroupCollection(variantGroups if variantGroups else [])
        self.variants = VariantCollection(variants if variants else [])

    def __repr__(self) -> str:
        return (f"Allele(allele_name={self.allele_name}, "
                f"allele_type={self.allele_type}, "
                f"gene_symbol={self.gene_symbol}, "
                f"pv_id={self.pv_id})")
        
    
class AlleleCollection(BaseCollection):
    def __init__(self, data: List[Dict] = None):
        """
        Represents a list of alleles
        :param alleles: List of Allele objects
        """
        if not data:
            data = []
        # unpack dictionary if variant data is in dictionary format, otherwise it should be a Variant object
        self.alleles = [Allele(**allele) if isinstance(allele, dict) else allele for allele in data]

    def _get_items(self):
        return self.alleles

    def __len__(self) -> int:
        return len(self.alleles)

    def __getitem__(self, index: int) -> Variant:
        return self.alleles[index]

    def __iter__(self):
        return iter(self.alleles)

    def __repr__(self) -> str:
        return f"AlleleList(alleles={self.alleles})"