from typing import Optional
import logging
from exceptions.exceptions import InvalidArgumentError
from api import RestAdapter, VariantEndpoint, AlleleEndPoint
from models import Variant, VariantCollection, AlleleCollection
from config import APIConfig
class PharmVarApi:
    """
    PharmVarApi is a class that provides methods to interact with the PharmVar API.
    It uses a RestAdapter object to make HTTP requests to the API server.
    The class provides methods to get variants, alleles and genes from the PharmVar database.
    """    
    # Initialize the PharmVarApi object with a RestAdapter object
    def __init__(self, hostname: str = APIConfig.DEFAULT_HOST, api_key: str = '', ver: str = APIConfig.DEFAULT_VERSION, logger: logging.Logger = APIConfig.DEFAULT_LOGGER):
        self._rest_adapter = RestAdapter(hostname, api_key, ver, logger)
    
    # Variant methods

    def get_all_variants(self) -> VariantCollection:

        ### TODO: add parameters
        """
        Get all variants from the PharmVar database.
        return: VariantCollection object containing all variants
        """

        result = self._rest_adapter._do(http_method = 'GET', endpoint = VariantEndpoint.ALL.value)
        return VariantCollection(data = result.data)


    def get_variants_by_gene(self, gene_symbol: str) -> VariantCollection:
        """
        Get all variants for a given gene symbol.
        :param gene_symbol: The gene symbol to search for
        :return: VariantCollection object containing all variants for the gene
        """
        endpoint = VariantEndpoint.GENE.value.format(symbol = gene_symbol)
        result = self._rest_adapter._do(http_method = 'GET', endpoint = endpoint)
        return VariantCollection(data = result.data)

    def get_variants_by_allele(self, identifier: str) -> VariantCollection:
        """
        Get all variants for a given identifier.
        :param identifier: can be either a PharmVar ID or an allele name.
        :return: VariantCollection object containing all variants for the identifier.
        """
        endpoint = VariantEndpoint.ALLELE.value.format(identifier = identifier)
        result = self._rest_adapter._do(http_method = 'GET', endpoint = endpoint)
        return VariantCollection(data = result.data)
    
    def get_variants_by_rsid(self, rs_id: str) -> VariantCollection:
        """
        Get all variants for a given rsID.
        :param rs_id: The rsId of the variant
        :return: VariantCollection object containing all variants for the rsID
        """
        endpoint = VariantEndpoint.RSID.value.format(rsId = rs_id)
        result = self._rest_adapter._do(http_method = 'GET', endpoint = endpoint)
        return VariantCollection(data = result.data)
    
    def get_variants_by_spdi(self, spdi: str) -> VariantCollection:
        """
        Get all variants for a given SPDI.
        :param spdi: the SPDI designation of a variant: Reference Sequence:Position:Deletion:Insertion
        :return: VariantCollection object containing all variants for the SPDI
        """
        endpoint = VariantEndpoint.SPDI.value.format(spdi = spdi)
        result = self._rest_adapter._do(http_method = 'GET', endpoint = endpoint)
        return VariantCollection(data = result.data)
    
    def get_variant_impact(self, *, rs_id: str = None, spdi: str = None) -> str:
        """
        Get the impact of a variant for a given rsID or SPDI.
        :param rs_id: The rsId of the variant
        :param spdi: the SPDI designation of a variant: Reference Sequence:Position:Deletion:Insertion
        :return: VariantCollection object containing the impact of the variant for the rsID
        """
        if rs_id is not None and spdi is not None:
            raise InvalidArgumentError("Only one of rs_id or spdi can be provided, not both")

        if rs_id is not None:
            endpoint = VariantEndpoint.RSID_IMPACT.value.format(rsId=rs_id)
        elif spdi is not None:
            endpoint = VariantEndpoint.SPDI_IMPACT.value.format(spdi= spdi)
        else:
            raise InvalidArgumentError("Either rs_id or spdi must be provided.")

        result = self._rest_adapter._do(http_method = "GET", endpoint = endpoint, headers={"Accept": "*/*"})
        # Save in Variant instance with rs_id or spdi and impact
        allele = Variant(impact = result.data, rsId = rs_id,
                        referenceSequence = spdi.split(":")[0] if spdi else None,
                        position = spdi.split(":")[1] if spdi else None
                        )
        return allele.impact
    
        # Allele methods

    def get_all_alleles(
            self,
            exclude_sub_alleles: bool = False,
            function: Optional[str] = None,
            include_reference_variants: bool = False,
            include_retired_alleles: bool = False,
            include_retired_reference_sequences: bool = False,
            min_evidence_level: Optional[str] = None,
            position: Optional[int] = None,
            reference_base_sequence: Optional[str] = None,
            reference_collection: Optional[str] = None,
            reference_location_type: Optional[str] = None,
            reference_sequence: Optional[str] = None,
            variant_base_sequence: Optional[str] = None
        ) -> AlleleCollection:
            """
            Get all alleles from the PharmVar database with optional filtering parameters.

            Parameters:
                exclude_sub_alleles (bool): Exclude sub-allele definitions from results. Default is False.
                function (Function): Filter results by function (decreased function, function not assigned, increased function, normal function, possibly decreased, severely decreased, uncertain function, unknown function).
                include_reference_variants (bool): Include reference variants like A>A or C>C. Default is False.
                include_retired_alleles (bool): Include retired allele definitions from previous versions. Default is False.
                include_retired_reference_sequences (bool): Include variants from retired reference sequences. Default is False.
                min_evidence_level (EvidenceLevel): Filter by minimum evidence level (Definitive, Limited, Moderate).
                position (int): Filter results by variant position.
                reference_base_sequence (ReferenceBase): Filter by reference base sequence (A, C, G, T).
                reference_collection (ReferenceCollection): Filter by reference collection (GRCh37, GRCh38, etc.).
                reference_location_type (ReferenceLocationType): Filter by reference starting location.
                reference_sequence (str): Filter by reference sequence (e.g., "NG_008376.4").
                variant_base_sequence (ReferenceBase): Filter by observed variant base sequence (A, C, G, T).

            Returns:
                AlleleCollection: Collection of alleles matching the specified criteria
            """
            params = {
                "exclude-sub-alleles": exclude_sub_alleles,
                "include-reference-variants": include_reference_variants,
                "include-retired-alleles": include_retired_alleles,
                "include-retired-reference-sequences": include_retired_reference_sequences
            }

            # Add optional parameters if they are provided
            if function:
                params["function"] = function
            if min_evidence_level:
                params["min-evidence-level"] = min_evidence_level
            if position:
                params["position"] = position
            if reference_base_sequence:
                params["reference-base-sequence"] = reference_base_sequence
            if reference_collection:
                params["reference-collection"] = reference_collection
            if reference_location_type:
                params["reference-location-type"] = reference_location_type
            if reference_sequence:
                params["reference-sequence"] = reference_sequence
            if variant_base_sequence:
                params["variant-base-sequence"] = variant_base_sequence

            result = self._rest_adapter._do(
                http_method='GET', 
                endpoint=AlleleEndPoint.ALL.value,
                params=params
            )
            
            return AlleleCollection(data=result.data)