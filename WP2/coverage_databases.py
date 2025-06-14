import sys
sys.path.insert(1, 'pharmvar_api')

from client import PharmVarApi
import pandas as pd
import regex as re
import os

# Get all alleles in PharmVar database
pharmvar_api = PharmVarApi()
alleles = pharmvar_api.get_all_alleles(exclude_sub_alleles = True)

# Regular expression definitions #
star_allele_re = re.compile(r'\*\d+[A-Z]*') # most general star allele
sub_allele_re = re.compile(r'\*\d+\.\d+') # sub allele
legacy_allele_re = re.compile(r'\*\d+[A-Z]') # legacy allele
nonlegacy_allele_re = re.compile(r'\*\d+') # non-legacy allele

# Define function to obtain names of all alleles for any gene in PharmVar #
def getAlleleNames(gene: str, short = False) -> list:
    """
        Given the name of a pharmacogene in PharmVar, this function returns 
        a list containing the names of all core star alleles in the format "gene*X",
        with X the number of the star allele.
    """
    
    # Get all Allele objects from the gene
    gene_alleles = alleles.filter(gene_symbol = gene)
    
    # Extract allele names and return them in a list
    allele_names = [0] * gene_alleles.__len__()
    for i, allele in enumerate(gene_alleles):
        full_name = allele.allele_name
        if short:
            allele_names[i] = nonlegacy_allele_re.findall(full_name)[0]
        else:
            allele_names[i] = full_name
    return allele_names

# Legacy allele mappings for some genes
CYP1A1_legacy_allele_map = {"2A": "2", "2B": "2", "2C": "2"}
SLCO1B1_legacy_allele_map = {"1A": "1", "1B": "37", "17": "15", "18": "14", "20": "20", "21": "20", "22": "19", "35": "20"}
CYP1A2_legacy_allele_map = {"1A": "1", "1B": "1", "1F": "30", "1L": "30", "1M": "30", "1N": "30", "1P": "30", "1Q": "30", 
                            "1R": "30", "1S": "1", "1T": "1", "1U": "1"}
CYP2A6_legacy_allele_map = {"4D": "47"}
CYP2C19_legacy_allele_map = {"27": "1", "1A": "38"}
CYP2D6_legacy_allele_map = {"57": "36", "14A": "114"}
CYP3A4_legacy_allele_map = {"1A": "1", "1B": "1"}
CYP3A5_legacy_allele_map = {"2": "3", "3A": "3", "3B": "3", "3D": "3", "3F": "3", "3G": "3", "3J": "3",
                     "3K": "3", "3L": "3", "4": "3", "5": "3", "10": "3", "11": "3"}
NAT2_legacy_allele_map = {"12A" : "1", "12B" : "1", "12C" : "1", "12I" : "1", "11A" : "4", "13A" : "4", "20" : "4",
                          "5C" : "5", "5B" : "5", "5F" : "5", "5W" : "5", "5BB" : "5", "6B" : "6", "6A" : "6", "6D" : "6",
                          "6E" : "6", "6L" : "6", "7A" : "7", "7B" : "7", "14A" : "14", "14B" : "14", "14D" : "15",
                          "5D" : "16", "5A" : "16", "14F" : "29", "14C" : "29", "5E" : "30", "5L" : "31", "5O" : "32",
                          "5N" : "33", "6C" : "34", "6M" : "35",  "6P" : "36", "6K" : "37", "6H" : "38", "6O" : "39", 
                          "7C" : "40", "12E" : "41", "12F" : "42", "12G" : "43", "12H" : "44", "12J" : "45", "14E" : "46",
                          "14G" : "46", "14I" : "46", "14H" : "47", "12D" : "48",  "12N" : "51", "6G" : "52"}

# Define dictionary containing the above dictionaries
legacy_allele_dict = {"CYP1A1": CYP1A1_legacy_allele_map, "SLC01B1": SLCO1B1_legacy_allele_map,
                      "CYP1A2": CYP1A2_legacy_allele_map, "CYP2A6": CYP2A6_legacy_allele_map,
                      "CYP2C19": CYP2C19_legacy_allele_map, "CYP2D6": CYP2D6_legacy_allele_map,
                      "CYP3A4": CYP3A4_legacy_allele_map, "CYP3A5": CYP3A5_legacy_allele_map,
                      "NAT2": NAT2_legacy_allele_map}

# Define function that maps known legacy alleles to current alleles
def mapLegacyAllele(gene: str, legacy_allele: str) -> str:
    """
        Given the name of a (pharmaco)gene and the formula of a star allele of that gene,
        this function will return the current star allele that the given allele corresponds to
        (if it is a legacy allele). Otherwise, the inputted star allele is returned
    """

    # Check if gene has predefined legacy alleles   
    if gene in legacy_allele_dict:
        legacy_alleles = legacy_allele_dict[gene]

        # Check if given allele is a predefined legacy allele for the gene
        if legacy_allele in legacy_alleles:
            mapping = legacy_alleles[legacy_allele]
            return f"*{mapping}"
    
    # Check if allele is another (basic) legacy allele
    if legacy_allele_re.match(f"*{legacy_allele}"):
        mapping = nonlegacy_allele_re.findall(f"*{legacy_allele}")[0]
        return mapping
        
    return f"*{legacy_allele}"

# Define function to obtain allele coverage
def get_allele_coverage(gene: str, database: pd.DataFrame) -> list[float, list[str]]:
    """
    Given the name of a pharmacogene, this function returns a list containing
        1) its star allele coverage in a given data set, and
        2) a list containing the star alleles of the gene in the given database
    """
    # Extract gene data from given data set
    gene_data = database[gene]
    gene_data = gene_data.to_list()

    # Extract all star alleles that are encountered
    star_allele_list = []
    for x in range(len(gene_data)):
        entry = gene_data[x]

        # Check if entry is NaN:
        if entry != entry:
            continue

        # Otherwise, extract all star alleles from the diplotype
        else:
            all_alleles = star_allele_re.findall(gene_data[x])
            star_allele_list.extend(all_alleles)

    # Convert list to set to get unique star alleles
    star_alleles = list(set(star_allele_list))

    # Make sure that legacy alleles are replaced by the correct allele and only keep unique core alleles
    star_alleles_correct = list(set([mapLegacyAllele(gene, x[1:]) for x in star_alleles]))
    
    # Verify that star alleles are recognized by PharmVar
    valid_alleles = getAlleleNames(gene)
    star_alleles_correct2 = [x for x in star_alleles_correct if f"{gene}{x}" in valid_alleles]

    # Get number of unique star alleles found
    num_SA_database = len(star_alleles_correct2)

    # Get total number of star alleles as defined in PharmVar
    test_gene_all_alleles = alleles.filter(gene_symbol = gene)
    num_SA_all = test_gene_all_alleles.__len__()

    # Get percentage of star alleles found in given data set
    coverage = num_SA_database / num_SA_all
    return [coverage, star_alleles_correct2]

# Define list for all pharmacogenes available in GeT-RM and Star Allele Search
pharmacogenes_getRM = ["CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6", 
                       "CYP3A4", "CYP3A5", "CYP4F2", "SLCO1B1", "NAT2", "CYP1A2"]

pharmacogenes_starAlleleSearch = ["CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6", 
                                  "CYP3A4", "CYP3A5", "CYP4F2", "SLCO1B1", "NUDT15", "CYP2A13"]

genes_together = list(set(pharmacogenes_getRM).union(set(pharmacogenes_starAlleleSearch)))

# Read in data sets
data_getRM = pd.read_excel("Datasets/updated_getrm_calls2.xlsx")
data_starAlleleSearch = pd.read_excel("Datasets/Star-allele-search-db.xlsx")

# Initiate list to store star alleles that are encountered in both datasets
covered_alleles = {}

# Get allele coverages for all PharmVar pharmacogenes that are present in GeT-RM
coverages_getRM = [0] * len(pharmacogenes_getRM)
for i, gene in enumerate(pharmacogenes_getRM):

    coverage = get_allele_coverage(pharmacogenes_getRM[i], data_getRM)
    
     # First obtain percentage coverage of the gene in GeT-RM
    coverages_getRM[i] = coverage[0]

    # Initiate dictionary and add list of star alleles covered in GeT-RM
    covered_alleles[gene] = {}
    covered_alleles[gene]["GeT-RM"] = coverage[1]

# Get allele coverages for all PharmVar pharmacogenes that are present in Star Allele Search
coverages_starAlleleSearch = [0] * len(pharmacogenes_starAlleleSearch)
for i, gene in enumerate(pharmacogenes_starAlleleSearch):

    coverage = get_allele_coverage(pharmacogenes_starAlleleSearch[i], data_starAlleleSearch)

    # First obtain percentage coverage of the gene in Star Allele Search
    coverages_starAlleleSearch[i] = coverage[0]

    # Add list of star alleles covered in Star Allele Search to gene dictionary (if it exists)
    if gene not in covered_alleles.keys(): covered_alleles[gene] = {}
    covered_alleles[gene]["StarAlleleSearch"] = coverage[1]

results_getRM = pd.DataFrame({"Gene": pharmacogenes_getRM, "GeT-RM": coverages_getRM})
results_starAlleleSearch = pd.DataFrame({"Gene": pharmacogenes_starAlleleSearch, "Star allele search": coverages_starAlleleSearch})
result_all = results_getRM.merge(results_starAlleleSearch, how = "outer")
print(result_all)

# Fill up dictionary with non-covered genes
covered_alleles["NUDT15"].update({"GeT-RM": []})
covered_alleles["CYP2A13"].update({"GeT-RM": []})
covered_alleles["NAT2"].update({"StarAlleleSearch": []})
covered_alleles["CYP1A2"].update({"StarAlleleSearch": []})

# Convert dictionary to dataframe to be exported:
covered_alleles_df = pd.DataFrame({"Database": ["GeT-RM", "StarAlleleSearch"]})
for i, gene in enumerate(genes_together):

    # Extract callable alleles for gene
    gene_data = pd.DataFrame(covered_alleles[gene].items())
    gene_data.columns = ["Database", f"{gene}"]

    # Add to dataframe
    covered_alleles_df = covered_alleles_df.merge(gene_data)
    
covered_alleles_df.to_csv("Covered_Alleles_Databases.csv", sep = "\t", index = False, header = True)

# Export results to .csv file
result_all.to_csv("Database_Allele_Coverage.csv", sep = "\t", index = False, header = True)

# For each gene, get list of star alleles recognized by PharmVar and export them
defined_alleles_df = pd.DataFrame({"Gene": genes_together})
defined_alleles = []
for i, gene in enumerate(genes_together):
    
    # Extract names of valid star alleles
    valid_alleles = getAlleleNames(gene, short = True)
    defined_alleles.append(valid_alleles)

defined_alleles_df["Alleles"] = defined_alleles

# Export results to .csv file
defined_alleles_df.to_csv("Defined_alleles_all.csv", sep = "\t", index = False, header = True)
