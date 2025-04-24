from client import PharmVarApi
import pandas as pd
from yaml import safe_load
import numpy as np
import regex as re
import os

# Set path to working directory
os.chdir("/mnt/c/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP2/pharmvar_api")

# Get all core alleles in PharmVar database
pharmvar_api = PharmVarApi()
alleles = pharmvar_api.get_all_alleles(exclude_sub_alleles = True)

# Get all core and sub alleles in PharmVar database
subAlleles = pharmvar_api.get_all_alleles(exclude_sub_alleles = False)

# Define list of PharmVar pharmacogenes to investigate
pgx_genes = ["CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4", 
             "CYP3A5", "CYP4F2", "SLCO1B1", "NUDT15", "CYP2A13", "NAT2", "CYP1A2"]

# To be added later: DPYD

# Obtain total number of star alleles from PharmVar for all pharmacogenes:
total_SA = [0] * len(pgx_genes)
for i, gene in enumerate(pgx_genes):
    all_alleles = alleles.filter(gene_symbol = gene)
    total_SA[i] = all_alleles.__len__()

total_SA_dict = dict(zip(pgx_genes, total_SA))
print(total_SA_dict)

# Define most general star allele regular expression
star_allele_re = re.compile(r'\*\d+[A-Z]*')

# Define sub allele regular expression
sub_allele_re = re.compile(r'\*\d+\.\d+')

# Define legacy allele regular expression
legacy_allele_re = re.compile(r'\*\d+[A-Z]')

# Define non-legacy allele regular expression
nonlegacy_allele_re = re.compile(r'\*\d+')

# Define function to obtain names of all alleles for any gene in PharmVar
def getAlleleNames(gene: str) -> list:
    """
        Given the name of a pharmacogene in PharmVar, this function returns 
        a list containing the names of all core star alleles in the format "gene*X"
        with X the number of the star allele.
    """
    
    # Get all Allele objects from the gene
    gene_alleles = alleles.filter(gene_symbol = gene)
    
    # Extract allele names and return them in a list
    allele_names = [0] * gene_alleles.__len__()
    for i, allele in enumerate(gene_alleles):
        allele_names[i] = allele.allele_name
    return allele_names

# Define function to obtain names of all alleles for any gene in PharmVar
def getAlleleNamesShort(gene: str) -> list:
    """
        Given the name of a pharmacogene in PharmVar, this function returns 
        a list containing the names of all core star alleles in the format "*X"
        with X the number of the star allele.
    """
    
    # Get all Allele objects from the gene
    gene_alleles = alleles.filter(gene_symbol = gene)
    
    # Extract star alleles from allele names and return them in a list
    allele_names = [0] * gene_alleles.__len__()
    for i, allele in enumerate(gene_alleles):
        full_name = allele.allele_name
        allele_names[i] = nonlegacy_allele_re.findall(full_name)[0]
    return allele_names

# Define function to obtain names of all alleles for any gene in PharmVar
def getSubAlleleNames(gene: str) -> list:
    """
        Given the name of a pharmacogene in PharmVar, this function returns 
        a list containing the names of all star alleles in the format "gene*X"
        with X the number of the star allele, including all sub alleles.
    """
      
    # Get all Allele objects from the gene
    gene_alleles = subAlleles.filter(gene_symbol = gene)
    
    # Extract allele names and return them in a list
    allele_names = [0] * gene_alleles.__len__()
    for i, allele in enumerate(gene_alleles):
        allele_names[i] = allele.allele_name
    return allele_names

# Legacy alleles and legacy allele mappings for some genes
CYP1A1_legacy_alleles = ["2A", "2B", "2C"]
CYP1A1_legacy_allele_map = {"2A": "2", "2B": "2", "2C": "2"}

SLC01B1_legacy_alleles= ["1A", "1B", "17", "18", "20", "21", "22", "35"]
SLCO1B1_legacy_allele_map = {"1A": "1", "1B": "37", "17": "15", "18": "14", "20": "20", "21": "20", "22": "19", "35": "20"}

CYP1A2_legacy_alleles = ["1A", "1B", "1F", "1L", "1M", "1N", "1P", "1Q", "1R", "1S", "1T", "1U"]
CYP1A2_legacy_allele_map = {"1A": "1", "1B": "1", "1F": "30", "1L": "30", "1M": "30", "1N": "30", "1P": "30", "1Q": "30", "1R": "30", "1S": "1", "1T": "1", "1U": "1"}

CYP2A6_legacy_alleles = ["4D"]
CYP2A6_legacy_allele_map = {"4D": "47"}

CYP2C19_legacy_alleles = ["27", "1A"]
CYP2C19_legacy_allele_map = {"27": "1", "1A": "38"}

CYP2D6_legacy_alleles = ["57", "14A"]
CYP2D6_legacy_allele_map = {"57": "36", "14A": "114"}

CYP3A4_legacy_alleles = ["1A", "1B"]
CYP3A4_legacy_allele_map = {"1A": "1", "1B": "1"}

CYP3A5_legacy_alleles = ["2","3A", "3B", "3D", "3F", "3G", "3J", "3K", "3L", "4", "5", "10", "11"]
CYP3A5_legacy_allele_map = {"2": "3", "3A": "3", "3B": "3", "3D": "3", "3F": "3", "3G": "3", "3J": "3",
                     "3K": "3", "3L": "3", "4": "3", "5": "3", "10": "3", "11": "3"}

NAT2_legacy_alleles = ["12A", "12B", "12C", "12I", "11A", "13A", "20", "5C", "5B", "5F", "5W", "5BB", "6B", "6A", 
                       "6D", "6E", "6L", "7A", "7B", "14A", "14B", "14D", "5D", "5A", "14F", "14C", "5E", "5L", "5O", 
                       "5N", "6C", "6M", "6P", "6K", "6H", "6O", "7C", "12E", "12F", "12G", "12H", "12J", "14E", "14G", 
                       "14I", "14H", "12D", "12N", "6G"]
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

# For any other legacy allele of the form *(number)(LETTER), the allele is simply *(number)

def isLegacyAllele(gene: str, legacy_allele: str) -> str:
    """
        Given the name of a (pharmaco)gene and the formula of a star allele of that gene,
        this function will return True if that star allele is a legacy allele and False
        otherwise.
    """

    # Check if allele is a known legacy allele for the gene
    if gene in legacy_allele_dict and legacy_allele in legacy_allele_dict[gene]:
        return True

    # Check if allele is another legacy allele of the form *(number)(letter)
    if legacy_allele_re.match(f"*{legacy_allele}"):
        return True

    # If the allele is not a legacy allele, return False
    return False

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

## Initiate dictionary to store all callable star alleles per gene and per caller:
callable_alleles = {}

#############
### PyPGx ###
#############

# Read in star allele data
data_pypgx = pd.read_csv('Datasets/allele-table_pypgx.csv')

# Get all genes from pgx_genes for which at least one star allele can be called by PyPGx
genes_pypgx = [x for x in pgx_genes if x in data_pypgx["Gene"].unique()]

# For each of the selected genes, derive the number of star alleles that can be called by PyPGx
num_SA_pypgx = [0] * len(genes_pypgx)
for i, gene in enumerate(genes_pypgx):
    # Get part of dataframe that contains the current gene
    data_gene = data_pypgx.loc[data_pypgx['Gene'] == gene]

    # Extract unique star alleles
    star_alleles = data_gene['StarAllele'].unique().tolist()

    # Verify that star alleles have correct syntax
    star_alleles_correct = [x for x in star_alleles if star_allele_re.match(x)]

    # Make sure that legacy alleles are replaced by the correct allele
    star_alleles_correct2 = [mapLegacyAllele(gene, x[1:]) for x in star_alleles_correct]
    
    # Verify that star alleles are recognized by PharmVar
    valid_alleles = getAlleleNames(gene)
    star_alleles_correct3 = [x for x in star_alleles_correct2 if f"{gene}{x}" in valid_alleles]

    # Only keep unique star alleles
    star_alleles_unique = list(set(star_alleles_correct3))

    # Count number of unique verified star alleles
    num_SA_pypgx[i] = len(star_alleles_unique)
    
    # Save list of unique callable core star alleles
    callable_alleles[gene] = {}
    callable_alleles[gene]["PyPGx"] = star_alleles_unique

# Get total number of star alleles from PharmVar for all pharmacogenes in PyPGx
total_SA_pypgx = [total_SA_dict[x] for x in genes_pypgx if x in total_SA_dict]

# Get allele coverages for all genes
coverage_pypgx = np.array(num_SA_pypgx) / np.array(total_SA_pypgx)
results_pypgx = pd.DataFrame({"Gene": genes_pypgx, "PyPGx_Coverage": coverage_pypgx})
results = results_pypgx

## PROBLEEM: binnen de PyPGx data set staan meer star alleles geregistreerd voor CYP3A5 dan in PharmVar!!
## ook voor CYP2D6 en NAT2 is de coverage hoger als er niet gecheckt wordt of de star alleles ook in PharmVar staan!

## CONCLUSIE: voor CYP3A5, CYP2D6 en NAT2 kunnen bepaalde star alleles gecalled worden door PyPGx die geen deel uitmaken van PharmVar

##################
### StellarPGx ###
##################

# Genes that are not called by StellarPGx: DPYD, CYP2A13
genes_stellarPGx = [x for x in pgx_genes if x not in ["DPYD", "CYP2A13"]]

# For each gene, find number of star alleles that can be called by StellarPGx
num_SA_stellarPGx = [0] * len(genes_stellarPGx)
for i, gene in enumerate(genes_stellarPGx):

    # Read in data from StellarPGx
    data_gene = pd.read_csv(f"Datasets/stellarPGx_{gene}.txt", sep='\t', header=None)

    # Extract unique (core) star alleles from StellarPGx data
    star_alleles = list(set(data_gene[0]))

    # Verify that star alleles have correct syntax
    star_alleles_correct = [x for x in star_alleles if star_allele_re.match(x)]

    # Make sure that legacy alleles are replaced by the correct allele
    star_alleles_correct2 = [mapLegacyAllele(gene, x[1:]) for x in star_alleles_correct]
    
    # Verify that star alleles are recognized by PharmVar
    valid_alleles = getAlleleNames(gene)
    star_alleles_correct3 = [x for x in star_alleles_correct2 if f"{gene}{x}" in valid_alleles]

    # Only keep unique star alleles
    star_alleles_unique = list(set(star_alleles_correct3))

    # Count number of unique verified star alleles
    num_SA_stellarPGx[i] = len(star_alleles_unique)

    # Save list of unique callable core star alleles
    callable_alleles[gene]["StellarPGx"] = star_alleles_unique

# Get total number of star alleles from PharmVar for all pharmacogenes in StellarPGx
total_SA_stellarPGx = [total_SA_dict[x] for x in genes_stellarPGx if x in total_SA_dict]

# Get allele coverages for all genes
coverage_stellarPGx = np.array(num_SA_stellarPGx) / np.array(total_SA_stellarPGx)
results_stellarPGx = pd.DataFrame({"Gene": genes_stellarPGx, "StellarPGx_Coverage": coverage_stellarPGx})
results = results.merge(results_stellarPGx, how = "outer")

# Genes for which star alleles can be called by StellarPGx that are not registered in PharmVar:
# CYP2D6 (78% -> 76.8%), CYP3A4 (75.6% -> 73.3%), CYP3A5 (150% -> 100%)

# Some star alleles were included multiple times in the list of StellarPGx (different versions) =>
# ensure uniqueness of the star alleles!

################
### PharmCAT ###
################

# Get data on star alleles from pharmCAT:
data_pharmcat = pd.read_csv("Datasets/pharmCAT_alleles.tsv", sep='\t')

# Get names of all PharmVar genes that can be called by pharmCAT:
genes_pharmcat = [x for x in pgx_genes if x in set(data_pharmcat['Gene'])]

# Compute allele coverage for all genes in genes_pharmcat
num_SA_pharmcat = [0] * len(genes_pharmcat)
for i, gene in enumerate(genes_pharmcat):

    # Get data from pharmCAT for specific gene
    data_gene = data_pharmcat.loc[data_pharmcat['Gene'] == gene]    

    # Extract unique star alleles
    star_alleles = data_gene["Named Alleles"].to_list()[0]
    star_alleles_correct = star_allele_re.findall(star_alleles)
    star_alleles_correct = list(set(star_alleles_correct))

    # Make sure that legacy alleles are replaced by the correct allele
    star_alleles_correct2 = [mapLegacyAllele(gene, x[1:]) for x in star_alleles_correct]
    
    # Verify that star alleles are recognized by PharmVar
    valid_alleles = getAlleleNames(gene)
    star_alleles_correct3 = [x for x in star_alleles_correct2 if f"{gene}{x}" in valid_alleles]

    # Only keep unique star alleles
    star_alleles_unique = list(set(star_alleles_correct3))

    # Count number of unique verified star alleles
    num_SA_pharmcat[i] = len(star_alleles_unique)

    # Save list of unique callable core star alleles
    callable_alleles[gene]["PharmCAT"] = star_alleles_unique

# Get total number of star alleles from PharmVar for all pharmacogenes in PharmCAT
total_SA_pharmcat = [total_SA_dict[x] for x in genes_pharmcat if x in total_SA_dict]

# Get allele coverages for all genes
coverage_pharmcat = np.array(num_SA_pharmcat) / np.array(total_SA_pharmcat)
results_pharmcat = pd.DataFrame({"Gene": genes_pharmcat, "PharmCAT_Coverage": coverage_pharmcat})
results = results.merge(results_pharmcat, how = "outer")

###############
### ursaPGx ###
###############

# Note: the *1 allele is not registered among the callable alleles of ursaPGx, but it is always callable!
#       => Add 1 to number of callable alleles!

# SPECIAL CASE: CYP2D6 => called by Cyrius

# Get data on star alleles from ursaPGx:
data_ursapgx = pd.read_csv("Datasets/callable-vs-defined-alleles.tsv", sep='\t')

# Get names of all PharmVar genes that can be called by ursaPGx:
genes_ursapgx = [x for x in pgx_genes if x in set(data_ursapgx['gene'])]
genes_ursapgx.append("CYP2D6")

# Compute allele coverage for all genes in genes_ursapgx
num_SA_ursapgx = [0] * len(genes_ursapgx)
for i, gene in enumerate(genes_ursapgx):

    # Treat CYP2D6 separately
    if gene == "CYP2D6":
        
        # Read in data from Cyrius (sep = ',' ensures that a single string is returned per row)
        data_gene = pd.read_csv("Datasets/star_table_Cyrius.txt", sep=",", header=None)

        # Extract star alleles from Cyrius data
        star_alleles = data_gene[0]
        star_alleles_correct = [star_allele_re.findall(x)[0] for x in star_alleles if star_allele_re.match(x)]

    else:

        # Get data from ursaPGx for specific gene
        data_gene = data_ursapgx.loc[data_ursapgx['gene'] == gene]    

        # Extract star alleles from urszPGx data
        star_alleles = data_gene["callable"].to_list()[0]
        star_alleles_correct = star_allele_re.findall(star_alleles)

        # Add *1 allele
        star_alleles_correct.append('*1')    

    # Only keep unique star alleles
    star_alleles_correct = list(set(star_alleles_correct))

    # Make sure that legacy alleles are replaced by the correct allele
    star_alleles_correct2 = [mapLegacyAllele(gene, x[1:]) for x in star_alleles_correct]
    
    # Verify that star alleles are recognized by PharmVar
    valid_alleles = getAlleleNames(gene)
    star_alleles_correct3 = [x for x in star_alleles_correct2 if f"{gene}{x}" in valid_alleles]

    # Only keep unique star alleles
    star_alleles_unique = list(set(star_alleles_correct3))

    if gene == "CYP2D6":
        print(star_alleles_unique)

    # Count number of unique verified star alleles 
    num_SA_ursapgx[i] = len(star_alleles_unique)

    # Save list of unique callable core star alleles
    callable_alleles[gene]["ursaPGx"] = star_alleles_unique

# Get total number of star alleles from PharmVar for all pharmacogenes in ursaPGx
total_SA_ursapgx = [total_SA_dict[x] for x in genes_ursapgx if x in total_SA_dict]

# Get allele coverages for all genes
coverage_ursapgx = np.array(num_SA_ursapgx) / np.array(total_SA_ursapgx)
results_ursapgx = pd.DataFrame({"Gene": genes_ursapgx, "ursaPGx_Coverage": coverage_ursapgx})
results = results.merge(results_ursapgx, how = "outer")

#############
### PAnno ###
#############

# Read in star allele information from PAnno
data_panno = pd.read_json('Datasets/PAnno_pgx_diplotypes.json')

# Get names of all PharmVar genes that can be called by PAnno:
genes_panno = [x for x in pgx_genes if x in data_panno.columns]

# Compute allele coverage for all genes in genes_ursapgx
num_SA_panno = [0] * len(genes_panno)
for i, gene in enumerate(genes_panno):
    
    # Obtain all star alleles of gene that can be called by PAnno
    data_gene = data_panno[gene].loc['haplotype_mutated_loci']
    star_alleles = list(data_gene.keys())

    # Verify that star alleles have correct syntax
    star_alleles_correct = [x for x in star_alleles if star_allele_re.match(x)]

    # Make sure that legacy alleles are replaced by the correct allele
    star_alleles_correct2 = [mapLegacyAllele(gene, x[1:]) for x in star_alleles_correct]
    
    # Verify that star alleles are recognized by PharmVar
    valid_alleles = getAlleleNames(gene)
    star_alleles_correct3 = [x for x in star_alleles_correct2 if f"{gene}{x}" in valid_alleles]

    # Only keep unique star alleles
    star_alleles_unique = list(set(star_alleles_correct3))

    # Count number of unique verified star alleles
    num_SA_panno[i] = len(star_alleles_unique)

    # Save list of unique callable core star alleles
    callable_alleles[gene]["PAnno"] = star_alleles_unique

# Get total number of star alleles from PharmVar for all pharmacogenes in PAnno
total_SA_panno = [total_SA_dict[x] for x in genes_panno if x in total_SA_dict]

# Get allele coverages for all genes
coverage_panno = np.array(num_SA_panno) / np.array(total_SA_panno)
results_panno = pd.DataFrame({"Gene": genes_panno, "PAnno_Coverage": coverage_panno})
results = results.merge(results_panno, how = "outer")

############
### Aldy ###
############

# DPYD: gebruikt nog sterallel systeem!

# Get all Aldy genes
genes_aldy = pgx_genes

# For each gene, find number of star alleles that can be called by Aldy
num_SA_aldy = [0] * len(genes_aldy)
for i, gene in enumerate(genes_aldy):

    # Get data from Aldy for gene
    with open(f'Datasets/{gene.lower()}.yml', 'r') as f:
        loaded_data = safe_load(f)
    
    # Extract all star alleles
    allele_data = loaded_data["alleles"]
    star_alleles_full = allele_data.keys()

    # Differentiate between alleles that are labeled and those that are not
    star_alleles_label = [allele_data[x]["label"] for x in star_alleles_full if "label" in allele_data[x].keys()]
    star_alleles_no_label = [x for x in star_alleles_full if "label" not in allele_data[x].keys()]

    # Find star alleles in core/legacy allele format (e.g. *1 or *1B)
    star_alleles_label = [star_allele_re.findall(x)[0] for x in star_alleles_label if star_allele_re.search(x)]

    # Find star alleles in the sub allele format (e.g. *1.002)
    star_alleles_no_label = [sub_allele_re.findall(x)[0] for x in star_alleles_no_label if sub_allele_re.search(x)]

    # Verify that sub star alleles are recognized by PharmVar
    valid_sub_alleles = getSubAlleleNames(gene)
    sub_alleles_correct = [x for x in star_alleles_no_label if f"{gene}{x}" in valid_sub_alleles]

    # Convert sub alleles to core alleles (i.e. *1; *2 instead of *1.001; *2.003)
    sub_alleles_correct2 = list(set([nonlegacy_allele_re.findall(x)[0] for x in sub_alleles_correct if nonlegacy_allele_re.search(x)]))

    # Combine converted sub alleles with core/legacy alleles
    star_alleles_correct = list(set(star_alleles_label + sub_alleles_correct2))

    # Convert legacy alleles (if any)
    star_alleles_correct2 = [mapLegacyAllele(gene, x[1:]) for x in star_alleles_correct]

     # Verify that star alleles are recognized by PharmVar
    valid_alleles = getAlleleNames(gene)
    star_alleles_correct3 = [x for x in star_alleles_correct2 if f"{gene}{x}" in valid_alleles]

    # Get unique star alleles
    star_alleles_unique = list(set(star_alleles_correct3))

    # Count number of unique verified star alleles
    num_SA_aldy[i] = len(star_alleles_unique)

    # Save list of unique callable core star alleles
    callable_alleles[gene]["Aldy"] = star_alleles_unique

# Get total number of star alleles from PharmVar for all pharmacogenes in Aldy
total_SA_aldy = [total_SA_dict[x] for x in genes_aldy if x in total_SA_dict]

# Get allele coverages for all genes
coverage_aldy = np.array(num_SA_aldy) / np.array(total_SA_aldy)
results_aldy = pd.DataFrame({"Gene": genes_aldy, "Aldy_Coverage": coverage_aldy})
results = results.merge(results_aldy, how = "outer")

## CYP2D6 coverage ging omhoog na toevoeging van sub allele search => nog uitzoeken waarom dit zo is!

### Print final results ###
print(results)

# Export results to .csv file
#results.to_csv("DiplotypeCallerAlleleCoverage2.csv", sep="\t", index=False, header=True)

## Legacy alleles omzetten heeft voorlopig enkel een impact of CYP2D6 voor StellarPGx:
## StellarPGx: het legacy allel *57 is aanwezig, maar het sterallel *36 (waarnaar *57 mapt) niet => hogere coverage

#########################
### COMBINING CALLERS ###
#########################

# Get list of callers
callers = ["PyPGx", "StellarPGx", "PharmCAT", "ursaPGx", "PAnno", "Aldy"]

# Complete callable_alleles dictionary by adding empty lists when a gene is not included
callable_alleles["CYP2A13"].update({"StellarPGx": [], "PharmCAT": [], "PAnno": []})
callable_alleles["CYP2A6"].update({"PharmCAT": [], "PAnno": []})
callable_alleles["CYP2C8"].update({"PharmCAT": []})
callable_alleles["NAT2"].update({"PharmCAT": [], "ursaPGx": [], "PAnno": []})
callable_alleles["CYP1A2"].update({"PharmCAT": [], "ursaPGx": [], "PAnno": []})

# Convert callable_alleles dictionary to data frame, to be exported
callable_alleles_df = pd.DataFrame({"Caller": callers})
for i, gene in enumerate(pgx_genes):
    
    callable_alleles_gene = callable_alleles[gene]
    callable_alleles_list = [callable_alleles_gene[x] for x in callers]
    
    temp_df = pd.DataFrame({"Caller": callers, f"{gene}": callable_alleles_list})
    callable_alleles_df = callable_alleles_df.merge(temp_df, how = "outer")
    
print(callable_alleles_df)
callable_alleles_df.to_csv("Callable_alleles.csv", sep = "\t", index = False, header = True)

# For each gene, get list of star alleles recognized by PharmVar and export them
defined_alleles_df = pd.DataFrame({"Gene": pgx_genes})
defined_alleles = []
for i, gene in enumerate(pgx_genes):
    
    # Extract names of valid star alleles
    valid_alleles = getAlleleNamesShort(gene)
    defined_alleles.append(valid_alleles)

defined_alleles_df["Alleles"] = defined_alleles
#defined_alleles_df.to_csv("Defined_alleles.csv", sep = "\t", index = False, header = True)

# First, check for each of these genes whether all star alleles can be called by at least one of the callers
# Get list of genes that don't have 100% coverage for any star allele caller
genes_low_coverage = ["CYP2A6", "CYP2B6", "CYP2D6", "CYP4F2", "NAT2"]

is_covered = [False] * len(genes_low_coverage)
not_covered = {}
for i, gene in enumerate(genes_low_coverage):
    
    # Extract lists of callable star alleles and convert to sets
    callable_alleles_gene = list(callable_alleles[gene].values())

    # Obtain union of all lists
    callable_alleles_gene_all = set().union(*callable_alleles_gene)

    # Get total number of core star alleles registered in PharmVar
    total_alleles = total_SA_dict[gene]

    # Check if same number of alleles can be called by any combination of callers
    if len(callable_alleles_gene_all) == total_alleles: 
        is_covered[i] = True
        not_covered[gene] = []

    # Otherwise, find out which star alleles are not covered
    else:
        all_alleles = getAlleleNamesShort(gene)
        alleles_not_found = set(all_alleles).difference(callable_alleles_gene_all)
        not_covered[gene] = list(alleles_not_found)

print(is_covered)
print(not_covered)

# Conclusion: only CYP2B6 is covered by a combination of callers