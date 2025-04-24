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

# Obtain total number of star alleles from PharmVar for all pharmacogenes:
total_SA = [0] * len(pgx_genes)
for i, gene in enumerate(pgx_genes):
    all_alleles = subAlleles.filter(gene_symbol = gene)
    total_SA[i] = all_alleles.__len__()

total_SA_dict = dict(zip(pgx_genes, total_SA))
print(total_SA_dict)