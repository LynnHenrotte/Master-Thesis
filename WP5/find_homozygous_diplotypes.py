import sys
import os
import re
from pathlib import Path

in_filename = sys.argv[1]

# Diplotype regex pattern
diplotype_pattern = re.compile(r"\*(\d+)/\*(\d+)")

# Set counter for homozygous diplotypes
count_homo = 0

# Scan file for homozygous diplotypes
with open(in_filename, "r") as file:
    for line in file:
        alleles = diplotype_pattern.search(line).groups()
        if alleles[0] == alleles[1]:
            count_homo += 1
            print(alleles)
    print(count_homo)