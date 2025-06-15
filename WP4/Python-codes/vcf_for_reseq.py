import sys
import re
from pathlib import Path

vcf_filename = sys.argv[1]
vcf_new_filename = f"{Path(vcf_filename).stem}_reseq.vcf"

chr10_contig = re.compile(r"##contig=<ID=chr10,length=133797422>")

with open(vcf_filename, mode='r') as vcf:
    with open(vcf_new_filename, mode = 'w') as new_vcf:
        for i, line in enumerate(vcf):
            # Duplicate first two lines
            if i < 2:
                new_vcf.write(line)

            # Remove any contig line, except for the chr10 contig
            elif line.startswith("##contig"):
                if chr10_contig.match(line):
                    new_vcf.write(line)

            # Write other meta-information lines to new vcf
            elif line.startswith("##"):
                new_vcf.write(line)

            # If line contains column names, add new columns
            elif line.startswith("#CHROM"):
                new_vcf.write(line.strip() + "\tFORMAT\tunknown\n")

            # If the line contains a variant definition, copy it and add values in the new columns
            else:
                new_vcf.write(line.strip() + "\tGT\t1\n")