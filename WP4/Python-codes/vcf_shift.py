import sys
import re
from pathlib import Path

vcf_filename = sys.argv[1]
#new_start_pos = int(sys.argv[2])
new_start_pos = 94682000

vcf_new_filename = f"{Path(vcf_filename).stem}_shift.vcf"

# Define variant line regular expression
variant_pattern = re.compile(r'(\S+\t)(\d+)(\t.+\n)')

with open(vcf_filename, mode='r') as vcf:
    with open(vcf_new_filename, mode = 'w') as new_vcf:
        for line in vcf:
            # If the line starts with '#', just write it into the new file
            if line.startswith("#"):
                new_vcf.write(line)
            # Otherwise, adjust the position on the line
            else:
                before, pos, after = variant_pattern.findall(line)[0]
                pos = int(pos)
                new_pos = pos - new_start_pos + 1
                new_line = f"{before}{new_pos}{after}"
                new_vcf.write(new_line)
                