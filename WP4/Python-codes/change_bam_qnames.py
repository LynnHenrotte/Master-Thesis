from pysam import AlignmentFile
from pathlib import Path
import sys

bam_in_filename = sys.argv[1] # First argument: BAM file to be modified (haplotype!!)
bam_out_filename = f"{Path(bam_in_filename).stem}_mod.bam"
modification = sys.argv[2] # This string will be added at the end of each QNAME in the input BAM file

bam_in = AlignmentFile(bam_in_filename)
bam_out = AlignmentFile(bam_out_filename, mode = "w", template = bam_in) # template makes sure to copy the header

for read in bam_in.fetch(until_eof = True):
    read.query_name += modification
    bam_out.write(read)
bam_in.close()