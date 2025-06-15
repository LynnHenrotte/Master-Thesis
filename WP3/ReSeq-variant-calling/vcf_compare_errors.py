import sys
import os
import re
from pathlib import Path

vcf_filename = sys.argv[1]
errors_reverse = sys.argv[2]
errors_forward = sys.argv[3]

out_filename = f"{Path(vcf_filename).stem}.txt"
print(out_filename)
new_start_pos = 94682000 # for adjustment
adjust_start_pos = 1 # this is the number we add to pos - new_start_pos
seq_length = 389001 # for reverse sequence, can change due to indels

# Obtain error sequence from reverse strand
with open(errors_reverse, "r") as rev:
    rev_seq = rev.readlines()[1].strip("\n")
    
# Obtain error sequence from forward strand
with open(errors_forward, "r") as forw:
    for_seq = forw.readlines()[1].strip("\n")

# Write function that computes a nucleotide's complement
def comp(nucleotide: str) -> str:
    """ Takes as argument one of the four nucleotides A, C, T and G, or N. 
        If N is given, the function returns D. Otherwise, the complementary
        nucleotide is returned. 
    """
    if nucleotide not in ["A", "C", "G", "T", "N"]:
        raise Exception("Not a valid nucleotide! Choose one of: A, C, G, T or N.")
    if nucleotide == "N":
        return "D"
    if nucleotide == "A":
        return "T"
    if nucleotide == "T":
        return "A"
    if nucleotide == "C":
        return "G"
    if nucleotide == "G":
        return "C"

# Open vcf file and new text file for output
with open(vcf_filename, "r") as vcf, open(out_filename, "w") as out:

    # Write header for output file
    out.write(f"#POS\t#REF\t#ALT\t#QUAL\t#TYPE\t#ERROR\t#STRAND\n")

    for line in vcf:
        
        ## Skip over meta-information lines
        if line.startswith("#"):
            continue
        
        ## For each variant line,
        else:
            
            ### Obtain necessary information about variant
            entries = re.split(r'\t+', line.rstrip('\t'))
            pos = int(entries[1])
            ref = entries[3]
            alt = entries[4]
            qual = float(entries[5])
            #print(f"Position: {pos}, REF: {ref}, ALT: {alt}, quality: {qual}.") 

            ### Check if error is indel or substitution
            insertion = False
            deletion = False
            snp = False
            if len(ref) > 1:
                deletion = True
                out_message = f"{pos}\t{ref}\t{alt}\t{qual}\tDEL\t.\t.\n"
            elif len(alt) > 1:
                insertion = True
                out_message = f"{pos}\t{ref}\t{alt}\t{qual}\tINS\t.\t.\n"
            else:
                snp = True
                out_message = f"{pos}\t{ref}\t{alt}\t{qual}\tSNP\t"
            
            # Complete out message later
            end_out_message = f""

            ### If an indel was found, don't do anything yet
            if insertion or deletion:
                out.write(out_message)

            ### Scan forward fastq
            
            #### Find position, adjusted to CYP2C locus
            pos_adj = pos - new_start_pos + 1

            #### For now, skip over any variants that are outside of the CYP2C locus
            if pos_adj > len(for_seq) or pos_adj < 1:
                continue

            if snp == True:

                #### Extract nucleotide at position
                nucleotide = for_seq[pos_adj - 1] # -1 because Python is 0-based

                #### Check if nucleotide corresponds to ALT 
                if nucleotide == alt:
                    out_message += f"YES\tFORWARD\n"
                    out.write(out_message)
                    continue
                
                #### Check if some other nucleotide different from N is found 
                elif nucleotide != "N":
                    end_out_message = f"INCONCLUSIVE\tFORWARD\n"                    

            ### If not found in forward, scan reverse fastq

            #### Further adjust position, from end of fastq
            pos_adj = seq_length - pos_adj
            
            if snp == True:

                #### Extract nucleotide at position
                nucleotide = rev_seq[pos_adj]

                #### If nucleotide corresponds to complement of ALT, write YES to output
                if nucleotide == comp(alt):
                    out_message += f"YES\tREVERSE\n" 
                    out.write(out_message)
                    continue

                #### Check if another nucleotide different from N is found
                if nucleotide != "N":
                    out_message += f"INCONCLUSIVE\tREVERSE\n"

                # Check if an error was found on the forward strand
                elif end_out_message == f"INCONCLUSIVE\tFORWARD\n":
                    out_message += end_out_message

                #### Otherwise, write NO to output, since error was not found
                else:
                    out_message += f"NO\t.\n"
                
                out.write(out_message)
                
# Conclude by sorting the output file by decreasing quality
out_sorted_filename = f"{Path(vcf_filename).stem}_sorted.txt"

## Define function for sorting
def sort_out(line: str) -> int:
    entries = re.split(r'\t+', line.rstrip('\t'))
    qual = float(entries[3])
    return qual

with open(out_filename, "r") as out, open(out_sorted_filename, "w") as out_sort:
    
    # Obtain list of lines in current output file
    out_lines = out.readlines()

    # Store and remove header before sorting
    header = out_lines[0]
    out_lines = out_lines[1:]

    # Sort the remaining lines based on descending quality
    out_lines.sort(key = sort_out, reverse = True)

    # Write header to sorted output file
    out_sort.write(header)

    # Write sorted lines to sorted output file
    for line in out_lines:
        out_sort.write(line)

# Finally, remove the unsorted output file
os.remove(out_filename)