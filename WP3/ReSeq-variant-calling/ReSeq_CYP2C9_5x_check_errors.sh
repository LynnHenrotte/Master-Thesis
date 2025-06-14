#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-check-errors-CYP2C9-5x_%a # Name of job
#SBATCH --output=stdout-ReSeq-check-errors-CYP2C9-5x_%a   # Standard output name
#SBATCH --error=stderror-ReSeq-check-errors-CYP2C9-5x_%a    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=4     # Number of CPU cores
#SBATCH --time=3-00:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
#SBATCH --array=3-85        # First run this with array=2; when finished, run with array=3-85
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

# Set paths and parameters
ALLELE=$SLURM_ARRAY_TASK_ID
DATA=/data/leuven/373/vsc37367/
REF=${DATA}references/chr10.fa
OUT=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_5x/
VARCALL_VCF=${OUT}ReSeq_CYP2C9_${ALLELE}_VarCall.vcf
ALLELE_VCFS=${DATA}CYP2C9-6.2.3/GRCh38/
ALLELE_VCF=${ALLELE_VCFS}CYP2C9_${ALLELE}.vcf
ERRORS=${OUT}ReSeq_CYP2C9_${ALLELE}_errors.fq
PYTHON_FUNC=${DATA}Simulations_CYP2C9/Python-codes/vcf_compare_errors.py

source activate bcftools

## Prepare star allele VCF ##

if ! [ -f ${ALLELE_VCFS}CYP2C9_${ALLELE}_norm.vcf ]
then
    # Compress original vcf file
    bgzip -k -o $ALLELE_VCF.gz $ALLELE_VCF

    # Normalize vcf file
    ALLELE_VCF_NORM=${ALLELE_VCFS}CYP2C9_${ALLELE}_norm.vcf
    bcftools norm -f $REF -o $ALLELE_VCF_NORM $ALLELE_VCF.gz

    # Compress normalized vcf file
    bgzip -k -o $ALLELE_VCF_NORM.gz $ALLELE_VCF_NORM

    # Index vcf file
    bcftools index -o $ALLELE_VCF_NORM.gz.csi $ALLELE_VCF_NORM.gz

    # Remove non-normalized compressed and indexed vcf file
    rm -f $ALLELE_VCF.gz*
else 
    ALLELE_VCF_NORM=${ALLELE_VCFS}CYP2C9_${ALLELE}_norm.vcf
fi

## Prepare variant call VCF ##

# Compress original vcf file
bgzip -k -o $VARCALL_VCF.gz $VARCALL_VCF

# Remove called variants that are outside of the CYP2C locus in variant call vcf
bcftools index -o $VARCALL_VCF.gz.csi $VARCALL_VCF.gz
bcftools view -r chr10:94682000-95071000 -o $VARCALL_VCF $VARCALL_VCF.gz
bgzip -k -o $VARCALL_VCF.gz $VARCALL_VCF
rm -f $VARCALL_VCF.gz.csi

# Normalize vcf file
VARCALL_VCF_NORM=${OUT}ReSeq_CYP2C9_${ALLELE}_VarCall_norm.vcf
bcftools norm -f $REF -o $VARCALL_VCF_NORM $VARCALL_VCF.gz

# Compress normalized vcf file
bgzip -k -o $VARCALL_VCF_NORM.gz $VARCALL_VCF_NORM

# Index vcf file
bcftools index -o $VARCALL_VCF_NORM.gz.csi $VARCALL_VCF_NORM.gz

# Remove non-normalized compressed and indexed vcf file
rm -f $VARCALL_VCF.gz*

## Verify that the required variants have been simulated/called ##

# Find variants that are shared by both vcf files
bcftools isec -n =2 -p ${OUT}ReSeq_CYP2C9_${ALLELE}_compare_vcf $ALLELE_VCF_NORM.gz $VARCALL_VCF_NORM.gz

# Find called variants that are not part of the allele
bcftools isec -C -w 1 -o ${OUT}ReSeq_CYP2C9_${ALLELE}_VarCallErrors_norm.vcf $VARCALL_VCF_NORM.gz $ALLELE_VCF_NORM.gz
VARCALL_ERROR_VCF=${OUT}ReSeq_CYP2C9_${ALLELE}_VarCallErrors_norm.vcf

conda deactivate
    
## Verify that all remaining called variants are due to sequencing errors introduced by the simulation ##

# Split ReSeq error file   
source activate seqtk

seqtk split -n 2 ${OUT}ReSeq_CYP2C9_${ALLELE}_errors $ERRORS

conda deactivate

# Check if remaining called variants are due to errors
source activate python

ERRORS_REVERSE=${OUT}ReSeq_CYP2C9_${ALLELE}_errors.00001.fa
ERRORS_FORWARD=${OUT}ReSeq_CYP2C9_${ALLELE}_errors.00002.fa

python $PYTHON_FUNC $VARCALL_ERROR_VCF $ERRORS_REVERSE $ERRORS_FORWARD

conda deactivate

# Some file management
OUT_ERRORS=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_${ALLELE}_VarCallErrors_norm_sorted.txt
rm -f $ERRORS_FORWARD $ERRORS_REVERSE
mv $OUT_ERRORS ${OUT}
OUT_ERRORS=${OUT}ReSeq_CYP2C9_${ALLELE}_VarCallErrors_norm_sorted.txt

## Summarize information on allele ##

source activate gatk4

# Count variants of the allele
NUM_VAR="$( gatk CountVariants -V $ALLELE_VCF | sed -n 2p )"

# Count total variants called by HaplotypeCaller
TOT_CALLED="$( gatk CountVariants -V $VARCALL_VCF_NORM | sed -n 2p )"

# Count total non-allelic variants called by HaplotypeCaller
NON_ALLELE_CALLED="$( gatk CountVariants -V $VARCALL_ERROR_VCF | sed -n 2p )"

# Count number of allelic variants called by HaplotypeCaller
COMP_FILE=${OUT}ReSeq_CYP2C9_${ALLELE}_compare_vcf/sites.txt
ALLELE_CALLED=$( wc -l < "${COMP_FILE}" )

conda deactivate

# Compute binary indicator for allelic variants called
if [ $ALLELE_CALLED == $NUM_VAR ]
then 
    ALL_CALLED="YES"
else
    ALL_CALLED="NO"
fi

# Count number of called variants that were linked to simulated errors
ERRORS_FOUND="$( grep -o YES ${OUT_ERRORS} | wc -l )"

# Count number of called variants that were inconclusively linked to simulated errors
ERRORS_FOUND_INC="$( grep -o INCONCLUSIVE ${OUT_ERRORS} | wc -l )"

# Compute percentage of non-allelic called variants that are linked to confirmed simulated errors
ERROR_PERCENT="$( echo "$ERRORS_FOUND" "$NON_ALLELE_CALLED" | awk '{printf "%.1f%%", $1 / $2 * 100 }' )"

# Compute percentage of non-allelic called variants that are linked to simulated errors, confirmed or inconclusive
ERROR_PERCENT_INC="$( echo "$ERRORS_FOUND" "$ERRORS_FOUND_INC" "$NON_ALLELE_CALLED" | awk '{printf "%.1f%%", ($1 + $2) / $3 * 100 }' )"

# Compute average quality score of allelic variant calls
source activate bcftools

ALLELE_CALLS=${OUT}ReSeq_CYP2C9_${ALLELE}_compare_vcf/0001.vcf
bcftools query -f '%QUAL\n' $ALLELE_CALLS > qual_scores_temp1_${ALLELE}.txt

conda deactivate

source activate python

MEAN_QUAL_PY=${DATA}Simulations_CYP2C9/Python-codes/vcf_average_qual.py
MEAN_QUAL_ALLELIC="$( python ${MEAN_QUAL_PY} qual_scores_temp1_${ALLELE}.txt )"

# Compute minimal quality score of allelic variant calls

EXTREMUM_PY=${DATA}Simulations_CYP2C9/Python-codes/extremum.py
MIN_QUAL_ALLELIC="$( python ${EXTREMUM_PY} min qual_scores_temp1_${ALLELE}.txt )"

rm -f qual_scores_temp1_${ALLELE}.txt

# Compute average quality score of (confirmed) erroneous variant calls

awk -F'\t' '$6 == "YES" {print $4}' ${OUT_ERRORS} > qual_scores_temp2_${ALLELE}.txt

MEAN_QUAL_PY=${DATA}Simulations_CYP2C9/Python-codes/vcf_average_qual.py
MEAN_QUAL_ERR="$( python ${MEAN_QUAL_PY} qual_scores_temp2_${ALLELE}.txt )"

# Compute maximal quality score of (confirmed) erroneous variant calls

MAX_QUAL_ERR="$( python ${EXTREMUM_PY} max qual_scores_temp2_${ALLELE}.txt )"

rm -f qual_scores_temp2_${ALLELE}.txt

# Compute average read depth of called variants
source activate bcftools

DEPTH="$(bcftools query -f "[%DP]\n]" $VARCALL_VCF | awk '{sum+=$1}END{printf "%.1f\n", sum/NR}')"

conda deactivate

# Initiate stats file, if necessary
STATS=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_5x_sim_stats.tsv
if [ $ALLELE == 2 ]
then
    COLUMNS="ALLELE\t#ALLELIC_VARIANTS\t#ALLELIC_CALLED\tALL_ALLELIC_CALLED\tTOTAL_CALLED\tERROR_CALLED\t
    ERRORS_CALLED_CONF\tERRORS_CALLED_ALL\tMEAN_QUAL_ALLELIC\tMEAN_QUAL_ERRORS\tMIN_QUAL_ALLELIC\t
    MAX_QUAL_ERRORS\tMEAN_DEPTH\n"
    echo -e -n $COLUMNS > $STATS 
fi

# Write to stats file
ENTRY="${ALLELE}\t${NUM_VAR}\t${ALLELE_CALLED}\t${ALL_CALLED}\t${TOT_CALLED}\t${ERRORS_FOUND}\t
${ERROR_PERCENT}\t${ERROR_PERCENT_INC}\t${MEAN_QUAL_ALLELIC}\t${MEAN_QUAL_ERR}\t
${MIN_QUAL_ALLELIC}\t${MAX_QUAL_ERR}\t${DEPTH}\n"
echo -e -n $ENTRY >> ${STATS}
