#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-CYP2C19-15x-SNP-calling_%a # Name of job
#SBATCH --output=stdout-ReSeq-CYP2C19-15x-SNP-calling_%a   # Standard output name
#SBATCH --error=stderror-ReSeq-CYP2C19-15x-SNP-calling_%a    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=4     # Number of CPU cores
#SBATCH --time=3-00:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
#SBATCH --array=2-19,22-26,28-35,39        # First run this with array=1; when finished, run with array=2-19,22-26,28-35,39
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

# Set paths and parameters
ALLELE=$SLURM_ARRAY_TASK_ID
DATA=/data/leuven/373/vsc37367/
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
REF=${DATA}references/chr10.fa
OUT=${STAGING}Simulations_CYP2C19/ReSeq_CYP2C19_15x_SNP_calls/
ALLELE_VCFS=${DATA}CYP2C19-6.2.3/GRCh38/
ALLELE_VCF=${ALLELE_VCFS}CYP2C19_${ALLELE}.vcf

## Variant calling with HaplotypeCaller
VARCALL_VCF1=${OUT}ReSeq_CYP2C19_${ALLELE}_VarCall.vcf
VARCALL_VCF2=${OUT}ReSeq_CYP2C19_${ALLELE}_VarCall2.vcf
BAM1=${STAGING}Simulations_CYP2C19/ReSeq_CYP2C19_15x/ReSeq_CYP2C19_${ALLELE}_aligned
BAM2=${STAGING}Simulations_CYP2C19/ReSeq_CYP2C19_15x/ReSeq_CYP2C19_${ALLELE}_aligned2

source activate gatk4

gatk HaplotypeCaller -I ${BAM1}.bam -O $VARCALL_VCF1 -R $REF
gatk HaplotypeCaller -I ${BAM2}.bam -O $VARCALL_VCF2 -R $REF

conda deactivate 

## Prepare star allele VCF ##
source activate bcftools

if ! [ -f ${ALLELE_VCFS}CYP2C19_${ALLELE}_norm.vcf ]
then
    # Compress original vcf file
    bgzip -k -o $ALLELE_VCF.gz $ALLELE_VCF

    # Normalize vcf file
    ALLELE_VCF_NORM=${ALLELE_VCFS}CYP2C19_${ALLELE}_norm.vcf
    bcftools norm -f $REF -o $ALLELE_VCF_NORM $ALLELE_VCF.gz

    # Compress normalized vcf file
    bgzip -k -o $ALLELE_VCF_NORM.gz $ALLELE_VCF_NORM

    # Index vcf file
    bcftools index -o $ALLELE_VCF_NORM.gz.csi $ALLELE_VCF_NORM.gz

    # Remove non-normalized compressed and indexed vcf file
    rm -f $ALLELE_VCF.gz*
else 
    ALLELE_VCF_NORM=${ALLELE_VCFS}CYP2C19_${ALLELE}_norm.vcf
fi

## Prepare variant call VCFs ##

# Compress original vcf file
bgzip -k -o $VARCALL_VCF1.gz $VARCALL_VCF1
bgzip -k -o $VARCALL_VCF2.gz $VARCALL_VCF2

# Remove called variants that are outside of the CYP2C locus in variant call vcf
#bcftools index -o $VARCALL_VCF1.gz.csi $VARCALL_VCF1.gz
#bcftools view -r chr10:94682000-95071000 -o $VARCALL_VCF1 $VARCALL_VCF1.gz
#bgzip -k -o $VARCALL_VCF1.gz $VARCALL_VCF1
#rm -f $VARCALL_VCF1.gz.csi
#bcftools index -o $VARCALL_VCF2.gz.csi $VARCALL_VCF2.gz
#bcftools view -r chr10:94682000-95071000 -o $VARCALL_VCF2 $VARCALL_VCF2.gz
#bgzip -k -o $VARCALL_VCF2.gz $VARCALL_VCF2
#rm -f $VARCALL_VCF2.gz.csi

# Normalize vcf file
VARCALL_VCF_NORM1=${OUT}ReSeq_CYP2C19_${ALLELE}_VarCall_norm.vcf
bcftools norm -f $REF -o $VARCALL_VCF_NORM1 $VARCALL_VCF1.gz
VARCALL_VCF_NORM2=${OUT}ReSeq_CYP2C19_${ALLELE}_VarCall_norm2.vcf
bcftools norm -f $REF -o $VARCALL_VCF_NORM2 $VARCALL_VCF2.gz

# Compress normalized vcf file
bgzip -k -o $VARCALL_VCF_NORM1.gz $VARCALL_VCF_NORM1
bgzip -k -o $VARCALL_VCF_NORM2.gz $VARCALL_VCF_NORM2

# Index vcf file
bcftools index -o $VARCALL_VCF_NORM1.gz.csi $VARCALL_VCF_NORM1.gz
bcftools index -o $VARCALL_VCF_NORM2.gz.csi $VARCALL_VCF_NORM2.gz

# Remove non-normalized compressed and indexed vcf file
rm -f $VARCALL_VCF1.gz*
rm -f $VARCALL_VCF2.gz*

## Verify that the required variants have been simulated/called ##

# Find variants that are shared by both vcf files
bcftools isec -n =2 -p ${OUT}ReSeq_CYP2C19_${ALLELE}_compare_vcf $ALLELE_VCF_NORM.gz $VARCALL_VCF_NORM1.gz
bcftools isec -n =2 -p ${OUT}ReSeq_CYP2C19_${ALLELE}_compare_vcf2 $ALLELE_VCF_NORM.gz $VARCALL_VCF_NORM2.gz

# Find called variants that are not part of the allele
bcftools isec -C -w 1 -o ${OUT}ReSeq_CYP2C19_${ALLELE}_VarCallErrors_norm.vcf $VARCALL_VCF_NORM1.gz $ALLELE_VCF_NORM.gz
VARCALL_ERROR_VCF1=${OUT}ReSeq_CYP2C19_${ALLELE}_VarCallErrors_norm.vcf
bcftools isec -C -w 1 -o ${OUT}ReSeq_CYP2C19_${ALLELE}_VarCallErrors_norm2.vcf $VARCALL_VCF_NORM2.gz $ALLELE_VCF_NORM.gz
VARCALL_ERROR_VCF2=${OUT}ReSeq_CYP2C19_${ALLELE}_VarCallErrors_norm2.vcf

conda deactivate
    
## Summarize information on allele ##

source activate gatk4

# Count variants of the allele
NUM_VAR="$( gatk CountVariants -V $ALLELE_VCF | sed -n 2p )"

# Count total variants called by HaplotypeCaller
TOT_CALLED1="$( gatk CountVariants -V $VARCALL_VCF_NORM1 | sed -n 2p )"
TOT_CALLED2="$( gatk CountVariants -V $VARCALL_VCF_NORM2 | sed -n 2p )"

# Count number of allelic variants called by HaplotypeCaller
COMP_FILE1=${OUT}ReSeq_CYP2C19_${ALLELE}_compare_vcf/sites.txt
ALLELE_CALLED1=$( wc -l < "${COMP_FILE1}" )
COMP_FILE2=${OUT}ReSeq_CYP2C19_${ALLELE}_compare_vcf2/sites.txt
ALLELE_CALLED2=$( wc -l < "${COMP_FILE2}" )

# Count number of erroneous variants called by HaplotypeCaller
ERROR_CALLED1="$( gatk CountVariants -V $VARCALL_ERROR_VCF1 | sed -n 2p )"
ERROR_CALLED2="$( gatk CountVariants -V $VARCALL_ERROR_VCF2 | sed -n 2p )"

conda deactivate

# Compute binary indicator for allelic variants called
if [ $ALLELE_CALLED1 == $NUM_VAR ]; then 
    ALL_CALLED1="YES"
else
    ALL_CALLED1="NO"
fi

if [ $ALLELE_CALLED2 == $NUM_VAR ]; then 
    ALL_CALLED2="YES"
else
    ALL_CALLED2="NO"
fi

# Compute average quality score of allelic variant calls
source activate bcftools

ALLELE_CALLS1=${OUT}ReSeq_CYP2C19_${ALLELE}_compare_vcf/0001.vcf
bcftools query -f '%QUAL\n' $ALLELE_CALLS1 > qual_scores_temp1_V1_${ALLELE}.txt
ALLELE_CALLS2=${OUT}ReSeq_CYP2C19_${ALLELE}_compare_vcf2/0001.vcf
bcftools query -f '%QUAL\n' $ALLELE_CALLS2 > qual_scores_temp1_V2_${ALLELE}.txt

conda deactivate

source activate python

MEAN_QUAL_PY=${STAGING}Simulations_CYP2C19/Python-codes/vcf_average_qual.py
MEAN_QUAL_ALLELIC1="$( python ${MEAN_QUAL_PY} qual_scores_temp1_V1_${ALLELE}.txt )"
MEAN_QUAL_ALLELIC2="$( python ${MEAN_QUAL_PY} qual_scores_temp1_V2_${ALLELE}.txt )"

# Compute minimal quality score of allelic variant calls

EXTREMUM_PY=${STAGING}Simulations_CYP2C19/Python-codes/extremum.py
MIN_QUAL_ALLELIC1="$( python ${EXTREMUM_PY} min qual_scores_temp1_V1_${ALLELE}.txt )"
MIN_QUAL_ALLELIC2="$( python ${EXTREMUM_PY} min qual_scores_temp1_V2_${ALLELE}.txt )"

rm -f qual_scores_temp1_*_${ALLELE}.txt

conda deactivate

# Compute average quality score of erroneous variant calls
source activate bcftools

ERROR_CALLS1=${VARCALL_ERROR_VCF1}
bcftools query -f '%QUAL\n' $ERROR_CALLS1 > qual_scores_temp2_V1_${ALLELE}.txt
ERROR_CALLS2=${VARCALL_ERROR_VCF2}
bcftools query -f '%QUAL\n' $ERROR_CALLS2 > qual_scores_temp2_V2_${ALLELE}.txt

conda deactivate
source activate python

MEAN_QUAL_PY=${STAGING}Simulations_CYP2C19/Python-codes/vcf_average_qual.py
MEAN_QUAL_ERR1="$( python ${MEAN_QUAL_PY} qual_scores_temp2_V1_${ALLELE}.txt )"
MEAN_QUAL_ERR2="$( python ${MEAN_QUAL_PY} qual_scores_temp2_V2_${ALLELE}.txt )"

# Compute maximal quality score of erroneous variant calls
MAX_QUAL_ERR1="$( python ${EXTREMUM_PY} max qual_scores_temp2_V1_${ALLELE}.txt )"
MAX_QUAL_ERR2="$( python ${EXTREMUM_PY} max qual_scores_temp2_V2_${ALLELE}.txt )"

rm -f qual_scores_temp2_*_${ALLELE}.txt
conda deactivate

# Compute average read depth of called variants
#source activate bcftools

#DEPTH="$(bcftools query -f "[%DP]\n]" $VARCALL_VCF | awk '{sum+=$1}END{printf "%.1f\n", sum/NR}')"

#conda deactivate

# Initiate stats file, if necessary
STATS=${STAGING}Simulations_CYP2C19/ReSeq_CYP2C19_15x_SNP_call_stats.tsv
if [ $ALLELE == 1 ]
then
    COLUMNS="Allele\t#Allelic-Variants\t#Allelic-Called\tAll-Allelic-Called\tTotal-Called\tError_Called\t
    Mean-Allelic-Quality\tMean-Error-Quality\tMin-Allelic-Quality\tMax-Error-Quality\n"
    echo -e -n $COLUMNS > $STATS 
fi

# Write to stats file
ENTRY1="${ALLELE}\t${NUM_VAR}\t${ALLELE_CALLED1}\t${ALL_CALLED1}\t${TOT_CALLED1}\t${ERROR_CALLED1}\t
${MEAN_QUAL_ALLELIC1}\t${MEAN_QUAL_ERR1}\t${MIN_QUAL_ALLELIC1}\t${MAX_QUAL_ERR1}\n"
ENTRY2="${ALLELE}\t${NUM_VAR}\t${ALLELE_CALLED2}\t${ALL_CALLED2}\t${TOT_CALLED2}\t${ERROR_CALLED2}\t
${MEAN_QUAL_ALLELIC2}\t${MEAN_QUAL_ERR2}\t${MIN_QUAL_ALLELIC2}\t${MAX_QUAL_ERR2}\n"
echo -e -n $ENTRY1 >> ${STATS}
echo -e -n $ENTRY2 >> ${STATS}
