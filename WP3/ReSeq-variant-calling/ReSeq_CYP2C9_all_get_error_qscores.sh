#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-CYP2C9-all-get-error-qscores # Name of job
#SBATCH --output=stdout-ReSeq-CYP2C9-all-get-error-qscores   # Standard output name
#SBATCH --error=stderror-ReSeq-CYP2C9-all-get-error-qscores    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=4     # Number of CPU cores
#SBATCH --time=3-00:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

## Set paths ##
DATA=/data/leuven/373/vsc37367/
IN_5X=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_5x/
IN_15X=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_15x/
IN_30X=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_30x/

## Initiate file ##
SUMMARY=${DATA}Simulations_CYP2C9/ReSeq/error_variant_qscores.tsv
echo -e "Coverage\tQScore" > ${SUMMARY}

## Get quality scores ##

source activate bcftools

# Get quality scores for 5x simulations #
for (( allele=1; allele<=85; allele++ )); do
    
    VCF=${IN_5X}ReSeq_CYP2C9_${allele}_VarCallErrors_norm.vcf
    bcftools query -f '5x\t%QUAL\n' ${VCF} >> ${SUMMARY}

done

# Get quality scores for 15x simulations #
for (( allele=1; allele<=85; allele++ )); do
    
    VCF=${IN_15X}ReSeq_CYP2C9_${allele}_VarCallErrors_norm.vcf
    bcftools query -f '15x\t%QUAL\n' ${VCF} >> ${SUMMARY}

done

# Get quality scores for 30x simulations #
for (( allele=1; allele<=85; allele++ )); do
    
    VCF=${IN_30X}ReSeq_CYP2C9_${allele}_VarCallErrors_norm.vcf
    bcftools query -f '30x\t%QUAL\n' ${VCF} >> ${SUMMARY}

done
