#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-diplotypes-CYP2C9-30x-to-60x_%a # Name of job
#SBATCH --output=stdout-ReSeq-diplotypes-CYP2C9-30x-to-60x_%a   # Standard output name
#SBATCH --error=stderror-ReSeq-diplotypes-CYP2C9-30x-to-60x_%a    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=4     # Number of CPU cores
#SBATCH --time=0-01:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
#SBATCH --array=1-85
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

# Set paths and parameters
DATA=/data/leuven/373/vsc37367/
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
OUT=${STAGING}ReSeq_diplotypes_60x/
ALLELE1=$SLURM_ARRAY_TASK_ID

source activate samtools

# Create diplotypes 
for ALLELE2 in $(seq ${ALLELE1} 85)
do

    # Get simulated bam files for chosen alleles
    if [ $ALLELE1 -eq $ALLELE2 ]
    then
        BAM1=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_30x/ReSeq_CYP2C9_${ALLELE1}_aligned.bam
        BAM2=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_30x/ReSeq_CYP2C9_${ALLELE1}_aligned2.bam
    else        
        BAM1=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_30x/ReSeq_CYP2C9_${ALLELE1}_aligned.bam
        BAM2=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_30x/ReSeq_CYP2C9_${ALLELE2}_aligned.bam
    fi

    # Merge haplotype bam files into diplotype bam file
    samtools merge -c -o ${OUT}CYP2C9_diplo_${ALLELE1}_${ALLELE2}_new.bam $BAM1 $BAM2

    # Index the merged bam file
    samtools index ${OUT}CYP2C9_diplo_${ALLELE1}_${ALLELE2}_new.bam

done
