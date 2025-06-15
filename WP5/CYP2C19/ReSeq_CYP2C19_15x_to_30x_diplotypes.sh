#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-diplotypes-CYP2C9-15x-to-30x_%a # Name of job
#SBATCH --output=stdout-ReSeq-diplotypes-CYP2C19-15x-to-30x_%a   # Standard output name
#SBATCH --error=stderror-ReSeq-diplotypes-CYP2C19-15x-to-30x_%a    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=4     # Number of CPU cores
#SBATCH --time=0-01:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
#SBATCH --array=1-19,22-26,28-39
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

# Set paths and parameters
DATA=/data/leuven/373/vsc37367/
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
BAM_IN=${STAGING}Simulations_CYP2C19/ReSeq_CYP2C19_15x_cut/
OUT=${STAGING}ReSeq_CYP2C19_diplotypes_30x/
ALLELE1=$SLURM_ARRAY_TASK_ID

source activate samtools

# Create diplotypes 
for ALLELE2 in $(seq ${ALLELE1} 39)
do
    # Verify that the allele exists
    if [[ $ALLELE2 -eq 20 || $ALLELE2 -eq 21 || $ALLELE2 -eq 27 ]]; then continue; fi

    if [[ $ALLELE1 -eq 36 && $ALLELE2 -eq 36 ]]; then
        for sub1 in {1..2}; do
            for ((sub2=$sub1; sub2<=2; sub2++)); do
                # Get simulated bam files for chosen alleles
                BAM1=${BAM_IN}ReSeq_CYP2C19_${ALLELE1}_00${sub1}_aligned.bam
                BAM2=${BAM_IN}ReSeq_CYP2C19_${ALLELE2}_00${sub2}_aligned2.bam

                # Merge haplotype bam files into diplotype bam file
                samtools merge -c -o ${OUT}CYP2C19_diplo_${ALLELE1}_00${sub1}_${ALLELE2}_00${sub2}.bam $BAM1 $BAM2

                # Index the merged bam file
                samtools index ${OUT}CYP2C19_diplo_${ALLELE1}_00${sub1}_${ALLELE2}_00${sub2}.bam
            done
        done
    elif [[ $ALLELE1 -eq 36 && $ALLELE2 -eq 37 ]]; then
        for sub1 in {1..2}; do
            for sub2 in {1..5}; do
                # Get simulated bam files for chosen alleles
                BAM1=${BAM_IN}ReSeq_CYP2C19_${ALLELE1}_00${sub1}_aligned.bam
                BAM2=${BAM_IN}ReSeq_CYP2C19_${ALLELE2}_00${sub2}_aligned2.bam

                # Merge haplotype bam files into diplotype bam file
                samtools merge -c -o ${OUT}CYP2C19_diplo_${ALLELE1}_00${sub1}_${ALLELE2}_00${sub2}.bam $BAM1 $BAM2

                # Index the merged bam file
                samtools index ${OUT}CYP2C19_diplo_${ALLELE1}_00${sub1}_${ALLELE2}_00${sub2}.bam
            done
        done
    elif [[ $ALLELE1 -eq 37 && $ALLELE2 -eq 37 ]]; then
        for sub1 in {1..5}; do
            for ((sub2=$sub1; sub2<=5; sub2++)); do
                # Get simulated bam files for chosen alleles
                BAM1=${BAM_IN}ReSeq_CYP2C19_${ALLELE1}_00${sub1}_aligned.bam
                BAM2=${BAM_IN}ReSeq_CYP2C19_${ALLELE2}_00${sub2}_aligned2.bam

                # Merge haplotype bam files into diplotype bam file
                samtools merge -c -o ${OUT}CYP2C19_diplo_${ALLELE1}_00${sub1}_${ALLELE2}_00${sub2}.bam $BAM1 $BAM2

                # Index the merged bam file
                samtools index ${OUT}CYP2C19_diplo_${ALLELE1}_00${sub1}_${ALLELE2}_00${sub2}.bam
            done
        done
    elif [[ $ALLELE1 -eq 36 ]]; then
        for sub1 in {1..2}; do
            # Get simulated bam files for chosen alleles
            BAM1=${BAM_IN}ReSeq_CYP2C19_${ALLELE1}_00${sub1}_aligned.bam
            BAM2=${BAM_IN}ReSeq_CYP2C19_${ALLELE2}_aligned2.bam

            # Merge haplotype bam files into diplotype bam file
            samtools merge -c -o ${OUT}CYP2C19_diplo_${ALLELE1}_00${sub1}_${ALLELE2}.bam $BAM1 $BAM2

            # Index the merged bam file
            samtools index ${OUT}CYP2C19_diplo_${ALLELE1}_00${sub1}_${ALLELE2}.bam
        done
    elif [[ $ALLELE2 -eq 36 ]]; then
        for sub2 in {1..2}; do
            # Get simulated bam files for chosen alleles
            BAM1=${BAM_IN}ReSeq_CYP2C19_${ALLELE1}_aligned.bam
            BAM2=${BAM_IN}ReSeq_CYP2C19_${ALLELE2}_00${sub2}_aligned2.bam

            # Merge haplotype bam files into diplotype bam file
            samtools merge -c -o ${OUT}CYP2C19_diplo_${ALLELE1}_${ALLELE2}_00${sub2}.bam $BAM1 $BAM2

            # Index the merged bam file
            samtools index ${OUT}CYP2C19_diplo_${ALLELE1}_${ALLELE2}_00${sub2}.bam
        done
    elif [[ $ALLELE1 -eq 37 ]]; then
        for sub1 in {1..5}; do
            # Get simulated bam files for chosen alleles
            BAM1=${BAM_IN}ReSeq_CYP2C19_${ALLELE1}_00${sub1}_aligned.bam
            BAM2=${BAM_IN}ReSeq_CYP2C19_${ALLELE2}_aligned2.bam

            # Merge haplotype bam files into diplotype bam file
            samtools merge -c -o ${OUT}CYP2C19_diplo_${ALLELE1}_00${sub1}_${ALLELE2}.bam $BAM1 $BAM2

            # Index the merged bam file
            samtools index ${OUT}CYP2C19_diplo_${ALLELE1}_00${sub1}_${ALLELE2}.bam
        done
    elif [[ $ALLELE2 -eq 37 ]]; then
        for sub2 in {1..5}; do
            # Get simulated bam files for chosen alleles
            BAM1=${BAM_IN}ReSeq_CYP2C19_${ALLELE1}_aligned.bam
            BAM2=${BAM_IN}ReSeq_CYP2C19_${ALLELE2}_00${sub2}_aligned2.bam

            # Merge haplotype bam files into diplotype bam file
            samtools merge -c -o ${OUT}CYP2C19_diplo_${ALLELE1}_${ALLELE2}_00${sub2}.bam $BAM1 $BAM2

            # Index the merged bam file
            samtools index ${OUT}CYP2C19_diplo_${ALLELE1}_${ALLELE2}_00${sub2}.bam
        done
    else
        # Get simulated bam files for chosen alleles
        BAM1=${BAM_IN}ReSeq_CYP2C19_${ALLELE1}_aligned.bam
        BAM2=${BAM_IN}ReSeq_CYP2C19_${ALLELE2}_aligned2.bam
        
        # Merge haplotype bam files into diplotype bam file
        samtools merge -c -o ${OUT}CYP2C19_diplo_${ALLELE1}_${ALLELE2}.bam $BAM1 $BAM2

        # Index the merged bam file
        samtools index ${OUT}CYP2C19_diplo_${ALLELE1}_${ALLELE2}.bam
    fi
done
