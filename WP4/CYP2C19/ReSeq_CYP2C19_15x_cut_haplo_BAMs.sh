#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-CYP2C9-15x-cut-haplo-BAMs_%a # Name of job
#SBATCH --output=stdout-ReSeq-CYP2C19-15x-cut-haplo-BAMs_%a   # Standard output name
#SBATCH --error=stderror-ReSeq-CYP2C19-15x-cut-haplo-BAMs_%a    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=4     # Number of CPU cores
#SBATCH --time=0-01:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
#SBATCH --array=36,37
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
BAM_FILES=${STAGING}Simulations_CYP2C19/ReSeq_CYP2C19_15x/
OUT=${STAGING}Simulations_CYP2C19/ReSeq_CYP2C19_15x_cut/
ALLELE=$SLURM_ARRAY_TASK_ID

source activate samtools

if [[ $ALLELE == 36 ]]
then
    for sub in 1 2; do
        # First set of simulated reads
        BAM_IN1=${BAM_FILES}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned.bam
        BAM_OUT1=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned.bam
        samtools view -h $BAM_IN1 "chr10:94400000-95071000" > ${BAM_OUT1}
        samtools sort ${BAM_OUT1} -o ${BAM_OUT1}
        samtools index ${BAM_OUT1}

        # Second set of simulated reads
        BAM_IN2=${BAM_FILES}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned2.bam
        BAM_OUT2=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned2.bam
        samtools view -h $BAM_IN2 "chr10:94400000-95071000" > ${BAM_OUT2}
        samtools sort ${BAM_OUT2} -o ${BAM_OUT2}
        samtools index ${BAM_OUT2}
    done
elif [[ $ALLELE == 37 ]]
then
    for sub in {1..5}; do
        # First set of simulated reads
        BAM_IN1=${BAM_FILES}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned.bam
        BAM_OUT1=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned.bam
        samtools view -h $BAM_IN1 "chr10:94400000-95071000" > ${BAM_OUT1}
        samtools sort ${BAM_OUT1} -o ${BAM_OUT1}
        samtools index ${BAM_OUT1}

        # Second set of simulated reads
        BAM_IN2=${BAM_FILES}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned2.bam
        BAM_OUT2=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned2.bam
        samtools view -h $BAM_IN2 "chr10:94400000-95071000" > ${BAM_OUT2}
        samtools sort ${BAM_OUT2} -o ${BAM_OUT2}
        samtools index ${BAM_OUT2}
    done
fi
