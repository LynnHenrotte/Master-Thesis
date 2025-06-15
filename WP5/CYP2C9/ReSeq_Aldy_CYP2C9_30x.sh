#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-Aldy-CYP2C9-30x-%a # Name of job
#SBATCH --output=stdout-ReSeq-Aldy-CYP2C9-30x-%a   # Standard output name
#SBATCH --error=stderror-ReSeq-Aldy-CYP2C9-30x-%a    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=8     # Number of CPU cores
#SBATCH --time=0-05:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
#SBATCH --array=1-3
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

# Set paths and parameters
DATA=/data/leuven/373/vsc37367/
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
BAM_FILES=${STAGING}ReSeq_diplotypes_30x_part${SLURM_ARRAY_TASK_ID}/
OUT=${STAGING}Diplotype_Callers/Aldy/Aldy_ReSeq_CYP2C9_diplotype_calls_30x/

# Index all bam files in the directory, if not done yet
source activate samtools

for bamfile in ${BAM_FILES}*.bam
do
    if ! test -f ${bamfile}.bai; then
        samtools index ${bamfile}
    fi
done

conda deactivate

# For each bam file in the directory, perform diplotype calling with aldy
source activate aldy

for bamfile in ${BAM_FILES}*.bam
do
    out_name=$( basename $bamfile .bam )
    if ! [ -f ${OUT}${out_name}.aldy ]; then
        aldy genotype -p wgs -g CYP2C9 -o ${OUT}${out_name}.aldy --genome hg38 -c 1,1 ${bamfile}    
    fi
done
