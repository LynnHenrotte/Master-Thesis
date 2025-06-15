#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-Aldy-CYP2C19-30x-%a # Name of job
#SBATCH --output=stdout-ReSeq-Aldy-CYP2C19-30x-%a   # Standard output name
#SBATCH --error=stderror-ReSeq-Aldy-CYP2C19-30x-%a    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=8     # Number of CPU cores
#SBATCH --time=0-05:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
#SBATCH --array=1-19,22-26,28-39
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

# Set paths and parameters
DATA=/data/leuven/373/vsc37367/
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
BAM_FILES=${STAGING}ReSeq_CYP2C19_diplotypes_30x/
RESULTS=${STAGING}Diplotype_Callers/Aldy/Aldy_ReSeq_CYP2C19_diplotype_calls_30x/
ALLELE1=$SLURM_ARRAY_TASK_ID

# Index all bam files in the directory, if not done yet
source activate samtools

for diplo_bam in ${BAM_FILES}*CYP2C19_diplo_${ALLELE1}_*.bam
do
    if ! test -f ${diplo_bam}.bai; then
        samtools index ${diplo_bam}
    fi
done

conda deactivate

# For each bam file in the directory, perform diplotype calling with aldy
source activate aldy

for diplo_bam in ${BAM_FILES}CYP2C19_diplo_${ALLELE1}_*.bam
do
    diplo_name=$( basename $diplo_bam .bam )

    if [[ $diplo_name == *36* || $diplo_name == *37* ]]; then
        aldy genotype -p wgs -g CYP2C19 -o ${RESULTS}${diplo_name}.aldy --genome hg38 --cn-neutral-region chr10:95036772-95069497 ${diplo_bam} 
    else
        if ! [ -f ${RESULTS}${diplo_name}.aldy ]; then
            aldy genotype -p wgs -g CYP2C19 -o ${RESULTS}${diplo_name}.aldy --genome hg38 -c 1,1 ${diplo_bam}    
        fi
    fi
done
