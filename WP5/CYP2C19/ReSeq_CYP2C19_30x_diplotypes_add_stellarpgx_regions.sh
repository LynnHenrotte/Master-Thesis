#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-CYP2C19_30x_diplotypes_add_stellarpgx_regions # Name of job
#SBATCH --output=stdout-ReSeq-CYP2C19_30x_diplotypes_add_stellarpgx_regions   # Standard output name
#SBATCH --error=stderror-ReSeq-CYP2C19_30x_diplotypes_add_stellarpgx_regions    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=8     # Number of CPU cores
#SBATCH --time=0-10:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

# Set paths
DATA=/data/leuven/373/vsc37367/
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
BAM_IN=${STAGING}ReSeq_CYP2C19_diplotypes_30x/
BAM_CHR7=${STAGING}Simulations_CYP2C19/stellarpgx_regions/stellarpgx_chr7_region.bam
BAM_CHR12=${STAGING}Simulations_CYP2C19/stellarpgx_regions/stellarpgx_chr12_region.bam
OUT=${STAGING}ReSeq_CYP2C19_stellarpgx_diplotypes_30x/

# Merge chr7 and chr12 region BAMs with diplotype BAMs
source activate samtools

for diplo_bam in ${BAM_IN}*.bam; do
    diplo_name=$( basename $diplo_bam .bam )
    samtools merge -c -o ${OUT}${diplo_name}.bam $diplo_bam $BAM_CHR7 $BAM_CHR12
done

# Sort and index merged diplotypes
for diplo_bam in ${OUT}*.bam; do
    samtools sort -o ${diplo_bam} ${diplo_bam}
    samtools index ${diplo_bam}
done

conda deactivate
