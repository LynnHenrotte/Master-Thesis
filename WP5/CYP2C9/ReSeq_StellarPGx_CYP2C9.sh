#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-StellarPGx-CYP2C9 # Name of job
#SBATCH --output=stdout-ReSeq-StellarPGx-CYP2C9   # Standard output name
#SBATCH --error=stderror-ReSeq-StellarPGx-CYP2C9    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=8     # Number of CPU cores
#SBATCH --time=0-10:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

# Set paths and parameters
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
MAIN=${STAGING}Diplotype_Callers/StellarPGx/main.nf

source activate stellarpgx

nextflow run ${MAIN} -profile standard --build hg38 --gene cyp2c9
