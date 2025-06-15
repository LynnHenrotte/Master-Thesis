#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=adjust-reference-CYP2C19-36-37  # Name of job
#SBATCH --output=stdout-adjust-reference-CYP2C19-36-37    # Standard output name
#SBATCH --error=stderror-adjust-reference-CYP2C19-36-37     # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=4     # Number of CPU cores
#SBATCH --time=0-05:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

DATA=/data/leuven/373/vsc37367/
BED_FILES=${DATA}Simulations_CYP2C19/SV_beds/
OUT_DIR=${DATA}Simulations_CYP2C19/VISOR_fastas/
REF=${DATA}references/chr10.fa

source activate visor

for bed in ${BED_FILES}*.bed
do
    # Adjust reference based on BED
    out_name=$( basename $bed .bed )
    VISOR HACk -b ${bed} -g ${REF} -o visor_${out_name}

    # File management
    mv visor_${out_name}/h1.fa ${OUT_DIR}${out_name}.fa
    mv visor_${out_name}/h1.fa.fai ${OUT_DIR}${out_name}.fa.fai
    rm -rf visor_${out_name}

done
