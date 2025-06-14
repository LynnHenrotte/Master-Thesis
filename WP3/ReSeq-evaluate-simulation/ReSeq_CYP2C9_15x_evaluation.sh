#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-CYP2C9-15x-evaluation # Name of job
#SBATCH --output=stdout-ReSeq-CYP2C9-15x-evaluation   # Standard output name
#SBATCH --error=stderror-ReSeq-CYP2C9-15x-evaluation    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=8     # Number of CPU cores
#SBATCH --time=0-05:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

DATA=/data/leuven/373/vsc37367/
REF=${DATA}references/chr10.fa
OUT=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_15x_evaluation/
IN=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_15x/

source activate bbmap

for ((ALLELE = 1; ALLELE <= 85; ALLELE++)); do

    bbmap.sh ref=$REF in=${IN}ReSeq_CYP2C9_${ALLELE}_left.fq in2=${IN}ReSeq_CYP2C9_${ALLELE}_right.fq \
        out=temp.sam statsfile=${OUT}CYP2C9_${ALLELE}_map_stats.txt \
        mhist=${OUT}CYP2C9_${ALLELE}_error_rates_per_base.txt \
        qhist=${OUT}CYP2C9_${ALLELE}_average_qscore_per_base.txt \
        gchist=${OUT}CYP2C9_${ALLELE}_GC_content.txt \
        ihist=${OUT}CYP2C9_${ALLELE}_insert_sizes.txt

done

rm -f temp.sam
rm -rf ref/
