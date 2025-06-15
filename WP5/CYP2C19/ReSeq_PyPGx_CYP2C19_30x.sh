#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-PyPGx-CYP2C19-30x_%a # Name of job
#SBATCH --output=stdout-ReSeq-PyPGx-CYP2C19-30x_%a   # Standard output name
#SBATCH --error=stderror-ReSeq-PyPGx-CYP2C19-30x_%a    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=8     # Number of CPU cores
#SBATCH --time=0-10:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
#SBATCH --array=1-19,22-26,28-39
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"
export PYPGX_BUNDLE="${STAGING}Diplotype_Callers/PyPGx/pypgx-bundle"

REF=/data/leuven/373/vsc37367/references/chr10.fa
OUT=${STAGING}Diplotype_Callers/PyPGx/
RESULTS=${OUT}ReSeq_PyPGx_CYP2C19_30x_calls/
BAMFILES=${STAGING}ReSeq_CYP2C19_diplotypes_30x/
ALLELE1=$SLURM_ARRAY_TASK_ID

source activate pypgx

for diplo_bam in ${BAMFILES}CYP2C19_diplo_${ALLELE1}_*.bam; do

    # Get filename without bam extension
    diplo_name="$( basename ${diplo_bam} .bam )"

    # Create the input VCF file from the BAM file
    pypgx create-input-vcf \
        PyPGx-${diplo_name}.vcf.gz ${REF} ${diplo_bam} \
        --assembly GRCh38 --genes CYP2C19 

    # Perform the diplotype calling with PyPGx
    pypgx run-ngs-pipeline CYP2C19 PyPGx-${diplo_name} --platform Targeted --assembly GRCh38 \
        --variants ${OUT}PyPGx-${diplo_name}.vcf.gz --force

    # Put away created files
    mv -f PyPGx-${diplo_name}.vcf.gz* ${STAGING}ReSeq_PyPGx_CYP2C19_30x_diplotype_VCFs/
    mv ${OUT}PyPGx-${diplo_name}/results.zip ${RESULTS}results_${diplo_name}.zip
    rm -rf ${OUT}PyPGx-${diplo_name}
    
done
