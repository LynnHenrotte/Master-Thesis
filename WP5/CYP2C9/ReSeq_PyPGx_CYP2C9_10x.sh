#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-PyPGx-CYP2C9-10x_%a # Name of job
#SBATCH --output=stdout-ReSeq-PyPGx-CYP2C9-10x_%a   # Standard output name
#SBATCH --error=stderror-ReSeq-PyPGx-CYP2C9-10x_%a    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=8     # Number of CPU cores
#SBATCH --time=0-04:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
#SBATCH --array=1-71
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"
export PYPGX_BUNDLE="${STAGING}Diplotype_Callers/PyPGx/pypgx-bundle"

REF=/data/leuven/373/vsc37367/references/chr10.fa
OUT=${STAGING}Diplotype_Callers/PyPGx/
RESULTS=${OUT}ReSeq_PyPGx_CYP2C9_10x_calls/
BAMFILES1=${STAGING}ReSeq_diplotypes_10x_part1/
BAMFILES2=${STAGING}ReSeq_diplotypes_10x_part2/
BAMFILES3=${STAGING}ReSeq_diplotypes_10x_part3/
ALLELE1=$SLURM_ARRAY_TASK_ID

source activate pypgx

for ((ALLELE2=${ALLELE1}; ALLELE2<=71; ALLELE2++))
do

    bamfile_name=CYP2C9_diplo_${ALLELE1}_${ALLELE2}
    if [ -f ${BAMFILES1}${bamfile_name}.bam ]
    then 
        BAMFILES=$BAMFILES1
    elif [ -f ${BAMFILES2}${bamfile_name}.bam ]
    then 
        BAMFILES=$BAMFILES2
    else
        BAMFILES=$BAMFILES3
    fi

    # Create the input VCF file from the BAM file
    pypgx create-input-vcf \
        PyPGx-CYP2C9-diplo-${ALLELE1}-${ALLELE2}.vcf.gz ${REF} ${BAMFILES}${bamfile_name}.bam \
        --assembly GRCh38 --genes CYP2C9 

    # Perform the diplotype calling with PyPGx
    pypgx run-ngs-pipeline CYP2C9 PyPGx_${bamfile_name} --platform Targeted --assembly GRCh38 \
        --variants ${OUT}PyPGx-CYP2C9-diplo-${ALLELE1}-${ALLELE2}.vcf.gz --force

    # Put away created files
    mv -f PyPGx-CYP2C9-diplo-${ALLELE1}-${ALLELE2}.vcf.gz* ${STAGING}ReSeq_PyPGx_CYP2C9_10x_diplotype_VCFs/
    mv ${OUT}PyPGx_CYP2C9_diplo_${ALLELE1}_${ALLELE2}/results.zip ${RESULTS}results_CYP2C9_${ALLELE1}_${ALLELE2}.zip
    rm -rf ${OUT}PyPGx_CYP2C9_diplo_${ALLELE1}_${ALLELE2}

done
