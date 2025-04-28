#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-Align-CYP2C9_%a # Name of job
#SBATCH --output=stdout-ReSeq-Align-CYP2C9_%a   # Standard output name
#SBATCH --error=stderror-ReSeq-Align-CYP2C9_%a    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=4     # Number of CPU cores
#SBATCH --time=3-00:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
#SBATCH --array=2-85
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

# Define paths and parameters
ALLELE=$SLURM_ARRAY_TASK_ID
REF=/data/leuven/373/vsc37367/references/chr10.fa
READS1=/data/leuven/373/vsc37367/Simulations/CYP2C9/ReSeq_CYP2C9/ReSeq_CYP2C9_${ALLELE}_left.fq
READS2=/data/leuven/373/vsc37367/Simulations/CYP2C9/ReSeq_CYP2C9/ReSeq_CYP2C9_${ALLELE}_right.fq
OUT_FILE=/data/leuven/373/vsc37367/Simulations/CYP2C9/ReSeq_CYP2C9/ReSeq_CYP2C9_${ALLELE}_aligned

# Alignment 
source activate bwa-mem

bwa mem -t 4 -R "@RG\tID:CYP2C9_${ALLELE}\tSM:CYP2C9_${ALLELE}" $REF $READS1 $READS2 > ${OUT_FILE}.sam

conda deactivate

# Convert SAM file to sorted BAM file
source activate samtools

samtools sort ${OUT_FILE}.sam -o ${OUT_FILE}.bam
samtools index ${OUT_FILE}.bam

# Obtain some alignment statistics
samtools flagstat ${OUT_FILE}.bam > ${OUT_FILE}_Stats.txt

conda deactivate

# Use GATK's haplotypecaller to verify whether the simulated variants are detected
OUT_VCF=/data/leuven/373/vsc37367/Simulations/CYP2C9/ReSeq_CYP2C9/ReSeq_CYP2C9_${ALLELE}_VarCall.vcf
source activate gatk4

gatk HaplotypeCaller -I ${OUT_FILE}.bam -O $OUT_VCF -R $REF

conda deactivate
