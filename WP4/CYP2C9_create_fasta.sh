#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=create-fasta-CYP2C9_%a  # Name of job
#SBATCH --output=stdout-create-fasta-CYP2C9_%a   # Standard output name
#SBATCH --error=stderror-create-fasta-CYP2C9_%a    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=4     # Number of CPU cores
#SBATCH --time=3-00:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
#SBATCH --array=2-85
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

# Define parameters and paths
ALLELE=$SLURM_ARRAY_TASK_ID
REF=/data/leuven/373/vsc37367/references/chr10.fa
VCF=/data/leuven/373/vsc37367/CYP2C9-6.2.3/GRCh38/CYP2C9_${ALLELE}.vcf
OUT=/data/leuven/373/vsc37367/Simulations/CYP2C9/Fasta_files/
SIMUG_OUT=/data/leuven/373/vsc37367/Simulations/CYP2C9/simuG_outfiles/
SIMUG=/data/leuven/373/vsc37367/Simulators/simuG/simuG.pl

# Modify wild type reference sequence
source activate simuG

perl $SIMUG -refseq $REF -snp_vcf $VCF -indel_vcf $VCF -prefix CYP2C9_${ALLELE}

conda deactivate

# Extract CYP2C locus from modified chr10 sequence and put into new FASTA file
source activate samtools

samtools faidx CYP2C9_${ALLELE}.simseq.genome.fa
samtools faidx CYP2C9_${ALLELE}.simseq.genome.fa "chr10:94682000-95071000" -o CYP2C9_${ALLELE}_short.fa -n 80

conda deactivate

# Organize files
mv CYP2C9_${ALLELE}.simseq.genome.fa CYP2C9_${ALLELE}_short.fa $OUT
mv CYP2C9_${ALLELE}.refseq2simseq* $SIMUG_OUT
rm -f CYP2C9_${ALLELE}.simseq.genome.fa.fai
