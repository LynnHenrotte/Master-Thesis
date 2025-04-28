#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=simulate-ReSeq-CYP2C9_%a  # Name of job
#SBATCH --output=stdout-simulate-ReSeq-CYP2C9_%a   # Standard output name
#SBATCH --error=stderror-simulate-ReSeq-CYP2C9_%a    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=4     # Number of CPU cores
#SBATCH --time=3-00:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
#SBATCH --array=2-95
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

# Set parameters and paths for VCF file preparation
ALLELE=$SLURM_ARRAY_TASK_ID
DATA=/data/leuven/373/vsc37367
VCF_FILES=/data/leuven/373/vsc37367/CYP2C9-6.2.3/GRCh38
OUT=/data/leuven/373/vsc37367/Simulations/CYP2C9/ReSeq_CYP2C9

# Remove unnecessary contig lines and add genotype columns
python ${DATA}/vcf_for_reseq.py ${VCF_FILES}/CYP2C9_${ALLELE}.vcf

# Shift position in VCF file, beginning at CYP2C locus
python ${DATA}/vcf_shift.py CYP2C9_${ALLELE}_reseq.vcf
rm -f CYP2C9_${ALLELE}_reseq.vcf
mv CYP2C9_${ALLELE}_reseq_shift.vcf $OUT

# Set paths for simulation
REF=/data/leuven/373/vsc37367/references/cyp2c.fa
STAT_FILE=/data/leuven/373/vsc37367/Simulators/ReSeq/Hs-Nova-TruSeq.reseq
PROB_FILE=/data/leuven/373/vsc37367/Simulators/ReSeq/Hs-Nova-TruSeq.reseq.ipf
VCF=${OUT}/CYP2C9_${ALLELE}_reseq_shift.vcf

# Simulate reads with ReSeq
source activate reseq

reseq illuminaPE -R $REF -s $STAT_FILE -p $PROB_FILE -V $VCF --coverage 30 --ipfIterations 0 \
    -1 ${OUT}/ReSeq_CYP2C9_${ALLELE}_left.fq -2 ${OUT}/ReSeq_CYP2C9_${ALLELE}_right.fq \
    --writeSysError ${OUT}/ReSeq_CYP2C9_${ALLELE}_errors.fq

conda deactivate