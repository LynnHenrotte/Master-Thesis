#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-CYP2C19_simulate_stellarpgx_regions_30x # Name of job
#SBATCH --output=stdout-ReSeq-CYP2C19_simulate_stellarpgx_regions_30x   # Standard output name
#SBATCH --error=stderror-ReSeq-CYP2C19_simulate_stellarpgx_regions_30x    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=8     # Number of CPU cores
#SBATCH --time=0-10:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

# Set paths
DATA=/data/leuven/373/vsc37367/
STAT_FILE=${DATA}Simulators/ReSeq/Hs-Nova-TruSeq.reseq
PROB_FILE=${DATA}Simulators/ReSeq/Hs-Nova-TruSeq.reseq.ipf

STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
REF=/staging/leuven/stg_00156/references/hg38.fa
OUT=${STAGING}Simulations_CYP2C19/stellarpgx_regions/

### Extract StellarPGx regions from hg38 ###
source activate samtools

samtools faidx $REF "chr12:47840000-47905100" -o ${OUT}stellarpgx_chr12.fa -n 80
samtools faidx $REF "chr7:55019000-55211700" -o ${OUT}stellarpgx_chr7.fa -n 80

conda deactivate

### Simulate reads from extracted regions ###

source activate reseq 

reseq illuminaPE -R ${OUT}stellarpgx_chr7.fa -s $STAT_FILE -p $PROB_FILE --coverage 30 \
    --ipfIterations 0 -1 ${OUT}stellarpgx_chr7_region_left.fq -2 ${OUT}stellarpgx_chr7_region_right.fq 

reseq illuminaPE -R ${OUT}stellarpgx_chr12.fa -s $STAT_FILE -p $PROB_FILE --coverage 30 \
    --ipfIterations 0 -1 ${OUT}stellarpgx_chr12_region_left.fq -2 ${OUT}stellarpgx_chr12_region_right.fq 

conda deactivate

### Alignment ###

source activate bwa-mem

bwa mem -t 4 -R "@RG\tID:CYP2C19_sim\tSM:CYP2C19_sim" $REF \
    ${OUT}stellarpgx_chr7_region_left.fq \
    ${OUT}stellarpgx_chr7_region_right.fq > ${OUT}stellarpgx_chr7_region.sam

bwa mem -t 4 -R "@RG\tID:CYP2C19_sim\tSM:CYP2C19_sim" $REF \
    ${OUT}stellarpgx_chr12_region_left.fq \
    ${OUT}stellarpgx_chr12_region_right.fq > ${OUT}stellarpgx_chr12_region.sam

conda deactivate

### Convert to BAM and index ###

source activate samtools

samtools view -h -o ${OUT}stellarpgx_chr7_region.bam ${OUT}stellarpgx_chr7_region.sam
samtools view -h -o ${OUT}stellarpgx_chr12_region.bam ${OUT}stellarpgx_chr12_region.sam

rm -f ${OUT}stellarpgx_chr*_region.sam

samtools sort ${OUT}stellarpgx_chr7_region.bam -o ${OUT}stellarpgx_chr7_region.bam
samtools index ${OUT}stellarpgx_chr7_region.bam

samtools sort ${OUT}stellarpgx_chr12_region.bam -o ${OUT}stellarpgx_chr12_region.bam
samtools index ${OUT}stellarpgx_chr12_region.bam

conda deactivate
