#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=simulate-ReSeq-CYP2C19_15x_%a  # Name of job
#SBATCH --output=stdout-simulate-ReSeq-CYP2C19_15x_%a   # Standard output name
#SBATCH --error=stderror-simulate-ReSeq-CYP2C19_15x_%a    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=8     # Number of CPU cores
#SBATCH --time=0-10:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
#SBATCH --array=1-19,22-26,28-39
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

##################
### SIMULATION ###
##################

# Set paths for simulation
ALLELE=$SLURM_ARRAY_TASK_ID
DATA=/data/leuven/373/vsc37367/
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
REF=${DATA}references/cyp2c_ext.fa
STAT_FILE=${DATA}Simulators/ReSeq/Hs-Nova-TruSeq.reseq
PROB_FILE=${DATA}Simulators/ReSeq/Hs-Nova-TruSeq.reseq.ipf
OUT=${STAGING}Simulations_CYP2C19/ReSeq_CYP2C19_15x/
FASTAS=${STAGING}Simulations_CYP2C19/VISOR_fastas/

source activate reseq

# CYP2C19*38 is the wild type!
if [ ${ALLELE} == 38 ]
then
    # Simulate two sets of reads from the reference without additional variants
    reseq illuminaPE -R $REF -s $STAT_FILE -p $PROB_FILE --coverage 15 --ipfIterations 0 \
        -1 ${OUT}ReSeq_CYP2C19_${ALLELE}_left.fq -2 ${OUT}ReSeq_CYP2C19_${ALLELE}_right.fq 
    
    reseq illuminaPE -R $REF -s $STAT_FILE -p $PROB_FILE --coverage 15 --ipfIterations 0 \
        -1 ${OUT}ReSeq_CYP2C19_${ALLELE}_left2.fq -2 ${OUT}ReSeq_CYP2C19_${ALLELE}_right2.fq 

elif [ ${ALLELE} == 36 ]
then
    # Simulate two sets of reads from the modified reference sequences with full gene deletions
    for sub in 1 2
    do
        fasta=${FASTAS}cyp2c19_36_${sub}.fa
        reseq illuminaPE -R $fasta -s $STAT_FILE -p $PROB_FILE --coverage 15 --ipfIterations 0 \
            -1 ${OUT}ReSeq_CYP2C19_36_00${sub}_left.fq -2 ${OUT}ReSeq_CYP2C19_36_00${sub}_right.fq
        reseq illuminaPE -R $fasta -s $STAT_FILE -p $PROB_FILE --coverage 15 --ipfIterations 0 \
            -1 ${OUT}ReSeq_CYP2C19_36_00${sub}_left2.fq -2 ${OUT}ReSeq_CYP2C19_36_00${sub}_right2.fq
    done

elif [ ${ALLELE} == 37 ]
then
    # Simulate two sets of reads from the modified reference sequences with partial gene deletions
    for sub in {1..5}
    do
        fasta=${FASTAS}cyp2c19_37_${sub}.fa
        reseq illuminaPE -R $fasta -s $STAT_FILE -p $PROB_FILE --coverage 15 --ipfIterations 0 \
            -1 ${OUT}ReSeq_CYP2C19_37_00${sub}_left.fq -2 ${OUT}ReSeq_CYP2C19_37_00${sub}_right.fq 

        reseq illuminaPE -R $fasta -s $STAT_FILE -p $PROB_FILE --coverage 15 --ipfIterations 0 \
            -1 ${OUT}ReSeq_CYP2C19_37_00${sub}_left2.fq -2 ${OUT}ReSeq_CYP2C19_37_00${sub}_right2.fq 
    done
        
else
    # Set paths for VCF file preparation
    VCF_FILES=${DATA}CYP2C19-6.2.3/GRCh38/
    PYTHON_CODES=${STAGING}Simulations_CYP2C19/Python-codes/
    
    # Remove unnecessary contig lines and add genotype columns
    python ${PYTHON_CODES}vcf_for_reseq.py ${VCF_FILES}CYP2C19_${ALLELE}.vcf

    # Shift position in VCF file, beginning at extended CYP2C locus
    python ${PYTHON_CODES}vcf_shift.py "CYP2C19_${ALLELE}_reseq.vcf" "94400000"
    rm -f CYP2C19_${ALLELE}_reseq.vcf
    mv CYP2C19_${ALLELE}_reseq_shift.vcf ${STAGING}Simulations_CYP2C19/ReSeq_CYP2C19_15x_vcfs/

    VCF=${STAGING}Simulations_CYP2C19/ReSeq_CYP2C19_15x_vcfs/CYP2C19_${ALLELE}_reseq.vcf     

    # Simulate reads with ReSeq
    reseq illuminaPE -R $REF -s $STAT_FILE -p $PROB_FILE -V $VCF --coverage 15 --ipfIterations 0 \
        -1 ${OUT}ReSeq_CYP2C19_${ALLELE}_left.fq -2 ${OUT}ReSeq_CYP2C19_${ALLELE}_right.fq 
    
    reseq illuminaPE -R $REF -s $STAT_FILE -p $PROB_FILE -V $VCF --coverage 15 --ipfIterations 0 \
        -1 ${OUT}ReSeq_CYP2C19_${ALLELE}_left2.fq -2 ${OUT}ReSeq_CYP2C19_${ALLELE}_right2.fq 
fi
conda deactivate

#################
### ALIGNMENT ###
#################

PYTHON_ADJUST=${STAGING}Simulations_CYP2C19/Python-codes/change_bam_qnames.py

if [ ${ALLELE} == 36 ]
then
    # Alignment
    source activate bwa-mem
    for sub in 1 2
    do
        # Set paths
        READS1_V1=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_left.fq; READS1_V2=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_left2.fq
        READS2_V1=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_right.fq; READS2_V2=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_right2.fq
        OUT_ALIGN_V1=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned; OUT_ALIGN_V2=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned2
        
        bwa mem -t 4 -R "@RG\tID:CYP2C19_sim\tSM:CYP2C19_sim" $REF $READS1_V1 $READS2_V1 > ${OUT_ALIGN_V1}.sam
        bwa mem -t 4 -R "@RG\tID:CYP2C19_sim\tSM:CYP2C19_sim" $REF $READS1_V2 $READS2_V2 > ${OUT_ALIGN_V2}.sam
    done 
    conda deactivate

    # Convert to BAM
    source activate samtools 
    for sub in 1 2
    do 
        OUT_ALIGN_V1=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned; OUT_ALIGN_V2=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned2
        samtools view -h -o ${OUT_ALIGN_V1}.bam ${OUT_ALIGN_V1}.sam
        samtools view -h -o ${OUT_ALIGN_V2}.bam ${OUT_ALIGN_V2}.sam
    done 
    conda deactivate

    # Adjust QNAME's in BAM files
    source activate pysam
    for sub in 1 2
    do
        OUT_ALIGN_V1=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned; OUT_ALIGN_V2=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned2
        NEW_ALIGN_V1=ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned_mod; NEW_ALIGN_V2=ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned2_mod
        python $PYTHON_ADJUST ${OUT_ALIGN_V1}.bam ":${ALLELE}_00${sub}_V1"
        python $PYTHON_ADJUST ${OUT_ALIGN_V2}.bam ":${ALLELE}_00${sub}_V2"

        mv -f ${NEW_ALIGN_V1}.bam ${OUT_ALIGN_V1}.bam
        mv -f ${NEW_ALIGN_V2}.bam ${OUT_ALIGN_V2}.bam
    done
    conda deactivate

    # Index the adjusted BAM files
    source activate samtools
    for sub in 1 2
    do
        OUT_ALIGN_V1=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned; OUT_ALIGN_V2=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned2
        samtools sort ${OUT_ALIGN_V1}.bam -o ${OUT_ALIGN_V1}.bam
        samtools index ${OUT_ALIGN_V1}.bam
        rm -f ${OUT_ALIGN_V1}.sam

        samtools sort ${OUT_ALIGN_V2}.bam -o ${OUT_ALIGN_V2}.bam
        samtools index ${OUT_ALIGN_V2}.bam
        rm -f ${OUT_ALIGN_V2}.sam
    done
    conda deactivate

elif [ ${ALLELE} == 37 ]
then
    # Alignment
    source activate bwa-mem
    for sub in {1..5}
    do
        # Set paths
        READS1_V1=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_left.fq; READS1_V2=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_left2.fq
        READS2_V1=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_right.fq; READS2_V2=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_right2.fq
        OUT_ALIGN_V1=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned; OUT_ALIGN_V2=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned2
        
        bwa mem -t 4 -R "@RG\tID:CYP2C19_sim\tSM:CYP2C19_sim" $REF $READS1_V1 $READS2_V1 > ${OUT_ALIGN_V1}.sam
        bwa mem -t 4 -R "@RG\tID:CYP2C19_sim\tSM:CYP2C19_sim" $REF $READS1_V2 $READS2_V2 > ${OUT_ALIGN_V2}.sam
    done 
    conda deactivate

    # Convert to BAM
    source activate samtools 
    for sub in {1..5}
    do 
        OUT_ALIGN_V1=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned; OUT_ALIGN_V2=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned2
        samtools view -h -o ${OUT_ALIGN_V1}.bam ${OUT_ALIGN_V1}.sam
        samtools view -h -o ${OUT_ALIGN_V2}.bam ${OUT_ALIGN_V2}.sam
    done 
    conda deactivate

    # Adjust QNAME's in BAM files
    source activate pysam
    for sub in {1..5}
    do
        OUT_ALIGN_V1=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned; OUT_ALIGN_V2=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned2
        NEW_ALIGN_V1=ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned_mod; NEW_ALIGN_V2=ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned2_mod
        python $PYTHON_ADJUST ${OUT_ALIGN_V1}.bam ":${ALLELE}_00${sub}_V1"
        python $PYTHON_ADJUST ${OUT_ALIGN_V2}.bam ":${ALLELE}_00${sub}_V2"

        mv -f ${NEW_ALIGN_V1}.bam ${OUT_ALIGN_V1}.bam
        mv -f ${NEW_ALIGN_V2}.bam ${OUT_ALIGN_V2}.bam
    done
    conda deactivate

    # Index the adjusted BAM files
    source activate samtools
    for sub in {1..5}
    do
        OUT_ALIGN_V1=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned; OUT_ALIGN_V2=${OUT}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned2
        samtools sort ${OUT_ALIGN_V1}.bam -o ${OUT_ALIGN_V1}.bam
        samtools index ${OUT_ALIGN_V1}.bam
        rm -f ${OUT_ALIGN_V1}.sam

        samtools sort ${OUT_ALIGN_V2}.bam -o ${OUT_ALIGN_V2}.bam
        samtools index ${OUT_ALIGN_V2}.bam
        rm -f ${OUT_ALIGN_V2}.sam
    done
    conda deactivate

else
    # Set paths
    READS1_V1=${OUT}ReSeq_CYP2C19_${ALLELE}_left.fq; READS1_V2=${OUT}ReSeq_CYP2C19_${ALLELE}_left2.fq
    READS2_V1=${OUT}ReSeq_CYP2C19_${ALLELE}_right.fq; READS2_V2=${OUT}ReSeq_CYP2C19_${ALLELE}_right2.fq
    OUT_ALIGN_V1=${OUT}ReSeq_CYP2C19_${ALLELE}_aligned; OUT_ALIGN_V2=${OUT}ReSeq_CYP2C19_${ALLELE}_aligned2

    # Alignment
    source activate bwa-mem
    bwa mem -t 4 -R "@RG\tID:CYP2C19_sim\tSM:CYP2C19_sim" $REF $READS1_V1 $READS2_V1 > ${OUT_ALIGN_V1}.sam
    bwa mem -t 4 -R "@RG\tID:CYP2C19_sim\tSM:CYP2C19_sim" $REF $READS1_V2 $READS2_V2 > ${OUT_ALIGN_V2}.sam
    conda deactivate

    # Convert to BAM
    source activate samtools 
    samtools view -h -o ${OUT_ALIGN_V1}.bam ${OUT_ALIGN_V1}.sam
    samtools view -h -o ${OUT_ALIGN_V2}.bam ${OUT_ALIGN_V2}.sam 
    conda deactivate

    # Adjust QNAME's in BAM files
    source activate pysam
    NEW_ALIGN_V1=ReSeq_CYP2C19_${ALLELE}_aligned_mod; NEW_ALIGN_V2=ReSeq_CYP2C19_${ALLELE}_aligned2_mod
    python $PYTHON_ADJUST ${OUT_ALIGN_V1}.bam ":${ALLELE}_V1"
    python $PYTHON_ADJUST ${OUT_ALIGN_V2}.bam ":${ALLELE}_V2"

    mv -f ${NEW_ALIGN_V1}.bam ${OUT_ALIGN_V1}.bam
    mv -f ${NEW_ALIGN_V2}.bam ${OUT_ALIGN_V2}.bam
    conda deactivate

    # Index the adjusted BAM files
    source activate samtools
    samtools sort ${OUT_ALIGN_V1}.bam -o ${OUT_ALIGN_V1}.bam
    samtools index ${OUT_ALIGN_V1}.bam
    rm -f ${OUT_ALIGN_V1}.sam

    samtools sort ${OUT_ALIGN_V2}.bam -o ${OUT_ALIGN_V2}.bam
    samtools index ${OUT_ALIGN_V2}.bam
    rm -f ${OUT_ALIGN_V2}.sam
    conda deactivate

fi
