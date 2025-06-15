#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=simulate-ReSeq-CYP2C9_5x_%a  # Name of job
#SBATCH --output=stdout-simulate-ReSeq-CYP2C9_5x_%a   # Standard output name
#SBATCH --error=stderror-simulate-ReSeq-CYP2C9_5x_%a    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=4     # Number of CPU cores
#SBATCH --time=0-05:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
#SBATCH --array=1-85
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

##################
### SIMULATION ###
##################

# Set paths
ALLELE=$SLURM_ARRAY_TASK_ID
DATA=/data/leuven/373/vsc37367/
REF=${DATA}references/cyp2c.fa
STAT_FILE=${DATA}Simulators/ReSeq/Hs-Nova-TruSeq.reseq
PROB_FILE=${DATA}Simulators/ReSeq/Hs-Nova-TruSeq.reseq.ipf

OUT=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_5x/

if [ ${ALLELE} == 1 ]
then
    # Simulate reads with ReSeq
    source activate reseq

    reseq illuminaPE -R $REF -s $STAT_FILE -p $PROB_FILE --coverage 5 --ipfIterations 0 \
        -1 ${OUT}ReSeq_CYP2C9_${ALLELE}_left.fq -2 ${OUT}ReSeq_CYP2C9_${ALLELE}_right.fq

    reseq illuminaPE -R $REF -s $STAT_FILE -p $PROB_FILE --coverage 5 --ipfIterations 0 \
        -1 ${OUT}ReSeq_CYP2C9_${ALLELE}_left2.fq -2 ${OUT}ReSeq_CYP2C9_${ALLELE}_right2.fq

    conda deactivate

else
    # Set path for VCF file preparation
    VCF_FILES=${DATA}CYP2C9-6.2.3/GRCh38
    PYTHON_CODES=${DATA}Simulations_CYP2C9/Python-codes/
    
    # Remove unnecessary contig lines and add genotype columns
    python ${PYTHON_CODES}vcf_for_reseq.py ${VCF_FILES}/CYP2C9_${ALLELE}.vcf

    # Shift position in VCF file, beginning at CYP2C locus
    python ${PYTHON_CODES}vcf_shift.py "CYP2C9_${ALLELE}_reseq.vcf" "94682000"
    rm -f CYP2C9_${ALLELE}_reseq.vcf
    mv CYP2C9_${ALLELE}_reseq_shift.vcf ${DATA}Simulations_CYP2C9/ReSeq/Shifted_ReSeq_vcfs/

    VCF=${DATA}Simulations_CYP2C9/ReSeq/Shifted_ReSeq_vcfs/CYP2C9_${ALLELE}_reseq_shift.vcf     

    # Simulate reads with ReSeq
    source activate reseq

    reseq illuminaPE -R $REF -s $STAT_FILE -p $PROB_FILE -V $VCF --coverage 5 --ipfIterations 0 \
        -1 ${OUT}ReSeq_CYP2C9_${ALLELE}_left.fq -2 ${OUT}ReSeq_CYP2C9_${ALLELE}_right.fq \
        --writeSysError ${OUT}ReSeq_CYP2C9_${ALLELE}_errors.fq
    
    reseq illuminaPE -R $REF -s $STAT_FILE -p $PROB_FILE -V $VCF --coverage 5 --ipfIterations 0 \
        -1 ${OUT}ReSeq_CYP2C9_${ALLELE}_left2.fq -2 ${OUT}ReSeq_CYP2C9_${ALLELE}_right2.fq \
        --writeSysError ${OUT}ReSeq_CYP2C9_${ALLELE}_errors2.fq

    conda deactivate
fi

###################
#### ALIGNMENT ####
###################

# Set paths
OUT=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_5x/
REF=${DATA}references/chr10.fa
READS1_V1=${OUT}ReSeq_CYP2C9_${ALLELE}_left.fq
READS2_V1=${OUT}ReSeq_CYP2C9_${ALLELE}_right.fq
OUT_ALIGN_V1=${OUT}ReSeq_CYP2C9_${ALLELE}_aligned
READS1_V2=${OUT}ReSeq_CYP2C9_${ALLELE}_left2.fq
READS2_V2=${OUT}ReSeq_CYP2C9_${ALLELE}_right2.fq
OUT_ALIGN_V2=${OUT}ReSeq_CYP2C9_${ALLELE}_aligned2

# Alignment
source activate bwa-mem

bwa mem -t 4 -R "@RG\tID:CYP2C9_sim\tSM:CYP2C9_sim" $REF $READS1_V1 $READS2_V1 > ${OUT_ALIGN_V1}.sam
bwa mem -t 4 -R "@RG\tID:CYP2C9_sim\tSM:CYP2C9_sim" $REF $READS1_V2 $READS2_V2 > ${OUT_ALIGN_V2}.sam

conda deactivate

# Convert SAM file to BAM file
source activate samtools 

samtools view -h -o ${OUT_ALIGN_V1}.bam ${OUT_ALIGN_V1}.sam
samtools view -h -o ${OUT_ALIGN_V2}.bam ${OUT_ALIGN_V2}.sam

conda deactivate

# Adjust the bam files' QNAMEs
PYTHON_ADJUST=${DATA}Simulations_CYP2C9/Python-codes/change_bam_qnames.py
BAM_IN_V1=${OUT}ReSeq_CYP2C9_${ALLELE}_aligned.bam
BAM_NEW_V1=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_${ALLELE}_aligned_mod.bam
BAM_IN_V2=${OUT}ReSeq_CYP2C9_${ALLELE}_aligned2.bam
BAM_NEW_V2=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_${ALLELE}_aligned2_mod.bam

source activate pysam

python $PYTHON_ADJUST "$BAM_IN_V1" ":${ALLELE}_V1"
python $PYTHON_ADJUST "$BAM_IN_V2" ":${ALLELE}_V2"

conda deactivate

# Replace original BAM files by adjusted BAM files
mv -f $BAM_NEW_V1 $BAM_IN_V1
mv -f $BAM_NEW_V2 $BAM_IN_V2

# Index the adjusted BAM files
source activate samtools

samtools sort ${OUT_ALIGN_V1}.bam -o ${OUT_ALIGN_V1}.bam
samtools index ${OUT_ALIGN_V1}.bam
rm -f ${OUT_ALIGN_V1}.sam

samtools sort ${OUT_ALIGN_V2}.bam -o ${OUT_ALIGN_V2}.bam
samtools index ${OUT_ALIGN_V2}.bam
rm -f ${OUT_ALIGN_V2}.sam

conda deactivate

#######################
### VARIANT CALLING ###
#######################

if [ ${ALLELE} != 1 ]
then
    OUT_VCF_V1=${OUT}ReSeq_CYP2C9_${ALLELE}_VarCall.vcf
    OUT_VCF_V2=${OUT}ReSeq_CYP2C9_${ALLELE}_VarCall2.vcf
    
    source activate gatk4
    gatk HaplotypeCaller -I ${OUT_ALIGN_V1}.bam -O ${OUT_VCF_V1} -R $REF
    gatk HaplotypeCaller -I ${OUT_ALIGN_V2}.bam -O ${OUT_VCF_V2} -R $REF
    conda deactivate
fi

################
### COVERAGE ###
################

# Specify paths
BED_PATH=${DATA}Simulations_CYP2C9/cyp2c.bed
OUT_COVERAGE=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_5x_coverage/
COVERAGE_SUMMARY=${DATA}Simulations_CYP2C9/ReSeq/ReSeq_CYP2C9_5x_coverage_summary.tsv

source activate mosdepth

mosdepth -t 3 -c chr10 -n -b ${BED_PATH} ${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V1 ${BAM_IN_V1}
mosdepth -t 3 -c chr10 -n -b ${BED_PATH} ${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V2 ${BAM_IN_V2}

SUMMARY_FILE_V1=${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V1.mosdepth.summary.txt
SUMMARY_FILE_V2=${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V2.mosdepth.summary.txt

LINE_V1=$( grep "chr10_region" ${SUMMARY_FILE_V1} )
LINE_V2=$( grep "chr10_region" ${SUMMARY_FILE_V2} )
IFS=$'\t'; read -ra elements_V1 <<< $LINE_V1
IFS=$'\t'; read -ra elements_V2 <<< $LINE_V2

echo -e "${ALLELE}\t${elements_V1[3]}\t${elements_V1[4]}\t${elements_V1[5]}" >> ${COVERAGE_SUMMARY}
echo -e "${ALLELE}\t${elements_V2[3]}\t${elements_V2[4]}\t${elements_V2[5]}" >> ${COVERAGE_SUMMARY}

# File management
grep "chr10" ${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V1.mosdepth.region.dist.txt > ${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V1.mosdepth.dist.txt
rm -f ${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V1.mosdepth.region.dist.txt 
rm -f ${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V1.mosdepth.global.dist.txt
rm -f ${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V1.regions*
rm -f ${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V1.mosdepth.summary.txt

grep "chr10" ${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V2.mosdepth.region.dist.txt > ${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V2.mosdepth.dist.txt
rm -f ${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V2.mosdepth.region.dist.txt 
rm -f ${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V2.mosdepth.global.dist.txt
rm -f ${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V2.regions*
rm -f ${OUT_COVERAGE}ReSeq_CYP2C9_${ALLELE}_V2.mosdepth.summary.txt

conda deactivate
