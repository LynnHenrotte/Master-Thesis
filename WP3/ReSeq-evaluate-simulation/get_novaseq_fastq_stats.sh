#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=get-novaseq-fastq-stats # Name of job
#SBATCH --output=stdout-get-novaseq-fastq-stats   # Standard output name
#SBATCH --error=stderror-get-novaseq-fastq-stats    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=10     # Number of CPU cores
#SBATCH --time=0-10:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=10gb            # Amount of memory (units: gb, mb, kb)
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
REF=/staging/leuven/stg_00156/references/hg38.fa
BAM_FILES=${STAGING}real_data/real_data_BAMs_cut/
FQ_FILES=${STAGING}real_data/real_data_FQs/
OUT=${STAGING}real_data/evaluation_output/

for bamfile in ${BAM_FILES}*.bam; do

    file_name="$( basename ${bamfile} .bam )"

    # Sort BAM file by read name
    source activate samtools

    samtools sort -n ${bamfile} -o ${bamfile}
    
    conda deactivate  

    # Convert the BAM file to FASTQ files.
    if [ ! -f ${FQ_FILES}${file_name}_left.fq ] && [ ! -f ${FQ_FILES}${file_name}_right.fq ]; then

        source activate bedtools

        bedtools bamtofastq -i ${bamfile} -fq ${FQ_FILES}${file_name}_left.fq \
            -fq2 ${FQ_FILES}${file_name}_right.fq

        conda deactivate

    fi

    # Obtain data characteristics
    source activate bbmap

    bbmap.sh ref=$REF in=${FQ_FILES}${file_name}_left.fq in2=${FQ_FILES}${file_name}_right.fq  \
        out=temp.sam statsfile=${OUT}${file_name}_map_stats.txt \
        mhist=${OUT}${file_name}_error_rates_per_base.txt \
        qhist=${OUT}${file_name}_data_average_qscore_per_base.txt \
        gchist=${OUT}${file_name}_GC_content.txt \
        ihist=${OUT}${file_name}_insert_sizes.txt

    rm -f temp.sam

    conda deactivate

done
