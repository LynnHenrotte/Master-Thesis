#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-CYP2C19-15x-SV-calls_%a  # Name of job
#SBATCH --output=stdout-ReSeq-CYP2C19-15x-SV-calls_%a   # Standard output name
#SBATCH --error=stderror-ReSeq-CYP2C19-15x-SV-calls_%a    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=8     # Number of CPU cores
#SBATCH --time=0-10:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
#SBATCH --array=37
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

ALLELE=$SLURM_ARRAY_TASK_ID
DATA=/data/leuven/373/vsc37367/
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
REF=${DATA}references/chr10.fa
MANTA=${DATA}miniconda3/envs/manta/bin/configManta.py
IN_BAM=${STAGING}Simulations_CYP2C19/ReSeq_CYP2C19_15x/
OUT=${STAGING}Simulations_CYP2C19/ReSeq_CYP2C19_15x_SV_calls/

# Set number of sub alleles
if [ $ALLELE == 36 ]; then
    num_sub=2
elif [ $ALLELE == 37 ]; then
    num_sub=5
else
    exit "This allele does not have structural variation"
fi

# Initiate output file
SUMMARY=${STAGING}Simulations_CYP2C19/ReSeq_CYP2C19_${ALLELE}_15x_SV_calling_results.tsv
echo -e "Allele\tNum-SV\tNum-SV-Called\tNum-TRUE-SV-Called\tAll-SV-Called" > $SUMMARY

for ((sub=1; sub<=$num_sub; sub++)); do
    
    ## Set paths ##
    OUT_DIR1=${STAGING}Simulations_CYP2C19/ReSeq_CYP2C19_15x_SV_calls/ReSeq_CYP2C19_15x_${ALLELE}_00${sub}_SV_call/
    OUT_DIR2=${STAGING}Simulations_CYP2C19/ReSeq_CYP2C19_15x_SV_calls/ReSeq_CYP2C19_15x_${ALLELE}_00${sub}_SV_call2/

    ## Structural Variant Calling ##
    
    source activate manta

    # Step 1: configuration
    $MANTA --bam "${IN_BAM}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned.bam" --referenceFasta ${REF} --runDir ${OUT_DIR1}
    $MANTA --bam "${IN_BAM}ReSeq_CYP2C19_${ALLELE}_00${sub}_aligned2.bam" --referenceFasta ${REF} --runDir ${OUT_DIR2}

    ## Step 2: Execution
    ${OUT_DIR1}runWorkflow.py -j 8
    ${OUT_DIR2}runWorkflow.py -j 8

    conda deactivate

    ## Verification ##
    
    # Get scored SVs found by manta
    scored_SV_results1=${OUT_DIR1}results/variants/diploidSV.vcf
    gzip -k -d ${scored_SV_results1}.gz
    scored_SV_results2=${OUT_DIR2}results/variants/diploidSV.vcf
    gzip -k -d ${scored_SV_results2}.gz

    # Convert diploidSV VCF to BED 
    source activate bcftools

    bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE' $scored_SV_results1 -o ${OUT_DIR1}CYP2C19_${ALLELE}_00${sub}_SV_call.bed
    bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE' $scored_SV_results2 -o ${OUT_DIR2}CYP2C19_${ALLELE}_00${sub}_SV_call2.bed

    conda deactivate

    # Find identical SVs between the true and called BED files
    BED_TRUE=${STAGING}Simulations_CYP2C19/SV_beds/cyp2c19_${ALLELE}_${sub}_VCF.bed
    BED_CALL1=${OUT_DIR1}CYP2C19_${ALLELE}_00${sub}_SV_call.bed
    BED_CALL2=${OUT_DIR2}CYP2C19_${ALLELE}_00${sub}_SV_call2.bed
    comm -1 -2 $BED_TRUE $BED_CALL1 > ${OUT_DIR1}CYP2C19_${ALLELE}_00${sub}_overlap.txt
    comm -1 -2 $BED_TRUE $BED_CALL2 > ${OUT_DIR2}CYP2C19_${ALLELE}_00${sub}_overlap2.txt

    # Count number of true and called SVs + number of overlapping SVs
    NUM_TRUE=$( wc -l $BED_TRUE | cut -f 1 -d' ' )
    NUM_CALLED1=$( wc -l $BED_CALL1 | cut -f 1 -d' ' )
    NUM_CALLED2=$( wc -l $BED_CALL2 | cut -f 1 -d' ' )
    NUM_BOTH1=$( wc -l ${OUT_DIR1}CYP2C19_${ALLELE}_00${sub}_overlap.txt | cut -f 1 -d' ' )
    NUM_BOTH2=$( wc -l ${OUT_DIR2}CYP2C19_${ALLELE}_00${sub}_overlap2.txt | cut -f 1 -d' ' )

    # Verify that all SVs were found
    if [ $NUM_TRUE == $NUM_BOTH1 ]; then
        ALL_CALLED1="Yes"
    else
        ALL_CALLED1="No"
    fi
    if [ $NUM_TRUE == $NUM_BOTH2 ]; then
        ALL_CALLED2="Yes"
    else
        ALL_CALLED2="No"
    fi

    # Write to output file
    echo -e "CYP2C19_${ALLELE}_00${sub}_V1\t${NUM_TRUE}\t${NUM_CALLED1}\t${NUM_BOTH1}\t${ALL_CALLED1}" >> $SUMMARY
    echo -e "CYP2C19_${ALLELE}_00${sub}_V2\t${NUM_TRUE}\t${NUM_CALLED2}\t${NUM_BOTH2}\t${ALL_CALLED2}" >> $SUMMARY

done
