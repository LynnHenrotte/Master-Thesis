#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-PyPGx-CYP2C9-10x-results # Name of job
#SBATCH --output=stdout-ReSeq-PyPGx-CYP2C9-10x-results   # Standard output name
#SBATCH --error=stderror-ReSeq-PyPGx-CYP2C9-10x-results    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=8     # Number of CPU cores
#SBATCH --time=0-05:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"
export PYPGX_BUNDLE="${STAGING}Diplotype_Callers/PyPGx/pypgx-bundle"

# Set paths and parameters
DATA=/data/leuven/373/vsc37367/
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
RESULTS=${STAGING}Diplotype_Callers/PyPGx/ReSeq_PyPGx_CYP2C9_10x_calls/
SUMMARY=${OUT}ReSeq_PyPGx_CYP2C9_10x_results.tsv

# Initiate results file
echo -e "True-Genotype\tCalled-Genotype\tCorrect\tMissed-Alleles\tCorrect-Alleles\tCandidate-Alleles\tAll-Alleles-Found\tAlleles-Not-Found" > ${SUMMARY}

source activate pypgx

for ((ALLELE1=1; ALLELE1<=71; ALLELE1++)); do 
    for ((ALLELE2=ALLELE1; ALLELE2<=71; ALLELE2++)); do

    # Obtain PyPGx results
    pypgx print-data ${RESULTS}results_CYP2C9_${ALLELE1}_${ALLELE2}.zip > ${OUT}PyPGx_temp_results_CYP2C9_${ALLELE1}_${ALLELE2}.txt

    # Extract relevant information
    grep "CYP2C9_sim" ${OUT}PyPGx_temp_results_CYP2C9_${ALLELE1}_${ALLELE2}.txt | cut -f 2,4,5 > ${OUT}PyPGx_temp_results_CYP2C9_${ALLELE1}_${ALLELE2}_cut.txt
    GENO_CALLED="$( cut -f 1 ${OUT}PyPGx_temp_results_CYP2C9_${ALLELE1}_${ALLELE2}_cut.txt )"
    GENO_TRUE="*${ALLELE1}/*${ALLELE2}"
    haplo1="$( cut -f 2 ${OUT}PyPGx_temp_results_CYP2C9_${ALLELE1}_${ALLELE2}_cut.txt )"
    haplo2="$( cut -f 3 ${OUT}PyPGx_temp_results_CYP2C9_${ALLELE1}_${ALLELE2}_cut.txt )"
    CANDIDATE_ALLELES="${haplo1}${haplo2}"

    if [[ "$GENO_TRUE" == "$GENO_CALLED" ]]
    then 
        CORRECT="100%"
        CORRECT_ALLELES="*${ALLELE1};*${ALLELE2}"
        ALL_ALLELES_FOUND="Yes"
        ALLELES_NOT_FOUND=""
        MISSED_ALLELES=""
    else
        # Verify whether one or none of the alleles was correctly called
        if [[ "$GENO_CALLED" == "*${ALLELE1}/"* ]] || [[ "$GENO_CALLED" == *"/*${ALLELE1}" ]]
        then
            CORRECT="50%"
            CORRECT_ALLELES="*${ALLELE1}"
            MISSED_ALLELES="*${ALLELE2}"
        elif [[ "$GENO_CALLED" == "*${ALLELE2}/"* ]] || [[ "$GENO_CALLED" == *"/*${ALLELE2}" ]]
        then
            CORRECT="50%"
            CORRECT_ALLELES="*${ALLELE2}"
            MISSED_ALLELES="*${ALLELE1}"
        else
            CORRECT="0%"
            CORRECT_ALLELES=""
            MISSED_ALLELES="*${ALLELE1};*${ALLELE2}"
        fi

        # Verify whether all correct alleles were found as candidates 
        if [[ "$CANDIDATE_ALLELES" == *"*${ALLELE1};"* ]] && [[ "$CANDIDATE_ALLELES" == *"*${ALLELE2};"* ]]
        then
            ALL_ALLELES_FOUND="Yes"
            ALLELES_NOT_FOUND=""
        else
            ALL_ALLELES_FOUND="No"

            # Verify which alleles were not found
            if [[ "$CANDIDATE_ALLELES" != *"*${ALLELE1};"* ]] && [[ "$CANDIDATE_ALLELES" != *"*${ALLELE2};"* ]]
            then
                ALLELES_NOT_FOUND="*${ALLELE1};*${ALLELE2}"
            elif [[ "$CANDIDATE_ALLELES" == *"*${ALLELE1};"* ]]
            then
                ALLELES_NOT_FOUND="*${ALLELE2}"
            else
                ALLELES_NOT_FOUND="*${ALLELE1}"
            fi
        fi
    fi

    # Write to output file

    echo -e "${GENO_TRUE}\t${GENO_CALLED}\t${CORRECT}\t${MISSED_ALLELES}\t${CORRECT_ALLELES}\t${CANDIDATE_ALLELES}\t${ALL_ALLELES_FOUND}\t${ALLELES_NOT_FOUND}" >> ${SUMMARY}

    # File management
    rm -f ${OUT}PyPGx_temp_results_CYP2C9_${ALLELE1}_${ALLELE2}*.txt
    
    done
done
