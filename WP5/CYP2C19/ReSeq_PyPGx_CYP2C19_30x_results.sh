#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-PyPGx-CYP2C9-30x-results # Name of job
#SBATCH --output=stdout-ReSeq-PyPGx-CYP2C19-30x-results   # Standard output name
#SBATCH --error=stderror-ReSeq-PyPGx-CYP2C19-30x-results    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=8     # Number of CPU cores
#SBATCH --time=0-10:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"
export PYPGX_BUNDLE="${STAGING}Diplotype_Callers/PyPGx/pypgx-bundle"

# Set paths and parameters
DATA=/data/leuven/373/vsc37367/
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
RESULTS=${STAGING}Diplotype_Callers/PyPGx/ReSeq_PyPGx_CYP2C19_30x_calls/
OUT=${STAGING}Diplotype_Callers/PyPGx/
BAMFILES=${STAGING}ReSeq_CYP2C19_diplotypes_30x/
SUMMARY=${OUT}ReSeq_PyPGx_CYP2C19_30x_results.tsv

# Initiate results file
echo -e -n "True-Genotype\tCalled-Genotype\tCorrect\tMissed-Alleles\tCorrect-Alleles" > ${SUMMARY}
echo -e "\tCandidate-Alleles\tAll-Alleles-Cand\tAlleles-Not-Cand\tAlt-Phase-Alleles\tAll-Alleles-Alt-Phase\tAlleles-Not-Alt-Phase" >> ${SUMMARY}

source activate pypgx

for diplo_bam in ${BAMFILES}CYP2C19_diplo_*.bam; do

    # Get filename without bam extension
    diplo_name="$( basename ${diplo_bam} .bam )"

    # Get true genotype
    geno_temp="$( echo ${diplo_name} | grep -oE "([0-9]{1,2}_00[0-9]|[0-9]{1,2})_([0-9]{1,2}_00[0-9]|[0-9]{1,2})" )"

    # Find first allele 
    ALLELE1=$( echo $geno_temp | grep -oE "^[0-9]{1,2}" )
    #if [[ $( echo $geno_temp | grep -oE "^[0-9]{2}_00[0-9]{1}" ) ]]; then
    #    ALLELE1=$( echo $geno_temp | grep -oE "^[0-9]{2}_00[0-9]{1}" )
    #else
    #    ALLELE1=$( echo $geno_temp | grep -oE "^[0-9]{1,2}" )
    #fi

    # Find second allele
    if [[ $( echo $geno_temp | grep -oE "[0-9]{2}_00[0-9]{1}$" ) ]]; then
        ALLELE2=$( echo $geno_temp | grep -oE "[0-9]{2}_00[0-9]{1}$" | grep -oE "^[0-9]{2}" )
    else
        ALLELE2=$( echo $geno_temp | grep -oE "[0-9]{1,2}$" )
    fi
    GENO_TRUE="*${ALLELE1}/*${ALLELE2}"

    # Obtain PyPGx results
    pypgx print-data ${RESULTS}results_${diplo_name}.zip > ${OUT}PyPGx_temp_results_${diplo_name}.txt

    # Extract relevant information
    grep "CYP2C19_sim" ${OUT}PyPGx_temp_results_${diplo_name}.txt | cut -f 2,4,5,6 > ${OUT}PyPGx_temp_results_${diplo_name}_cut.txt
    GENO_CALLED="$( cut -f1 ${OUT}PyPGx_temp_results_${diplo_name}_cut.txt )"

    haplo1="$( cut -f2 ${OUT}PyPGx_temp_results_${diplo_name}_cut.txt )"
    haplo2="$( cut -f3 ${OUT}PyPGx_temp_results_${diplo_name}_cut.txt )"
    CANDIDATE_ALLELES="${haplo1}${haplo2}"

    alt_phase="$( cut -f4 ${OUT}PyPGx_temp_results_${diplo_name}_cut.txt )"

    if [[ "$GENO_TRUE" == "$GENO_CALLED" ]]
    then 
        CORRECT="100%"
        CORRECT_ALLELES="*${ALLELE1};*${ALLELE2}"
        ALL_ALLELES_FOUND="N/A"
        ALLELES_NOT_FOUND=""
        MISSED_ALLELES=""
        ALT_PHASE="N/A"
        ALLELES_NOT_ALT_PHASE=""
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
            ALT_PHASE="N/A"
        else
            ALL_ALLELES_FOUND="No"

            # Verify which alleles were not found
            if [[ "$CANDIDATE_ALLELES" != *"*${ALLELE1};"* ]] && [[ "$CANDIDATE_ALLELES" != *"*${ALLELE2};"* ]]
            then
                ALLELES_NOT_FOUND="*${ALLELE1};*${ALLELE2}"
                
                # Verify whether the non-candidate alleles were suggested as an alternative phase:
                if [[ $alt_phase == *"*${ALLELE1};"* ]] && [[ $alt_phase == *"*${ALLELE2};"* ]]; then
                    ALT_PHASE="Yes"
                    ALLELES_NOT_ALT_PHASE=""
                elif [[ $alt_phase != *"*${ALLELE1};"* ]] && [[ $alt_phase != *"*${ALLELE2};"* ]]; then
                    ALT_PHASE="No"
                    ALLELES_NOT_ALT_PHASE="*${ALLELE1};*${ALLELE2}"
                elif [[ "$alt_phase" == *"*${ALLELE1};"* ]]; then
                    ALT_PHASE="No"
                    ALLELES_NOT_ALT_PHASE="*${ALLELE2}"
                else
                    ALT_PHASE="No"
                    ALLELES_NOT_ALT_PHASE="*${ALLELE1}"
                fi

            elif [[ "$CANDIDATE_ALLELES" == *"*${ALLELE1};"* ]]
            then
                ALLELES_NOT_FOUND="*${ALLELE2}"

                # Verify whether the non-candidate allele was suggested as an alternative phase:
                if [[ $alt_phase == *"*${ALLELE2};"* ]]; then
                    ALT_PHASE="Yes"
                    ALLELES_NOT_ALT_PHASE=""
                else
                    ALT_PHASE="No"
                    ALLELES_NOT_ALT_PHASE="*${ALLELE2}"
                fi
            else
                ALLELES_NOT_FOUND="*${ALLELE1}"

                # Verify whether the non-candidate allele was suggested as an alternative phase:
                if [[ $alt_phase == *"*${ALLELE1};"* ]]; then
                    ALT_PHASE="Yes"
                    ALLELES_NOT_ALT_PHASE=""
                else
                    ALT_PHASE="No"
                    ALLELES_NOT_ALT_PHASE="*${ALLELE1}"
                fi
            fi
        fi
    fi

    # Write to output file

    echo -e "${GENO_TRUE}\t${GENO_CALLED}\t${CORRECT}\t${MISSED_ALLELES}\t${CORRECT_ALLELES}\t${CANDIDATE_ALLELES}\t${ALL_ALLELES_FOUND}\t${ALLELES_NOT_FOUND}\t${alt_phase}\t${ALT_PHASE}\t${ALLELES_NOT_ALT_PHASE}" >> ${SUMMARY}

    # File management
    rm -f ${OUT}PyPGx_temp_results_${diplo_name}*.txt
    
done
