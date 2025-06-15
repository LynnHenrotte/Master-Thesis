#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-Aldy-CYP2C9-10x-results # Name of job
#SBATCH --output=stdout-ReSeq-Aldy-CYP2C9-10x-results   # Standard output name
#SBATCH --error=stderror-ReSeq-Aldy-CYP2C9-10x-results    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=8     # Number of CPU cores
#SBATCH --time=0-01:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

# Set paths and parameters
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
RESULTS=${STAGING}Diplotype_Callers/Aldy/Aldy_ReSeq_CYP2C9_diplotype_calls_10x/
SUMMARY=${STAGING}Diplotype_Callers/Aldy/Aldy_Calls_CYP2C9_10x_summary.tsv
PYTHON_FUNC=${STAGING}Simulations_CYP2C19/Python-codes/aldy_verify_diplotype.py

# Initiate results file
echo -e "True-diplotype\tCalled-diplotype\tCorrect\tCorrect-alleles\tMissed-alleles" > ${SUMMARY}

source activate python

for ((allele1=1; allele1<=85; allele1++))
do for ((allele2=allele1; allele2<=85; allele2++))
    do

    # Extract called diplotype using Python script
    aldy_output="${RESULTS}CYP2C9_diplo_${allele1}_${allele2}.aldy"
    aldy_dip_line=$( grep "#Solution 1:" $aldy_output )
    result=$( python $PYTHON_FUNC "$aldy_dip_line" "CYP2C9_diplo_${allele1}_${allele2}" )
    IFS=","; read -r true_dip called_dip percent <<< $result

    # True diplotype
    true_dip="*${allele1}/*${allele2}"
    true_dip2="*${allele2}/*${allele1}"

    # Verify how correct the diplotype was called
    if [[ "$true_dip" == "$called_dip" ]] || [[ "$true_dip2" == "$called_dip" ]]
    then 
        correct="100%"
        correct_alleles="*${allele1};*${allele2}"
        missed_alleles=""

    else
        
        if [[ "$called_dip" == "*${allele1}/"* ]] || [[ "$called_dip" == *"/*${allele1}" ]]
        then
            correct="50%"
            correct_alleles="*${allele1}"
            missed_alleles="*${allele2}"

        elif [[ "$called_dip" == "*${allele2}/"* ]] || [[ "$called_dip" == *"/*${allele2}" ]]
        then
            correct="50%"
            correct_alleles="*${allele2}"
            missed_alleles="*${allele1}"
        else
            correct="0%"
            correct_alleles=""
            missed_alleles="*${allele1};*${allele2}"
        fi
    fi

    # Write to summary file
    echo -e "${true_dip}\t${called_dip}\t${correct}\t${correct_alleles}\t${missed_alleles}" >> ${SUMMARY}

    done
done
