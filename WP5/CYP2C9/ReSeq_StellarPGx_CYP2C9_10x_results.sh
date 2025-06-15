#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-StellarPGx-CYP2C9-10x-results # Name of job
#SBATCH --output=stdout-ReSeq-StellarPGx-CYP2C9-10x-results   # Standard output name
#SBATCH --error=stderror-ReSeq-StellarPGx-CYP2C9-10x-results    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=8     # Number of CPU cores
#SBATCH --time=0-01:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

# Set paths and parameters
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
RESULTS=${STAGING}Diplotype_Callers/StellarPGx/ReSeq_CYP2C9_diplotype_calls_10x/cyp2c9/alleles/
SUMMARY=${STAGING}Diplotype_Callers/StellarPGx/ReSeq_StellarPGx_CYP2C9_10x_summary.tsv

# Initiate results file
echo -e "True-diplotype\tCalled-diplotype\tCorrect\tCorrect-alleles\tMissed-alleles" > ${SUMMARY}

for allele1 in $(seq 1 85)
do
    for allele2 in $(seq ${allele1} 85)
    do
    
    # True diplotype
    true_dip="*${allele1}/*${allele2}"
    true_dip2="*${allele2}/*${allele1}"

    # Extract called diplotype
    if [ -f ${RESULTS}CYP2C9_diplo_${allele1}_${allele2}_cyp2c9.alleles ]
    then 
        called_dip=$( grep "*" ${RESULTS}CYP2C9_diplo_${allele1}_${allele2}_cyp2c9.alleles | sed 's/[][]//g' )
    else
        called_dip=""
    fi

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

