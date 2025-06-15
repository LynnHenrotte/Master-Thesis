#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=ReSeq-Aldy-CYP2C19-30x-results # Name of job
#SBATCH --output=stdout-ReSeq-Aldy-CYP2C19-30x-results   # Standard output name
#SBATCH --error=stderror-ReSeq-Aldy-CYP2C19-30x-results    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=8     # Number of CPU cores
#SBATCH --time=0-10:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=8gb            # Amount of memory (units: gb, mb, kb)
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

# Set paths and parameters
STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
RESULTS=${STAGING}Diplotype_Callers/Aldy/Aldy_ReSeq_CYP2C19_diplotype_calls_30x/
SUMMARY=${STAGING}Diplotype_Callers/Aldy/Aldy_Calls_CYP2C19_30x_summary.tsv
PYTHON_FUNC=${STAGING}Simulations_CYP2C19/Python-codes/aldy_verify_diplotype2.py

# Initiate results file
echo -e "True-diplotype\tCalled-diplotype\tCorrect\tCorrect-alleles\tMissed-alleles" > ${SUMMARY}

source activate python

for diplo_bam in ${BAMFILES}CYP2C19_diplo_*.bam; do

    # Get filename without bam extension
    diplo_name="$( basename ${diplo_bam} .bam )"

    # Get true genotype
    geno_temp="$( echo ${diplo_name} | grep -oE "([0-9]{1,2}_00[0-9]|[0-9]{1,2})_([0-9]{1,2}_00[0-9]|[0-9]{1,2})" )"
    
    # Find first allele 
    allele1=$( echo $geno_temp | grep -oE "^[0-9]{1,2}" )

    # Find second allele
    if [[ $( echo $geno_temp | grep -oE "[0-9]{2}_00[0-9]{1}$" ) ]]; then
        allele2=$( echo $geno_temp | grep -oE "[0-9]{2}_00[0-9]{1}$" | grep -oE "^[0-9]{2}" )
    else
        allele2=$( echo $geno_temp | grep -oE "[0-9]{1,2}$" )
    fi
    true_dip="*${allele1}/*${allele2}"
    true_dip2="*${allele2}/*${allele1}"
    echo "$geno_temp: $allele1; $allele2; -> $true_dip"

    # Extract called diplotype using Python script
    aldy_output="${RESULTS}${diplo_name}.aldy"
    aldy_dip_line=$( grep "#Solution 1:" $aldy_output )
    result=$( python $PYTHON_FUNC "$aldy_dip_line" "$allele1;$allele2" )
    IFS=","; read -r true_dip_temp called_dip percent <<< $result

    # Verify how correct the diplotype was called
    if [[ "$true_dip" == "$called_dip" ]] || [[ "$true_dip2" == "$called_dip" ]]
    then 
        correct="100%"
        correct_alleles="*${allele1};*${allele2}"
        missed_alleles=""

    else
        
        if [[ *"$called_dip"* == *"*${allele1}/"* ]] || [[ *"$called_dip"* == *"/*${allele1}"* ]]
        then
            correct="50%"
            correct_alleles="*${allele1}"
            missed_alleles="*${allele2}"

        elif [[ "$called_dip" == *"*${allele2}/"* ]] || [[ "$called_dip" == *"/*${allele2}"* ]]
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
