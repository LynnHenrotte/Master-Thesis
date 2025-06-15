#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=venn-diagram-data # Name of job
#SBATCH --output=stdout-venn-diagram-data   # Standard output name
#SBATCH --error=stderror-venn-diagram-data   # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=12     # Number of CPU cores
#SBATCH --time=0-10:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=10gb            # Amount of memory (units: gb, mb, kb)
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/
ALDY=${STAGING}Diplotype_Callers/Aldy/
PYPGX=${STAGING}Diplotype_Callers/PyPGx/
STELLARPGX=${STAGING}Diplotype_Callers/StellarPGx/
OUT=${STAGING}Diplotype_Callers/Venn_diagram_data/

#### Aldy ####

### CYP2C9 ###

## 60x ##

# Perfect calls
cut -f 1 ${ALDY}Aldy_Calls_CYP2C9_60x_100%.tsv | grep -e "^*" > ${OUT}Aldy_CYP2C9_60x_100%.txt

# Partially correct calls
cut -f 1 ${ALDY}Aldy_Calls_CYP2C9_60x_50%.tsv | grep -e "^*" > ${OUT}Aldy_CYP2C9_60x_50%.txt

# Wrong calls
cut -f 1 ${ALDY}Aldy_Calls_CYP2C9_60x_0%.tsv | grep -e "^*" > ${OUT}Aldy_CYP2C9_60x_0%.txt

## 30x ##

# Perfect calls
cut -f 1 ${ALDY}Aldy_Calls_CYP2C9_30x_100%.tsv | grep -e "^*" > ${OUT}Aldy_CYP2C9_30x_100%.txt

# Partially correct calls
cut -f 1 ${ALDY}Aldy_Calls_CYP2C9_30x_50%.tsv | grep -e "^*" > ${OUT}Aldy_CYP2C9_30x_50%.txt

# Wrong calls
cut -f 1 ${ALDY}Aldy_Calls_CYP2C9_30x_0%.tsv | grep -e "^*" > ${OUT}Aldy_CYP2C9_30x_0%.txt

## 10x ##

# Perfect calls
cut -f 1 ${ALDY}Aldy_Calls_CYP2C9_10x_100%.tsv | grep -e "^*" > ${OUT}Aldy_CYP2C9_10x_100%.txt

# Partially correct calls
cut -f 1 ${ALDY}Aldy_Calls_CYP2C9_10x_50%.tsv | grep -e "^*" > ${OUT}Aldy_CYP2C9_10x_50%.txt

# Wrong calls
cut -f 1 ${ALDY}Aldy_Calls_CYP2C9_10x_0%.tsv | grep -e "^*" > ${OUT}Aldy_CYP2C9_10x_0%.txt

### CYP2C19 ###

## 30x ##

# Perfect calls
cut -f 1 ${ALDY}Aldy_Calls_CYP2C19_30x_100%.tsv | grep -e "^*" > ${OUT}Aldy_CYP2C19_30x_100%.txt

# Partially correct calls
cut -f 1 ${ALDY}Aldy_Calls_CYP2C19_30x_50%.tsv | grep -e "^*" > ${OUT}Aldy_CYP2C19_30x_50%.txt

# Wrong calls
cut -f 1 ${ALDY}Aldy_Calls_CYP2C19_30x_0%.tsv | grep -e "^*" > ${OUT}Aldy_CYP2C19_30x_0%.txt

#### PyPGx ####

### CYP2C9 ###

## 60x ##

# Perfect calls
cut -f 1 ${PYPGX}ReSeq_PyPGx_results_CYP2C9_60x_100%.tsv | grep -e "^*" > ${OUT}PyPGx_CYP2C9_60x_100%.txt

# Partially correct calls
cut -f 1 ${PYPGX}ReSeq_PyPGx_results_CYP2C9_60x_50%.tsv | grep -e "^*" > ${OUT}PyPGx_CYP2C9_60x_50%.txt

# Wrong calls
cut -f 1 ${PYPGX}ReSeq_PyPGx_results_CYP2C9_60x_0%.tsv | grep -e "^*" > ${OUT}PyPGx_CYP2C9_60x_0%.txt

## 30x ##

# Perfect calls
cut -f 1 ${PYPGX}ReSeq_PyPGx_results_CYP2C9_30x_100%.tsv | grep -e "^*" > ${OUT}PyPGx_CYP2C9_30x_100%.txt

# Partially correct calls
cut -f 1 ${PYPGX}ReSeq_PyPGx_results_CYP2C9_30x_50%.tsv | grep -e "^*" > ${OUT}PyPGx_CYP2C9_30x_50%.txt

# Wrong calls
cut -f 1 ${PYPGX}ReSeq_PyPGx_results_CYP2C9_30x_0%.tsv | grep -e "^*" > ${OUT}PyPGx_CYP2C9_30x_0%.txt

## 10x ##

# Perfect calls
cut -f 1 ${PYPGX}ReSeq_PyPGx_results_CYP2C9_10x_100%.tsv | grep -e "^*" > ${OUT}PyPGx_CYP2C9_10x_100%.txt

# Partially correct calls
cut -f 1 ${PYPGX}ReSeq_PyPGx_results_CYP2C9_10x_50%.tsv | grep -e "^*" > ${OUT}PyPGx_CYP2C9_10x_50%.txt

# Wrong calls
cut -f 1 ${PYPGX}ReSeq_PyPGx_results_CYP2C9_10x_0%.tsv | grep -e "^*" > ${OUT}PyPGx_CYP2C9_10x_0%.txt

### CYP2C19 ###

## 30x ##

# Perfect calls
cut -f 1 ${PYPGX}ReSeq_PyPGx_CYP2C19_30x_results_100%.tsv | grep -e "^*" > ${OUT}PyPGx_CYP2C19_30x_100%.txt

# Partially correct calls
cut -f 1 ${PYPGX}ReSeq_PyPGx_CYP2C19_30x_results_50%.tsv | grep -e "^*" > ${OUT}PyPGx_CYP2C19_30x_50%.txt

# Wrong calls
cut -f 1 ${PYPGX}ReSeq_PyPGx_CYP2C19_30x_results_0%.tsv | grep -e "^*" > ${OUT}PyPGx_CYP2C19_30x_0%.txt

#### StellarPGx ####

### CYP2C9 ###

## 60x ##

# Perfect calls
cut -f 1 ${STELLARPGX}ReSeq_StellarPGx_CYP2C9_60x_summary_100%.tsv | grep -e "^*" > ${OUT}StellarPGx_CYP2C9_60x_100%.txt

# Partially correct calls
cut -f 1 ${STELLARPGX}ReSeq_StellarPGx_CYP2C9_60x_summary_50%.tsv | grep -e "^*" > ${OUT}StellarPGx_CYP2C9_60x_50%.txt

# Wrong calls
cut -f 1 ${STELLARPGX}ReSeq_StellarPGx_CYP2C9_60x_summary_0%.tsv | grep -e "^*" > ${OUT}StellarPGx_CYP2C9_60x_0%.txt

## 30x ##

# Perfect calls
cut -f 1 ${STELLARPGX}ReSeq_StellarPGx_CYP2C9_30x_summary_100%.tsv | grep -e "^*" > ${OUT}StellarPGx_CYP2C9_30x_100%.txt

# Partially correct calls
cut -f 1 ${STELLARPGX}ReSeq_StellarPGx_CYP2C9_30x_summary_50%.tsv | grep -e "^*" > ${OUT}StellarPGx_CYP2C9_30x_50%.txt

# Wrong calls
cut -f 1 ${STELLARPGX}ReSeq_StellarPGx_CYP2C9_30x_summary_0%.tsv | grep -e "^*" > ${OUT}StellarPGx_CYP2C9_30x_0%.txt

## 10x ##

# Perfect calls
cut -f 1 ${STELLARPGX}ReSeq_StellarPGx_CYP2C9_10x_summary_100%.tsv | grep -e "^*" > ${OUT}StellarPGx_CYP2C9_10x_100%.txt

# Partially correct calls
cut -f 1 ${STELLARPGX}ReSeq_StellarPGx_CYP2C9_10x_summary_50%.tsv | grep -e "^*" > ${OUT}StellarPGx_CYP2C9_10x_50%.txt

# Wrong calls
cut -f 1 ${STELLARPGX}ReSeq_StellarPGx_CYP2C9_10x_summary_0%.tsv | grep -e "^*" > ${OUT}StellarPGx_CYP2C9_10x_0%.txt

### CYP2C19 ###

# Perfect calls
cut -f 1 ${STELLARPGX}ReSeq_StellarPGx_CYP2C19_30x_results_100%_V2.tsv | grep -e "^*" > ${OUT}StellarPGx_CYP2C19_30x_100%.txt

# Partially correct calls
cut -f 1 ${STELLARPGX}ReSeq_StellarPGx_CYP2C19_30x_results_50%_V2.tsv | grep -e "^*" > ${OUT}StellarPGx_CYP2C19_30x_50%.txt

# Wrong calls
cut -f 1 ${STELLARPGX}ReSeq_StellarPGx_CYP2C19_30x_results_0%_V2.tsv | grep -e "^*" > ${OUT}StellarPGx_CYP2C19_30x_0%.txt
