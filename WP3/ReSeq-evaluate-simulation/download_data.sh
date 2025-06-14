#!/bin/bash

#SBATCH --account=lp_h_master_thesis_volders_2025 
#SBATCH --clusters=genius
#SBATCH --partition=batch   # Name of Partition
#SBATCH --job-name=download-data # Name of job
#SBATCH --output=stdout-download-data   # Standard output name
#SBATCH --error=stderror-download-data    # Standard Errorname
#SBATCH --ntasks=1            # Number of tasks
#SBATCH --cpus-per-task=12     # Number of CPU cores
#SBATCH --time=0-10:00:00       # Wall time (format: d-hh:mm:ss)
#SBATCH --mem=10gb            # Amount of memory (units: gb, mb, kb)
export PATH="${VSC_DATA}/miniconda3/bin:${PATH}"

STAGING=/staging/leuven/stg_00156/thesis_projects_2025/lynnhenrotte/

source activate samtools

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793898/NA02325.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793898/NA02325.bam.bai

samtools view -h --fetch-pairs NA02325.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut/NA02325_cut.bam
rm -f NA02325.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793910/NA21939.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793910/NA21939.bam.bai

samtools view -h --fetch-pairs NA21939.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut/NA21939_cut.bam
rm -f NA21939.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793906/NA20381.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793906/NA20381.bam.bai

samtools view -h --fetch-pairs NA20381.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut/NA20381_cut.bam
rm -f NA20381.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793905/NA20273_20221005Run.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793905/NA20273_20221005Run.bam.bai

samtools view -h --fetch-pairs NA20273_20221005Run.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut/NA20273_20221005Run_cut.bam
rm -f NA20273_20221005Run.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793914/ND01039.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793914/ND01039.bam.bai

samtools view -h --fetch-pairs ND01039.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut/ND01039_cut.bam
rm -f ND01039.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793900/NA17221.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793900/NA17221.bam.bai

samtools view -h --fetch-pairs NA17221.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut/NA17221_cut.bam
rm -f NA17221.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR148/ERR14800375/NA05117_20221005Run.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR148/ERR14800375/NA05117_20221005Run.bam.bai

samtools view -h --fetch-pairs NA05117_20221005Run.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut/NA05117_20221005Run_cut.bam
rm -f NA05117_20221005Run.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793901/NA18668.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793901/NA18668.bam.bai

samtools view -h --fetch-pairs NA18668.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut2/NA18668_cut.bam
rm -f NA18668.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR148/ERR14800374/NA04520.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR148/ERR14800374/NA04520.bam.bai

samtools view -h --fetch-pairs NA04520.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut2/NA04520_cut.bam
rm -f NA04520.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793908/NA21075.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793908/NA21075.bam.bai

samtools view -h --fetch-pairs NA21075.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut2/NA21075_cut.bam
rm -f NA21075.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793909/NA21698.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793909/NA21698.bam.bai

samtools view -h --fetch-pairs NA21698.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut2/NA21698_cut.bam
rm -f NA21698.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR148/ERR14800371/NA03330.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR148/ERR14800371/NA03330.bam.bai

samtools view -h --fetch-pairs NA03330.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut2/NA03330_cut.bam
rm -f NA03330.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793903/NA19239.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793903/NA19239.bam.bai

samtools view -h --fetch-pairs NA19239.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut2/NA19239_cut.bam
rm -f NA19239.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793913/NA23710_20221005Run.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793913/NA23710_20221005Run.bam.bai

samtools view -h --fetch-pairs NA23710_20221005Run.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut2/NA23710_20221005Run_cut.bam
rm -f NA23710_20221005Run.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793897/NA00343.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793897/NA00343.bam.bai

samtools view -h --fetch-pairs NA00343.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut2/NA00343_cut.bam
rm -f NA00343.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793899/NA14626.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR147/ERR14793899/NA14626.bam.bai

samtools view -h --fetch-pairs NA14626.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut2/NA14626_cut.bam
rm -f NA14626.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR148/ERR14800372/NA04372.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR148/ERR14800372/NA04372.bam.bai

samtools view -h --fetch-pairs NA04372.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut2/NA04372_cut.bam
rm -f NA04372.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR148/ERR14800377/NA10283.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR148/ERR14800377/NA10283.bam.bai

samtools view -h --fetch-pairs NA10283.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut2/NA10283_cut.bam
rm -f NA10283.bam*

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR148/ERR14800380/NA13480.bam
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR148/ERR14800380/NA13480.bam.bai

samtools view -h --fetch-pairs NA13480.bam "chr10:94400000-95071000" > ${STAGING}real_data/real_data_BAMs_cut2/NA13480_cut.bam
rm -f NA13480.bam*
