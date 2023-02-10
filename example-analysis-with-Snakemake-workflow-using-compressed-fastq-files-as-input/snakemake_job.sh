#!bin/bash
#
#$ -S /bin/bash
#$ -N snakejob
#$ -cwd

module load anaconda3/2020.02
conda activate snakemake
module load sratoolkit/2.10.8

##provide one thread less to the workflow, because snakemake uses one for job monitoring
a=1; b=$NSLOTS
THREADS=$((b-a))

snakemake --cores ${THREADS} 
