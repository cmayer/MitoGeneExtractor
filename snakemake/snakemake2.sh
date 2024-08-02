#!/bin/bash
#SBATCH --job-name=MGE-test
#SBATCH --partition=hour
#SBATCH --output=MGE-test.out
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=20
#SBATCH --error=MGE-test.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate snakemake

snakemake \
   --snakefile Snakefile2 \
   --cores 20
