#!/bin/bash
#SBATCH --job-name=MGE
#SBATCH --partition=day
#SBATCH --output=%x.out
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=16
#SBATCH --error=%x.err

source ~/miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate mge_env


snakemake \
   --snakefile Snakefile \
   --configfile config.yaml \
   --cores 16
