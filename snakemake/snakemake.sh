#!/bin/bash
##These SBATCH parameters (parition, mem, and cpus-per-task) are solely for the main snakemake process.
##These may need some tweaking, although 5-20G of memory and 1-4 CPUs is likely more than enough.
#SBATCH --job-name=mge
#SBATCH --partition=long
#SBATCH --mem=15G
#SBATCH --cpus-per-task=4



## Conda environment

# Set path to conda explicitly
export PATH="/path/to/conda/bin:$PATH"

# Source conda.sh
source /path/to/conda/etc/profile.d/conda.sh

# Activate conda env
conda activate mge_env




##Snakemake

# Unlock directory
snakemake --snakefile ./workflow/Snakefile --configfile ./config/config.yaml --unlock

# Run snakemake workflow on cluster
snakemake \
          --cluster "sbatch --parsable --partition=day --signal=USR2@90 --mem={resources.mem_mb}MB --cpus-per-task={threads} --output=slurm-%j-%x.out --error=slurm-%j-%x.err" \
          --cluster-config ./config/cluster_config.yaml \
          --cores 28 \
          --jobs 15 \
          --snakefile ./workflow/Snakefile \
          --configfile ./config/config.yaml \
          --rerun-incomplete \
          --latency-wait 60
