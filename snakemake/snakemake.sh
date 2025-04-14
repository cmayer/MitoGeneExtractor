#!/bin/bash
##These SBATCH parameters (parition, mem, and cpus-per-task) are solely for the main snakemake process.
##These may need some tweaking, although 5-20G of memory and 1-4 CPUs is likely more than enough.
#SBATCH --job-name=mge
#SBATCH --partition=long
#SBATCH --mem=8G
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
snakemake --snakefile ./workflow/Snakefile \
  --configfile ./config/config.yaml \
  --cluster-config ./config/cluster_config.yaml \
  --cluster "sbatch --parsable --partition=[cluster.partition} --cpus-per-task={cluster.cpus-per-task} --mem={cluster.mem} --output=slurm-%j-%x.out --error=slurm-%j-%x.err" \
  --jobs 20 \
  --cores 32 \
  --rerun-incomplete \
  --latency-wait 60

