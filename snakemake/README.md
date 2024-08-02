# MGE via snakemake with added functionality for BGE

## Two running options:
**Snakefile:**
- Run using snakemake.sh
- Uses config.yaml (contains a path to BGE_test_samples.csv (lists sample names, and fwd and rev read paths))

**Snakefile2:** 
- Run using snakemake2.sh
- Uses config2.yaml (contains a path to BGE_test_samples.csv (lists sample names, and fwd and rev read paths) & BOLD_output_mge_fetch_BGE_test_data.csv (lists sample names, reference name (accession #) and path to .fasta reference sequence))

## To do
- Get Snakefile2 working (i.e. so MGE can take a samples_file and protein_references_file as input). If easier, the two .csv files can be merged.
- Get MGE to output to scratch space (output_dir = "/gpfs/nhmfsa/bulk/share/data/mbl/share/scratch/MGE")
