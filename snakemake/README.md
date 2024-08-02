# MGE via snakemake with added functionality for BGE

## Three running options:
**Snakefile:** **<-- WORKS**
- Run using snakemake.sh
- Uses config.yaml (contains a path to BGE_test_samples.csv (lists sample names, and fwd and rev read paths))

**Snakefile2:** 
- Run using snakemake2.sh
- Uses config2.yaml (contains a path to BGE_test_samples_all.csv (lists sample names, fwd and rev read paths) and BOLD_output_mge_fetch_BGE_test_data.csv (contains sample names, matched term (from mge_fetch.py), reference name (accession #), and path to .fasta reference sequence)

**Snakefile3:** 
- Run using snakemake3.sh
- Uses config3.yaml (contains a path to BGE_test_samples_all.csv (lists sample names, fwd and rev read paths, matched term (from mge_fetch.py), reference name (accession #), and path to .fasta reference sequence))

## To do
- Get Snakefile2 working (i.e. so MGE can take a samples_file and protein_references_file as input). OR get Snakefile3 working (i.e. so MGE can take a combined samples-protein_references_file as input).
- Get MGE to output to scratch space (output_dir = "/gpfs/nhmfsa/bulk/share/data/mbl/share/scratch/MGE")
