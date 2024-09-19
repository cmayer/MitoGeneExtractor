# MGE via snakemake with added functionality for BGE #
## Requirements: ##
- snakefile
- config.yaml (containing:)
  - path to'samples.csv' (example below - created via [sample-processing](https://github.com/bge-barcoding/sample-processing) workflow)
  - path to 'protein_references.csv' (example below - created using [1_gene_fetch.py](https://github.com/SchistoDan/MitoGeneExtractor/blob/main/snakemake/1_gene_fetch.py))
  - path to output directory (new directories will be created)
  - gene of interest (e.g. cox1)
  - r (Exonerate relative score threshold) and s (Exonerate minimum score threshold) parameter changes
- Activated conda env (with Snakemake, TrimGalore, Exonerate, Fastp, Biopython and Numpy installed). See mge_env.yaml
- Can be run on cluster using 'snakemake.sh'

**samples.csv example**
| ID | forward | reverse | taxid |
| --- | --- | --- | --- |
| BSNHM002-24  | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177658 |
| BSNHM038-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177627 |
| BSNHM046-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 3084599 |

**protein_references.csv example** 
| ID | matched_term | accession_number | reference_path | reference_name |
| --- | --- | --- | --- | --- |
| BSNHM002-24  | Apataniidae | YP_010586031.1 | abs/path/to/protein_references/BSNHM002-24.fasta | BSNHM002-24 |
| BSNHM038-24 | Trichoptera | YP_010894795.1 | abs/path/to/protein_references/BSNHM038-24.fasta | BSNHM038-24 |
| BSNHM046-24 | Polycentropodidae | YP_010426350.1 | abs/path/to/protein_references/BSNHM046-24.fasta | BSNHM046-24 |
  
  

## Running: ##
**snakefile:** 
- Uses config.yaml
  - Contains path to 'samples.csv', an output directory, and 'protein_references.csv'
  - Contains different parameter configurations for -s and -r
- Uses taxa-specific references from protein_references.csv (requires 1-gene_fetch.py to be run (test data protein_references.csv example provided)).
- 3_mge_tidy-snakemake.py and 4_mge_stats-snakemake.py functionality integrated into snakefile

**nakefile-fastpmerge:** 
- Utilises a fastp merge and post-merge concatenation of merged and unpaired reads as an alternative to the 'standard' MGE pipeline.

**snakefile_contam_refs**
- Uses a multi-fasta files containing target and non-target/contaminant protein reference sequences to map contaminant reads to protein references of common museum specimen contaminants (see 'add_contam_refs.py' in scripts  section below)

**snakefile_fastp_contam_refs**
- Same as above but employing fastp merge variant of MGE pipeline.

**Test run**
- Raw reads for 12 test samples can be downloaded [here](https://naturalhistorymuseum-my.sharepoint.com/personal/b_price_nhm_ac_uk/_layouts/15/onedrive.aspx?ct=1723035606962&or=Teams%2DHL&ga=1&LOF=1&id=%2Fpersonal%2Fb%5Fprice%5Fnhm%5Fac%5Fuk%2FDocuments%2F%5Ftemp%2F%5FBGEexamples4Felix%2F1%5Fraw%5Fdata). Each read pair must be in seperate subdirectories under a parent directory that can be called anything
- samples sheet (BGE_test_samples.csv) provided (paths to reads and references need to be altered to where you stored the reads)
- protein_references sheet (gene_fetch_BGE_test_data.csv) provided (in protein_references/).


## Scripts ##
- 1_gene_fetch.py = creates protein_references.csv
- add_contam_refs.py = If running Snakefile_contam_refs or Snakefile_fastp_contam_refs, protein reference files fetched ysing 1_gene_fetch.py need 'contaminant' reference sequences added to fasta files (creating multi-fasta files containing target and non-target/contaminant protein reference sequences (see common_contaminants.fasta)
- mge_contam_stats.py = If running Snakefile_contam_refs or Snakefile_fastp_contam_refs, statistics output by MGE are incorrect, and so cannot be parsed from .out files. This script (incorporated into 'rule extract_stats_to_csv' in Snakefile_contam_refs and Snakefile_fastp_contam_refs will use alignment fasta files to generate correct statistics).

## To do ##
- Integrate 1_gene_fetch.py into snakefile.
- Make Workflow Hub compatible.
- Generate RO-crates.
  
