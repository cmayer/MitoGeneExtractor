# MGE snakemake pipeline #
Snakemake workflow for extraction of barcoding regions built around MitoGeneExtractor and adapted for genome skims of museum speicmens. 

# Requirements: #
- Snakefile
- config.yaml (see below)
- Activated conda env - See mge_env.yaml
- Can be run on a cluster using 'snakemake.sh'





# Running: #
## 1. Set up conda environment and clone github repository ##
- Install conda.
```bash
conda env create -f mge_env.yaml
git clone https://github.com/bge-barcoding/MitoGeneExtractor-BGE.git [path/to/installation/dir]
cd [path/to/installation/dir]
git status
```

## 2. Generate samples.csv ###
- Can be created manually, or via [sample-processing](https://github.com/bge-barcoding/sample-processing) workflow.
- Column 1 can be named 'ID', 'process_id', 'Process ID', 'process id', 'Process id', 'PROCESS ID', 'sample', 'SAMPLE', or 'Sample'.
- Due to regex matching to aggregate statistics, the sample ID will be considered as the string before the first '_'. It is therefore recommended that sample names do not use '_' characters. E.g. BSNHM002-24 instead of BSNHM002_24, or P3-1-A10-2-G1 instead of P3_1_A10_2_G1.
- Reads must be PE, either be compressed or uncompressed.
  
**samples.csv example**
| ID | forward | reverse | taxid |
| --- | --- | --- | --- |
| BSNHM002-24  | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177658 |
| BSNHM038-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177627 |
| BSNHM046-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 3084599 |

## 3. Fetch sample-specific protein references using 1_gene_fetch.py ##
### Gene Fetch
A Python tool for retrieving protein and/or gene sequences from NCBI databases. The script can fetch both protein and nucleotide sequences for a given gene across multiple taxa, with support for traversing taxonomic hierarchies when sequences aren't available at the given taxonomic level (dictated by input taxid). See [gene_fetch](https://github.com/SchistoDan/gene_fetch/tree/main) repository for more information. 1_gene_fetch.py provided in /scripts.

## 4. Edit config.yaml for run ##
- Update config.yaml with neccessary paths and variables.
- 'merge' preprocessing = adapter- and poly g-trimming, deduplication and PE read merging (fastp)-> 'cleaning' of sequence headers -> MGE
- 'concat' preprocessing = gunzip and 'cleaning' of sequence headers -> adapter- and poly g-trimming, and deduplication (fastp) -> concatenation of PE reads -> read trimming (cutadapt) -> MGE
```
# MGE run name identifier
run_name: "mge_concat_r1_13_15_s50_100"
# Path to MGE installation
mge_path: "path/to/MitoGeneExtractor-v1.9.5"
# Path to samples.csv
samples_file: "path/to/samples.csv"
# Path to protein_references.csv
sequence_reference_file: "path/to/sequence_references.csv"
# Path to output directory. Will make final dir if it does not exist already.
output_dir: "path/to/results/mge_concat_r1_13_15_s50_100"
# Gene(s) of interest
genes:
  - cox1
# Exonerate relative score threshold parameter
r:
  - 1
  - 1.3
  - 1.5
# Exonerate minimum score threshold parameter
s:
  - 50
  - 100
#Read pre-processing method. Options: "merge" or "concat" (fastp-merge (i.e. 'merge') or standard PE fastq concatenation (i.e. 'concat')). 
preprocessing_mode: "concat"
```

## 5. Run snakefile ##
- Locally: snakemake --snakefile <Snakefile> --configfile <config.yaml>
  - Optional: '-n' for dry run. '-p' to print shell commands to log. '--unlock' to unlock directory after 'failed' run.
- Cluster: See snakemake.sh
  - Optional: '--rerun-incomplete' to resume a previously failed run.

## 6. Results structure ##
```
results/
├── alignment/
│   └── alignment_files.log
│   └── sample1_alignment.fasta
│   └── sample2_alignment.fasta
│   └── ...
├── consensus/
│   └── sample1_consensuys.fasta
│   └── sample2_consensus.fasta
│   └── ...
│   └── <run_name>.fasta
├── err/
│   └── sample1.err
│   └── sample2.err
│   └── ...
├── logs/
│   └── sample1.log
│   └── sample2.log
│   └── ...
├── out/
│   └── sample1.out
│   └── sample2.out
│   └── ...
├── trimmed_data/
│   └── ...
│   └── (reports/)
├── cleanup_complete.txt
├── <run_name>.csv
└── <run_name>-contaminants.csv
```

## Scripts ##
- **1_gene_fetch.py** = Fetches protein references for each sample using taxids from samples.csv to query NCBI DBs (via BioEntrez API). Fetches closest reference available to input taxid.
- **mge_stats.py** = This script (incorporated into 'rule extract_stats_to_csv') will use alignment fasta files and MGE.out files to generate summary statistics for each sample.
```
python ./scripts/mge_stats.py <log_file> <output_file> <out_file_dir>
  log_file: Path to a text file containing a list of FASTA file paths (one per line) # Produced by 'create_alignment_log' rule
  output_file: Name of the output CSV file (will also generate a -contaminants.csv variant)
  out_file_dir: Directory containing MGE .out files with barcode sequence statistics
```
- **fasta_compare.py** = Supplementary script that can be run after the MGE pipeline is finished. It will compare barcodes produced using different parameter combinations (from one run or multiple runs) for each sample, ranks each barcode 1-5 based on [BOLD BIN criteria](https://v3.boldsystems.org/index.php/resources/handbook?chapter=2_databases.html&section=bins), and select the 'best' (BOLD BIN compliant) barcode.
```
python ./scripts/fasta_compare_new.py OUTPUT_CSV OUTPUT_FASTA OUTPUT_BARCODE_FASTA INPUT_FILES --log-file LOG_FILE --verbose
  OUTPUT_CSV: Path where the analysis results CSV file will be saved
  OUTPUT_FASTA: Path where the best full sequences FASTA file will be saved
  OUTPUT_BARCODE_FASTA: Path where the best barcode sequences FASTA file will be saved
  INPUT_FILES: One or more input FASTA files to analyze (space-separated)
  
  Optional:
  --log-file LOG_FILE: Specify a custom path for the log file (default: creates timestamped log in current directory)
  --verbose, -v: Enable detailed debug logging
```
## Workflow ##
...tbc
## To do ##
- Integrate 1_gene_fetch.py into snakefile.
- Make Workflow Hub compatible.
- Generate RO-crates.
  
## Test run ##
- Raw reads for 12 test samples can be downloaded [here](https://naturalhistorymuseum-my.sharepoint.com/personal/b_price_nhm_ac_uk/_layouts/15/onedrive.aspx?ct=1723035606962&or=Teams%2DHL&ga=1&LOF=1&id=%2Fpersonal%2Fb%5Fprice%5Fnhm%5Fac%5Fuk%2FDocuments%2F%5Ftemp%2F%5FBGEexamples4Felix%2F1%5Fraw%5Fdata). Each read pair must be in seperate subdirectories under a parent directory that can be called anything
- samples sheet (BGE_test_samples.csv) provided (paths to reads and references need to be altered to where you stored the reads)
- protein_references sheet (gene_fetch_BGE_test_data.csv) provided (in protein_references/).


## Contributing ##

- Please gfeel free to submit issues, fork the repository, and create pull requests for any improvements.
