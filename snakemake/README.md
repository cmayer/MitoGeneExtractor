# MGE via snakemake with added functionality for BGE #
## Requirements: ##
- Snakefile
- config.yaml (containing:)
  - Name of MitoGeneExtractor run (will be used for naming summary stats output and concatenated consensus fasta file).
  - Path to 'MitoGeneExtractor-vx.y.z' file (see installation guidance on main readme)
  - Path to 'samples.csv' (example below)
  - Path to 'protein_references.csv' (example below)
  - Path to output directory (new directories will be created)
  - Gene of interest (e.g. cox1)
  - Parameters: r (Exonerate relative score threshold) and s (Exonerate minimum score threshold)
  - Read pre-processing mode. Options: 'merge' or 'conat'.
    - merge = adapter- and poly g-trimming, deduplication and PE read merging (fastp)-> 'cleaning' of sequence headers -> MGE
    - concat = gunzip and 'cleaning' of sequence headers -> adapter- and poly g-trimming, and deduplication (fastp) -> concatenation of PE reads -> read trimming (cutadapt) -> MGE
- Activated conda env - See mge_env.yaml
- Can be run on a cluster using 'snakemake.sh'





## Running: ##
### 1. Clone github repository ###
```bash
conda install anaconda::git # if not already installed
git clone https://github.com/SchistoDan/MitoGeneExtractor.git path/to/installation/dir
cd path/to/installation/dir
git status
```

### 2. Generate samples.csv ###
- Can be created via [sample-processing](https://github.com/bge-barcoding/sample-processing) workflow, or manually.
- Column 1 can be named 'ID', 'process_id', 'Process ID', 'process id', 'Process id', 'PROCESS ID', 'sample', 'SAMPLE', or 'Sample'.
- Reads must be PE, and can be gzipped or gunzipped.
  
**samples.csv example**
| ID | forward | reverse | taxid |
| --- | --- | --- | --- |
| BSNHM002-24  | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177658 |
| BSNHM038-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177627 |
| BSNHM046-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 3084599 |

### 3. Fetch sample-specific protein references using 1_gene_fetch.py ###
- [1_gene_fetch.py](https://github.com/SchistoDan/MitoGeneExtractor/blob/main/snakemake/1_gene_fetch.py) fetches sample-specific protein (pseudo)references using taxonomic ids and creates protein_references.csv (example below) required in config.yaml 
- 1_gene_fetch.py usage:
```
python 1_gene_fetch.py <gene_name> <output_directory> <samples.csv> 
  <gene_name>: Name of gene to search for in NCBI RefSeq database (e.g., cox1/COX1).
  <output_directory>: Path to directory to save output files (will save .fasta files and summary CSV in this directory). The directory will be created if it does not exist.
  <samples.csv>: Path to input CSV file containing Process IDs (ID column) and TaxIDs (taxid column).
```
- 'Checkpointing 'available: If the script fails during a run, it can be rerun using the same inputs and it will skip IDs with entries already in the protein_references.csv and with .fasta files already present in the output directory.
- Manually review the protein_references.csv after running as homonyms may lead to incorrect protein references being fetched on occasion.

**protein_references.csv example** 
| ID | taxid | accession_number | sequence_length | matched_rank | ncbi_taxonomy | reference_name | reference_path |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BSNHM002-24 | 177658 | AHF21732.1 | 510 | genus | Eukaryota, ..., Apataniinae, Apatania | BSNHM002-24 | abs/path/to/protein_references/BSNHM002-24.fasta | 


### 4. Edit config.yaml for run ###
- Update config.yaml with neccessary paths and variables.

### 5. Run snakefile ###
- Locally: snakemake --snakefile <Snakefile> --configfile <config.yaml>
  - Optional: '-n' for dry run. '-p' to print shell commands to log. '--unlock' to unlock directory after 'failed' run.
- Cluster: See snakemake.sh
  - Optional: '--rerun-incomplete' to resume a previously failed run.

### 6. Results structure ###
'''
results/
├── alignment/
│   └── alignment_files.log
├── consensus/
│   └── <run_name>.fasta
├── err/
├── logs/
├── out/
├── trimmed_data/
│   └── (reports/)
├── (raw_data)/
├── cleanup_complete.txt
├── <run_name>.csv
└── <run_name>-contaminants.csv
'''

N.B. 'raw_data/' only produced when running pipeline in 'concat' (i.e. 'standard') mode. 'reports/' only generated when running pipeline in 'merge' (i.e. fastp-merge) mode.


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

## To do ##
- Integrate BBsplit contam screen as an additional pre-processing mode for screening likely/known contaminant sequences.
- Integrate 1_gene_fetch.py into snakefile.
- Make Workflow Hub compatible.
- Generate RO-crates.
  
### Test run ###
- Raw reads for 12 test samples can be downloaded [here](https://naturalhistorymuseum-my.sharepoint.com/personal/b_price_nhm_ac_uk/_layouts/15/onedrive.aspx?ct=1723035606962&or=Teams%2DHL&ga=1&LOF=1&id=%2Fpersonal%2Fb%5Fprice%5Fnhm%5Fac%5Fuk%2FDocuments%2F%5Ftemp%2F%5FBGEexamples4Felix%2F1%5Fraw%5Fdata). Each read pair must be in seperate subdirectories under a parent directory that can be called anything
- samples sheet (BGE_test_samples.csv) provided (paths to reads and references need to be altered to where you stored the reads)
- protein_references sheet (gene_fetch_BGE_test_data.csv) provided (in protein_references/).
