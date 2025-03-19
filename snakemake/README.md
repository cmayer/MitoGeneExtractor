# MGE snakemake pipeline #
Snakemake workflow for extracting high-quality barcoding gene regions, built around MitoGeneExtractor and adapted for genome skims of museum speicmens. 

# Requirements: #
- Snakefile
- config.yaml
- Activated conda env (see mge_env.yaml)
- Can be run on SLURM cluster using 'snakemake.sh'
- cluster_config.yaml (if running on SLURM cluster)






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

## 3. Gathering sample-specific pseudo-references using Gene Fetch ##
A Python tool for retrieving protein and/or gene sequences from NCBI databases. The script can fetch both protein and nucleotide sequences for a given gene across multiple taxa, with support for traversing taxonomic hierarchies when sequences aren't available at the given taxonomic level (dictated by input taxid). See [gene_fetch](https://github.com/SchistoDan/gene_fetch/tree/main) repository for more information. 1_gene_fetch.py provided in scripts/.

## 4. Customising snakemake configuration file ##
- Update config.yaml with neccessary paths and variables.
- 'merge' preprocessing = adapter- and poly g-trimming, deduplication and PE read merging (fastp) -> 'cleaning' of sequence headers -> MGE
- 'concat' preprocessing = gunzip and 'cleaning' of sequence headers -> adapter- and poly g-trimming, and deduplication (fastp) -> concatenation of PE reads -> read trimming (Trim Galore (cutadapt)) -> MGE
![image](https://github.com/user-attachments/assets/b9c291e0-6146-4b5d-ab24-b4abefff2c0f)

```
## Run parameters
# MGE run name identifier
run_name: "mge_concat_r1_13_15_s50_100"

# Path to MGE installation (MitoGeneExtractor-vX.X.X file)
mge_path: "/absolute/path/to/MitoGeneExtractor/MitoGeneExtractor-v1.9.5"

# Path to samples.csv
samples_file: "/absolute/path/to/samples_file.csv"

# Path to references.csv
sequence_reference_file: "/absolute/path/to/sequence_reference_file.csv"

# Path to output directory. Will make final dir if it does not exist already
output_dir: "/absolute/path/to/output/directory"

## MGE parameters
# Exonerate relative score threshold parameter
r:
  - 1
  - 1.3
  - 1.5
# Exonerate minimum score threshold parameter
s:
  - 50
  - 100
  
## Read pre-processing method
# Options: "merge" or "concat" (fastp-merge (i.e. 'merge') or standard PE fastq concatenation (i.e. 'concat')). 
preprocessing_mode: "concat"

## Consensus sequence post-processing (using fasta_cleaner.py)
fasta_cleaner:
  consensus_threshold: 0.5
  human_threshold: 0.95
  at_difference: 0.1
  at_mode: "absolute"
  outlier_percentile: 90.0
  disable_human: false
  disable_at: false
  disable_outliers: false
  reference_dir: null 
```

## 5. Run workflow (via snakemake.sh) ##
- 
These SBATCH parameters (parition, mem, and cpus-per-task) are solely for the main snakemake process.
These may need some tweaking, although 5-20G of memory is likely more than enough.

The Master Job (Initial SBATCH Script)
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
This only allocates resources for the master job that runs MitoComp itself
This job primarily just coordinates and submits individual rule jobs to SLURM
It doesn't need much computational power - it's mostly managing the workflow

Snakemake Parallelism Settings
-s "--cores 28 --jobs 25"
--jobs 25: Up to 25 separate SLURM jobs can run simultaneously
--cores 28: Primarily affects any rules that might run locally (on the master node)

Cluster Config File
default:
   ntasks: 4
   mem: 50G
get_organelle:
   mem: 20G
   job-name: GETORG
Each individual rule job gets its own SLURM allocation based on this config
Rules use the specific settings if defined, otherwise fall back to "default"
GetOrganelle will get 20GB memory and 4 tasks (from default)

## 6. Results structure ##
```
results/
├── alignment/                      # MGE sequence alignments
│   └── sample1_alignment.fasta
│   └── sample2_alignment.fasta
│   └── ...
├── consensus/                      # MGE consensus sequences   
│   └── sample1_consensus.fasta
│   └── sample2_consensus.fasta  
│   └── ...
│   └── <run_name>.fasta            # Concatenated consensus multi-fasta
├── err/                            # MGE error/alignment logs
│   └── sample1.err
│   └── sample2.err
│   └── ...
├── out/                            # MGE raw outputs/run logs
│   └── sample1.out
│   └── sample2.out
│   └── ...
├── logs/
│   └── mge/
│   └── (clean_headers/)            # Per-sample clean_headers_merge logs (if run in 'merge' mode)
│   └── (gunzip/)                   # Per-sample fastq_concat logs (if run in 'concat' mode)
│   └── (concat/)                   # Per-sample clean_headers_merge logs (if run in 'merge' mode)
│   └── (trim_galore/)              # Per-sample qulity_trim logs (if run in 'merge' mode)
│   └── (gunzip.log)                # Aggregated logs from gunzip_and_clean_headers rule (if run in 'concat' mode)
│   └── (concat_reads.log)          # Aggregated logs from fastq_concat rule (if run in 'concat' mode)
│   └── (trim_galore.log)           # Aggregated logs from quality_trim rule (if run in 'concat' mode)
│   └── (clean_headers.log)         # Aggregated logs from clean_headers_merge rule (if run in 'merge' mode)
│   └── alignment_files.log         # List of alignment files in 'alignment/'
│   └── concat_consensus.log        # Log from concatenate_fasta rule 
│   └── rename_complete.txt         # Confirmation of consensus sequence filename and header renaming complete 
│   └── rename_fasta.log            # Log from concatenate_fasta rule 
│   └── fasta_cleaner.log           # Log from fasta_cleaner rule/fasta_cleaner.py
│   └── mge_stats.log               # Log from extract_stats_to_csv rule
│   └── cleaning_complete.txt       # Confirmation cleanup_files (final) rule has run
├── trimmed_data/
│   └── reports/                    # Per-sample fastp HTML and JSON reports (and trim_galore reports if run in 'concat' mode)
│   └── sample1_fastp.out
│   └── sample1_fastp.err
│   └── (sample1_merged_clean.fq)   # If run in 'merge' mode
│   └── (sample1_R1_trimmed.fastq)  # If run in 'concat' mode
│   └── (sample1_R2_trimmed.fastq)  # If run in 'concat' mode
│   └── sample1_R1_trimmed.fq.gz
│   └── sample1_R2_trimmed.fq.gz
│   └── ...
├── fasta_cleaner/
│   └── consensus_seqs/             # 'Cleaned' consensus sequences per sample
│   └── filter_annotated_seqs/      # Aligned and annotated (kept/removed) sequences per sample
│   └── filter_pass_seqs/           # Aligned sequences passing 'cleaning' parameters per sample
│   └── logs/                       # Fasta_cleaner.py logs containing running parameters and raw stats per sample
│   └── metrics/                    # Scores for each kept and removed read in the input MGE alignment 
│   └── removed_at_seqs/            # Sequences removed due to AT% threshold per sample
│   └── removed_human_seqs/         # Sequences removed due to Human COI similarity, per sample
│   └── removed_outlier_seqs/       # Sequences removed due to base outliers, per sample
│   └── removed_ref_comparison_seqs/# Sequences removed due to reference sequence similarity (if supplied), per sample
│   └── all_consensus_seqeunces.fasta
│   └── combined_statistics.csv     # Summary stats file of parsed stats from logs
└── <run_name>.csv                  # Summary stats file produced by extract_stats_to_csv rule (see below for example)
```

## Integrated and supplementary scripts ##
See scripts/.
- [**1_gene_fetch.py**](https://github.com/bge-barcoding/gene_fetch) = Supplementary script that fetches protein references for each sample using taxids from samples.csv to query NCBI DBs (via BioEntrez API). Fetches closest reference available to input taxid. See [1_gene_fetch.py](https://github.com/bge-barcoding/gene_fetch) github repository for more information.
- [**fasta_cleaner.py**](https://github.com/SchistoDan/MitoGeneExtractor/blob/main/snakemake/scripts/fasta_cleaner.py) = This script (incorproated into 'fasta_cleaner' rule) 'cleans' MGE alignment files using AT% thresholds, base consensus similarity, human COI similarity, and (if supplied) reference sequence similarity. Outputs 'cleaned' consensus sequences for each sample. Modified from [fasta_cleaner.py](https://github.com/bge-barcoding/fasta-cleaner), see original github repository for more information.
- [**mge_stats.py**](https://github.com/SchistoDan/MitoGeneExtractor/blob/main/snakemake/scripts/mge_stats.py) = This script (incorporated into 'rule extract_stats_to_csv') uses alignment fasta files and MGE.out files to generate summary statistics for each sample.
- [**fasta_compare.py**](https://github.com/SchistoDan/MitoGeneExtractor/blob/main/snakemake/scripts/fasta_compare.py) = Supplementary script that can be run after the MGE pipeline is finished. It will compare barcodes produced using different parameter combinations (from one run or multiple runs) for each sample, ranks each barcode 1-5 based on [BOLD BIN criteria](https://v3.boldsystems.org/index.php/resources/handbook?chapter=2_databases.html&section=bins), and select the 'best' (BOLD BIN compliant) barcode.

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

- Please feel free to submit issues, fork the repository, and create pull requests for any improvements.
