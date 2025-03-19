# MGE snakemake pipeline #
Snakemake workflow for extracting high-quality barcoding gene regions, built around MitoGeneExtractor and adapted for genome skims of museum speicmens. 

# Requirements: #
- Paired-end reads in .fastq.gz or .fastq format
- samples_file.csv
- sequence_references_file.csv
- Activated conda env (see mge_env.yaml)



# Running: #
## 1. Set up conda environment and clone github repository ##
- Install conda.
```bash
conda env create -f /path/to/mge_env.yaml
git clone https://github.com/bge-barcoding/MitoGeneExtractor-BGE.git [path/to/installation/dir]
cd [path/to/installation/dir]
git status
```

## 2. Generate samples.csv ###
- Must contain ID, forward (read paths), reverse (read paths), and taxid columns (see below for example). Column 1 can be named 'ID', 'process_id', 'Process ID', 'process id', 'Process id', 'PROCESS ID', 'sample', 'SAMPLE', or 'Sample'.
- Can be created manually, or via [sample-processing](https://github.com/bge-barcoding/sample-processing) workflow.
- Due to regex matching in order to aggregate statistics, the sample ID will be considered as the string before the first underscore. It is therefore recommended that sample names do not use '_' characters. E.g. BSNHM002-24 instead of BSNHM002_24, or P3-1-A10-2-G1 instead of P3_1_A10_2_G1.
- taxid's can be found manually by searching the expected species/genus/family of each sample in the [NCBI taxonomy database](https://www.ncbi.nlm.nih.gov/taxonomy).
  
**samples.csv example**
| ID | forward | reverse | taxid |
| --- | --- | --- | --- |
| BSNHM002-24  | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177658 |
| BSNHM038-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177627 |
| BSNHM046-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 3084599 |

## 3. Gathering sample-specific pseudo-references using Gene Fetch ##
- gene_fetch.py retrieves the protein pseudo-references for each sample using the samples taxonomic identifier (taxid).
- The tool can fetch both protein and nucleotide sequences from NCBI databases for a given gene. See [gene_fetch](https://github.com/SchistoDan/gene_fetch/tree/main) repository for more information.
- 'gene_fetch.py' also provided in scripts/.

## 4. Customising snakemake configuration file ##
- Pre-processing modes:
  - 'merge' = adapter- and poly g-trimming, deduplication and PE read merging (fastp) -> 'cleaning' of sequence headers -> MGE
  - 'concat' = gunzip and 'cleaning' of sequence headers -> adapter- and poly g-trimming, and deduplication (fastp) -> concatenation of PE reads -> read trimming (Trim Galore (cutadapt)) -> MGE
![image](https://github.com/user-attachments/assets/21ce71b2-42df-4442-bcde-d41ee89fa3c1)
- Update config.yaml with neccessary paths and variables.
  - See [MitoGeneExtractor README.md](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/README.md) for explanation of Exonernate run paramters.
  - See [fasta_cleaner.py repository](https://github.com/bge-barcoding/fasta-cleaner) for information on filtering variables and thresholds (default below suitable in most cases).
  - See [Gene Fetch repository](https://github.com/bge-barcoding/gene_fetch) for guidance on creating [sequence_reference_file.csv](https://github.com/bge-barcoding/gene_fetch?tab=readme-ov-file#normal-mode).
```
# config.yaml
## Run parameters
# MGE run name identifier
run_name: "mge_concat_r1_13_15_s50_100"

# Path to MGE installation (MitoGeneExtractor-vX.X.X file)
mge_path: "/absolute/path/to/MitoGeneExtractor/MitoGeneExtractor-v1.9.5"

# Path to samples_file.csv
samples_file: "/absolute/path/to/samples_file.csv"

# Path to sequence_reference_file.csv
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

## 5. Run snakemake workflow (via cluster) ##
- [snakemake.sh](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/snakemake/snakemake.sh) cluster submission script:
  - --cluster: Defines how jobs are submitted to SLURM.
    - --parsable: Tells sbatch to only return the job ID.
    - --signal=USR2@90: Sends a signal 90 seconds before job time limit (for clean termination).
    - --cluster-config: Lists path to cluster configuration file (see below for explanation) and enables use of rule-specific resource requirements.
    - --mem={resources.mem_mb}MB: Dynamically sets memory allocation by using the mem parameter from each snakemake rule.
    - --cpus-per-task={threads}: Uses the threads parameter from each snakemake rule.
    - --output=slurm-%j-%x.out & --error=slurm-%j-%x.err: Sets naming convention for .out and .err files of individual jobs. '%j' = Job ID. '%x' = Job name.
  - --snakefile: List path to Snakefile.
  - --configfile: List path to snakemake workflow configuration file.
  - --latency-wait: Required when working on a distributed filesystem (e.g. NFS/GPFS). Set at 60 seconds by default. May be necessary to increase if experiencing latency/'missing' file issues.
  - --rerun-incomplete: If the snakemake workflow fails or is stopped for any reason, adding this option to the run command will enable the workflow to carry on from where it stopped.

- [cluster_config.yaml]()
  - Enables job-specific resource allocation based on job requirements and system capability. 
  - Default: Sets the default/minimum parameters to fallback on if not listed for a specific rule
- The aforementioned files will need tweaking to run on your cluster set up.

### Resource allocation ###
- SBATCH scheduler job parameters:
  - 'cpus-per-task' and 'mem' only apply to the 'master' job that coordinates the workflow and submits individual jobs to the job scheduler. Specified resources are only allocated for this 'master' job. Therefore, only 5-15G of memory and 2-4 CPUs are likely needed. It is recommended to set a relatively 'long' partition (e.g. several days-week) for this 'master' job, as it will be active for the entire run.
- Rule-specific resources in Snakefile:
  - Each rule can specify threads and memory resources (in Mb). These are the base values Snakemake uses initially.
- Cluster config values:
  - The 'cpus-per-task' and 'mem' values override or supplement rule-specific values in the Snakefile. If a rule doesn't specify resources, it will fallback to the listed defaults.
- Global limits:
  - '--cores': Limits total cores used across all concurrent jobs in the workflow.
  - '--jobs': Maximum number of simultaneous cluster jobs that will be run. E.g., '--jobs 25' = Up to 25 separate SLURM jobs can run simultaneously. 100 parallel is the maximum allowe
 
  
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


## To do ##
- Split Snakefile into rules (goes towards Workflow Hub compatibility)
- Update test data and associated files.
- Integrate 1_gene_fetch.py into snakefile.
- Make Workflow Hub compatible.
- Generate RO-crates.
  
## Test run ##
- Raw reads for 12 test samples can be downloaded [here](https://naturalhistorymuseum-my.sharepoint.com/personal/b_price_nhm_ac_uk/_layouts/15/onedrive.aspx?ct=1723035606962&or=Teams%2DHL&ga=1&LOF=1&id=%2Fpersonal%2Fb%5Fprice%5Fnhm%5Fac%5Fuk%2FDocuments%2F%5Ftemp%2F%5FBGEexamples4Felix%2F1%5Fraw%5Fdata). Each read pair must be in seperate subdirectories under a parent directory that can be called anything
- samples sheet (BGE_test_samples.csv) provided (paths to reads and references need to be altered to where you stored the reads)
- protein_references sheet (gene_fetch_BGE_test_data.csv) provided (in protein_references/).


## Contributing ##
- Please feel free to submit issues, fork the repository, and create pull requests for any improvements.
