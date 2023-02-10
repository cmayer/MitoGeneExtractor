# Example analysis using a Snakemake workflow and MitoGeneExtractor:
Since MitoGeneExtractor version 1.9.5, several steps of this example analysis can be directly executed within MitoGeneExtractor. This comprises merging of multiple input files and the reconstruction of more than one gene sequence.
In addition, implementing quality trimming or the analysis with various parameter combinations in the snakemake analysis workflow can help to scale the analysis to large data sets while maintaining an organized file system. 

### Prerequisites:
Before starting the analyses, the user needs to provide all necessary input data. Please also have in mind that you might want to adjust the example workflow according to your directory structure. Further, if other versions of the used software are installed, some parameter names might have changed.
For the example workflow, you would need to provide the following files:

- Protein reference(s) for MitoGeneExtractor in FASTA format (See references in this folder). References are assumed to be in the directory ./protein_references/
- Snakefile (See example file in this folder)
- A quality trimming software such as TrimGalore! installed. TrimGalore!, which is a wrapper to the cutadapt software is used in the present workflow for this purpose, but any other software can be used. 
- The configuration file in .yaml format (See config.yaml file in this folder) 
- Raw sequencing data. Raw data is assumed to be in the directory ./raw_data/

For this example download the following SRA files from NCBI:
https://www.ncbi.nlm.nih.gov/sra/?term=SRR12554982 and https://www.ncbi.nlm.nih.gov/sra/?term=SRR12554985

The recommended way to download this data is to use the "prefetch" command from the [SRA toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software).

```{r, eval=TRUE}
prefetch -O raw_data/ SRR12554982 
prefetch -O raw_data/ SRR12554985
```
The .sra files are downloaded to the "raw_data" directory

### Workflow steps:
In our example workflow, the wildcard {sample} is replaced with the individual sample names found in the config.yaml file, followed by universal file extensions. The workflow from raw sequencing data to the desired consensus sequence will incorporate the following steps. Snakemake will determine dependencies between input and output file based on the user defined rules.

Convert the SRA files to the fastq format with fastq-dump from the SRA toolit. Retain unique reads IDs:
```{r, eval=TRUE}
fastq-dump --split-e --readids -O raw_data raw_data/SRR12554982/SRR12554982.sra
fastq-dump --split-e --readids -O raw_data raw_data/SRR12554985/SRR12554985.sra
```

Concatenate paired-end libraries:
```{r, eval=TRUE}
cat raw_data/SRR12554982_1.fastq raw_data/SRR12554982_2.fastq > raw_data/SRR12554982_concat.fastq
cat raw_data/SRR12554985_1.fastq raw_data/SRR12554985_2.fastq > raw_data/SRR12554985_concat.fastq
```

Remove adapters and low quality regions:
```{r, eval=TRUE}
perl TrimGalore-0.6.6/trim_galore --no_report_file --dont_gzip --output_dir trimmed_data/ raw_data/SRR12554982_concat.fastq #results will be written to "trimmed_data" dir 
perl TrimGalore-0.6.6/trim_galore --no_report_file --dont_gzip --output_dir trimmed_data/ raw_data/SRR12554985_concat.fastq #results will be written to "trimmed_data" dir 
```

If you want to organize your MitoGeneExtractor results in different directories without using snakemake, make sure that the path exists before running MitoGeneExtractor. For example, type:
```
mkdir COX1
```
Snakemake will handle this for you directly, so this does not need to be included in the Snakefile.

Finally, call MitoGeneExtractor in default mode to reconstruct e.g. the COX1 gene
```{r, eval=TRUE}
~/MGE_test/MitoGeneExtractor/MitoGeneExtractor-v1.9.5 --report_gaps_mode 1 -q trimmed_data/SRR12554982_concat_trimmed.fq -p protein_references/Passeriformes_COX1.fasta -o COX1/SRR12554982_out_alignment.fas -c COX1/SRR12554982_out_consensus.fas -V COX1/SRR12554982_vulgar.txt -e ~/bin/exonerate

~/MGE_test/MitoGeneExtractor/MitoGeneExtractor-v1.9.5 --report_gaps_mode 1 -q trimmed_data/SRR12554985_concat_trimmed.fq -p protein_references/Passeriformes_COX1.fasta -o COX1/SRR12554985_out_alignment.fas -c COX1/SRR12554985_out_consensus.fas -V COX1/SRR12554985_vulgar.txt -e ~/bin/exonerate
```
In this example, the protein reference file is stored in the protein_references/ directory and contains only the COX1 amino acid sequence. In the provided example workflow, the ND5 gene is additionally reconstructed, because it illustrates, how {wildcards} are used. Note, that it is **NOT** necessary to provide the references individually when you use MitoGeneExtractor version 1.9.5 or higher. 

### Execution:
There are different ways to execute Snakemake. Which one works best for you depends on the structure of your computing system (for example local computer vs. cluster).
Once you activated the Snakemake environment, you should perform a dry run in order to test the executability of your workflow by typing ```snakemake -n```
Snakemake will tell you whether you Snakefile is syntactically correct and whether all required input data is available. The provided example Snakefile has the 'rule all' defined, which tells Snakemake to produce the results for all samples present in the config file.

For upscaling of the data analysis, you can simultaneously run several jobs (execution of rules) via typing ```snakemake --jobs 10``` or  ```snakemake --cores 10```. Please have in mind that each job will need one core and Snakemake will need in addition one core to monitor and submit the jobs (so ```--jobs 10``` will actually need 11 cores).

An example command for Snakemake execution on a cluster based on a SGE queuing system would look like this:

```
module load anaconda3/2020.02 #load anaconda package manager
conda activate snakemake  #activate the snakemake environment
module load sratoolkit/2.10.8 #load the sratoolkit in this environment. 
```
You could execute snakemake via the command line. In this mode, snakemake will run on the head node (which is typically the node with submission rights) and submit one after another individual job scripts for each rule execution:
```
snakemake --cluster \"qsub\" --jobs 3 #"qsub" stands for the submission command for a SGE system. Snakemake will submit 3 jobs simultaneosuly.
```

Alternatively, you can run snakemake via submitting a jobscript with the following command:

```
snakemake --cores 3 #Snakemakem will distribute the rule execution over 3 cores. In our example, each rule claim one core, allowing to run 3 jobs simultaneosuly.
```

For further information please refer to:  
https://snakemake.github.io  
https://snakemake.readthedocs.io/en/stable/
