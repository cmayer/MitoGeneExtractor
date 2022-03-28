# Example analysis using a Snakemake workflow and MitoGeneExtractor to search for mitochondrial genes in multiple fastq files:

IMPORTANT: This example analysis is still under construction. The only example analysis that is fully implemented and documented can be found in the folder "example-analysis-for-MitoGeneExtractor".

Snakemake is a workflow management system which allows upscaling of data analyses in a reproducible way.
The usage of Snakemake is not required when using MitoGeneExtractor, but it is comparatively convenient to use if the user wants to extract mitochondrial genes from a large number of datasets. 

### Installation:
Snakemake relies on Python3, which must be installed. You can install Snakemake in various ways, but the usage of a package manager such as Anaconda is recommended. 
See https://snakemake.readthedocs.io/en/stable/getting_started/installation.html for further information. You can set up a Snakemake environment and install necessary software within the environment.

### Prerequisites:
Before starting the analyses, the user needs to provide all necessary input data. For the example workflow, you would need to provide the following files:

- A protein reference for MitoGeneExtractor in FASTA format  
- Snakefile (See example file in this folder)
- (Optional) If raw sequencing reads shall be preprocessed with the aim to remove sequencing adaptors and low quality regions, install a program such as TrimGalore. TrimGalore, which is a wrapper to the cutadapt software is used in the present workflow for this purpose.
<!-- Expert users can use the cutadapt.yaml file to install  file (Anaconda environment for cutadapt safed in .yaml file)  -->
The user can also incorporate different quality and adaptor trimming software. Please also have in mind that you might want to adjust the example workflow according to your directory structure. Further, if other versions of the used software are installed, some parameter names might have changed.
-The configuration file in .yaml format. See config.yaml file in this folder.  

- Input raw data, e.g. sequencing reads in SRA file format as in this example. For this example download the following SRA files from NCBI:
https://www.ncbi.nlm.nih.gov/sra/?term=SRR12554982
https://www.ncbi.nlm.nih.gov/sra/?term=SRR12554985

The recommended way to download this data is to use the "prefetch" command from the [SRA toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software).

```{r, eval=TRUE}
prefetch SRR12554982
prefetch SRR12554985
```

### Workflow steps:
In our example workflow, the wildcard {sample} is replaced with the individual sample names found in the config.yaml file, followed by universal file extensions. The workflow from raw sequencing data to the desired consensus sequence will incorporate the following steps. Snakemake will determine dependencies between input and output file based on the user defined rules.

Convert the SRA files to the fastq format while retaining unique reads IDs:
```{r, eval=TRUE}
fasterq-dump --split-e --readids SRR12554982.sra
fasterq-dump --split-e --readids SRR12554985.sra
```

Concatenate paired-end libraries:
```{r, eval=TRUE}
cat SRR12554982*.fastq > SRR12554982_concat.fastq
cat SRR12554985*.fastq > SRR12554985_concat.fastq
```

Remove adapters and low quality regions:
```{r, eval=TRUE}
perl TrimGalore-0.6.6/trim_galore --no_report_file --dont_gzip --output_dir ./ SRR12554982_concat.fastq #results in outputfile SRR12554982_concat_trimmed.fq
perl TrimGalore-0.6.6/trim_galore --no_report_file --dont_gzip --output_dir ./ SRR12554985_concat.fastq #results in outputfile SRR12554985_concat_trimmed.fq
```

Transform FASTQ to FASTA format:
```{r, eval=TRUE}
awk '(NR-1)%4 == 0 || (NR-2)%4==0' SRR12554982_concat_trimmed.fq | tr '@' '>' > SRR12554982_trimmed.fas
awk '(NR-1)%4 == 0 || (NR-2)%4==0' SRR12554985_concat_trimmed.fq | tr '@' '>' > SRR12554985_trimmed.fas
```

Replace spaces with underscores in the sequence headers of a given fasta file:
```{r, eval=TRUE}
fasta2fasta-v1.4 -r ' _' -u SRR12554982_trimmed.fas SRR12554982_trimmed_newname.fas
fasta2fasta-v1.4 -r ' _' -u SRR12554985_trimmed.fas SRR12554985_trimmed_newname.fas
```

Finally, call MitoGeneExtractor:
```{r, eval=TRUE}
MitoGeneExtractor -d SRR12554982_trimmed.fas -p protein_reference.fas -o SRR12554982_out_alignment.fas -c SRR12554982_out_consensus.fas -V SRR12554982_vulgar.txt -n 0 -t 0.5 -r 1 -e /home/usr/bin/exonerate
MitoGeneExtractor -d SRR12554985_trimmed.fas -p protein_reference.fas -o SRR12554985_out_alignment.fas -c SRR12554985_out_consensus.fas -V SRR12554985_vulgar.txt -n 0 -t 0.5 -r 1 -e /home/usr/bin/exonerate
```

### Execution:
There are different ways to execute Snakemake. Which one works best for you depends on the structure of your computing system (for example local computer vs. cluster).
Once you activated the Snakemake environment, you should perform a dry run in order to test the executability of your workflow by typing ```snakemake -n```
Snakemake will tell you whether you Snakefile is syntactically correct and whether all required input data is available. The provided example Snakefile has the 'rule all' defined, which tells Snakemake to produce the results for all samples present in the config file.

For upscaling of the data analysis, you can simultaneously run several jobs (execution of rules) via typing ```snakemake --jobs 10```. Please have in mind that each job will need one core and Snakemake will need in addition one core to monitor and submit the jobs (so ```--jobs 10``` will actually need 11 cores).

An example command for Snakemake execution on a cluster would look like this:

```
module load anaconda3/2020.02 #load anaconda package manager
conda activate snakemake  #activate the snakemake environment
module load sratoolkit/2.10.8 #load the sratoolkit in this environment. 
snakemake --cluster \"qsub\" --jobs 3 --use-conda #"qsub" stands for the submission command for a SGE queuing system
```


For further information please refer to:  
https://snakemake.github.io  
https://snakemake.readthedocs.io/en/stable/
