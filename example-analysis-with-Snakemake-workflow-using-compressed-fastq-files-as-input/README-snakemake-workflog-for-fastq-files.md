# Example analysis using a Snakemake workflow and MitoGeneExtractor to search for mitochondrial genes in multiple fastq files:

Snakemake is not a prerequisite to use MitoGeneExtractor, but it provides a convenient method to to analyse a large number of data sets.
Snakemake is a workflow management system which allows upscaling of data analyses in a reproducible way. 

### Installation:
Snakemake relies on Python3, which must be installed. You can install Snakemake in various ways, but the usage of a package manager such as Anaconda is recommended. 
See https://snakemake.readthedocs.io/en/stable/getting_started/installation.html for further information. You can set up a Snakemake environment and install necessary software within the environment.

### Prerequisites:
Before starting the analyses, the user needs to provide all necessary input data. For the example workflow, you would need to provide the following files:

- Snakefile (See example file in this folder)

- Input raw data, e.g. sequencing reads in SRA file format as in this example. For this example download the following SRA files from NCBI:
https://www.ncbi.nlm.nih.gov/sra/?term=SRR12554982
https://www.ncbi.nlm.nih.gov/sra/?term=SRR12554985

The recommended way to download this data is to use the "prefetch" command from the [SRA toolkidt](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software).

```{r, eval=TRUE}
prefetch SRR12554982
prefetch SRR12554985
```

Convert the SRA files to the fastq format:
```{r, eval=TRUE}
fasterq-dump --split-e SRR12554982.sra
fasterq-dump --split-e SRR12554985.sra
```

...

-The configuration file in .yaml format. See config.yaml file in this folder.  
-A protein reference for MitoGeneExtractor in FASTA format  

- (Optional) If raw sequencing reads shall be preprocessed with the aim to remove sequencing adaptors and low quality regions, install a program such as TrimGalore. TrimGalore, which is a wrapper to the cutadapt software is used in the present workflow for this purpose.
<!-- Expert users can use the cutadap.yaml file to install  file (Anaconda environment for cutadapt safed in .yaml file)  -->
The user can also incorporate different quality and adaptor trimming software. Please also have in mind that you might want to adjust the example workflow according to your directory structure. Further, if other versions of the used software are installed, some parameter names might have changed.

### Execution:
There are different ways to execute Snakemake. Which one works best for you depends on the structure of your computing system (for example local computer vs. cluster).
Once you activated the Snakemake environment, you should perform a dryrun in order to test the executability of your workflow by typing ```snakemake -n```
Snakemake will tell you whether you Snakefile is syntactically correct and whether all required input data is available. The provided example Snakefile has the 'rule all' defined, which tells Snakemake to produce the results for all samples present in the config file.

For upscaling of the data analysis, you can simultaneously run several jobs (execution of rules) via typing ```snakemake --jobs 10```. Please have in mind that each job will need one core and Snakemake will need in addition one core to monitor and submit the jobs (so ```--jobs 10``` will actually need 11 cores).

An example command for Snakemake execution on a cluster would look like this:

```
module load anaconda3/2020.02 #load anaconda package manager
conda activate snakemake  #activate the snakemake environment
module load sratoolkit/2.10.8 #load the sratoolkit in this environment. 
snakemake --cluster \"qsub\" --jobs 3 --use-conda #"qsub" stands for the submission command for a SGE queuing system
```

In our example workflow, the wildcard {sample} is replaced with the individual sample names, followed by universal file extensions. Snakemake will determine dependencies between input and output file based on the user defined rules.


For further information please refer to:  
https://snakemake.github.io  
https://snakemake.readthedocs.io/en/stable/
