## Use existing Snakemake workflows to easily extract genes from sequencing libraries

MitoGeneExtractor does not require the usage or the installation of Snakemake, but Snakemake provides a convenient way to analyse a large number of data sets. Using Snakemake allows to scale all steps of the analysis, including extraction and quality control steps such as sequence trimming.

Snakemake is a workflow management system which allows upscaling of data analyses in a reproducible way. In our example workflow, the wildcard {sample} is replaced with the individual sample names, followed by universal file extensions. Snakemake will determine dependencies between input and output file based on the user defined rules.

### Installation:
Snakemake relies on Python3, which must be installed. You can install Snakemake in various ways, but the usage of a package manager such as Anaconda is recommended. 
See https://snakemake.readthedocs.io/en/stable/getting_started/installation.html for further information. You can set up a Snakemake environment and install necessary software within the environment in order to minimize potential conflicts with other software packages installed with Anaconda.

### Prerequisites:
Before starting the analyses, the user needs to provide all necessary input data. For the example workflow, you would need to provide the following files:
- Snakefile. (See Snakefile in this folder as an example.)
- The configuration file in .yaml format, in this example called config.yaml. (See config.yaml in this folder as an example.)
- Input raw data, e.g. sequencing reads in FASTQ file format. In principle it is possible to include automatic download steps from NCBI into the workflow, e.g. to download SRA files from the NCBI SRA database.  
- A protein reference for MitoGeneExtractor in FASTA format.
- optional (but needed for the example workflow): cutadapt.yaml file (Anaconda environment for cutadapt safed in .yaml file)  

It is up to the user whether to incorporate analysis steps such as data trimming. Please also have in mind that you might want to adjust the example workflow according to your directory structure. Further, if other versions of the used software are installed, some parameter names might have changed.

### Execution:
There are different ways to execute Snakemake. In particular when computer that allow to use multiple CPU cores, Snakemake provides different ways to make use of these resources, e.g. on high performance clusters.
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

For further information please refer to:  
https://snakemake.github.io  
https://snakemake.readthedocs.io/en/stable/
