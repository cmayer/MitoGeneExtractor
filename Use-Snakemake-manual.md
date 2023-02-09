## Use existing Snakemake workflows to easily extract genes from sequencing libraries

MitoGeneExtractor does not require the usage or the installation of Snakemake, but Snakemake provides a convenient way to analyse a large number of data sets. Using Snakemake allows to scale all steps of the analysis, including extraction and quality control steps such as sequence trimming.

Snakemake is a workflow management system which allows upscaling of data analyses in a reproducible way. In our example workflow, the wildcard {sample} is replaced with the individual sample names, followed by universal file extensions. Snakemake will determine dependencies between input and output file based on the user defined rules.

### Installation:
Snakemake relies on Python3, which must be installed. You can install Snakemake in various ways, but the usage of a package manager such as Anaconda is recommended. 
See https://snakemake.readthedocs.io/en/stable/getting_started/installation.html for further information. You can set up a Snakemake environment and install necessary software within the environment in order to minimize potential conflicts with other software packages installed with Anaconda.
To install snakemake in a anaconda environment, run the following commands:
```
conda create -n snakemakeenv  #this creates the env; the name is arbitrary
conda activate snakemakeenv   #this activates the env; 
conda install snakemake       #this installs snakemake
```

### Prerequisites:
Before starting the analyses, the user needs to provide all necessary input data. For the example workflow, you would need to provide the following files:
- Snakefile
- The configuration file in .yaml format, in this example called config.yaml
- Input raw data, e.g. sequencing reads in FASTQ file format   
- A protein reference file for MitoGeneExtractor in FASTA format 

It is up to the user whether to incorporate analysis steps such as data trimming, but it is recommended.
The cutadapt wrapper script TrimGalore! can be used for that,  which is also included in the example Snakemake workflows, but any other trimming software can be used as well.
Please also have in mind that you might want to adjust the example workflow according to your directory structure. Further, if other versions of the used software are installed, some parameter names might have changed.


### Execution:
There are different ways to execute Snakemake. In particular when computer that allow to use multiple CPU cores, Snakemake provides different ways to make use of these resources, e.g. on high performance clusters.
Once you activated the Snakemake environment, you should perform a dryrun in order to test the executability of your workflow by typing ```snakemake -n```
Snakemake will tell you whether you Snakefile is syntactically correct and whether all required input data is available. The provided example Snakefile has the 'rule all' defined, which tells Snakemake to produce the results for all samples present in the config file.

For upscaling of the data analysis, you can simultaneously run several jobs (execution of rules) or provide multiple cores.
For running several jobs simultaneosly, type ```snakemake --jobs 10``` 
Several cores can be provided via ```snakemake --cores 10``` and Snakemake will scale the analysis according to the specifications in the Snakefile. Please have in mind that each job will need one core and Snakemake will need in addition one core to monitor and submit the jobs (so ```--jobs 10``` will actually need 11 cores).

An example command for running Snakemake on 10 cores could look like this:
```
conda activate snakemakeenv  #activate the snakemake environment
snakemake --cores 10 --snakefile Snakefile_example
```
If each rule in the workflow Snakefile_example claim only one core, then Snakemake scales the analysis accordingly, i.e. 10 jobs will be executed simultaneosly.

For further information please refer to:  
https://snakemake.github.io  
https://snakemake.readthedocs.io/en/stable/
