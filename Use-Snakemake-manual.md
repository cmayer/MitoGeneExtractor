## Use existing Snakemake workflows to easily extract genes from sequencing libraries

MitoGeneExtractor does not require the usage or the installation of Snakemake, but Snakemake provides a convenient way to analyse a large number of data sets. Using Snakemake allows to scale all steps of the analysis, including extraction and quality control steps such as sequence trimming.

Snakemake is a workflow management system which allows upscaling of data analyses in a reproducible way. In our example workflow, the wildcard {sample} is replaced with the individual sample names, followed by universal file extensions. Snakemake will determine dependencies between input and output file based on the user defined rules.

### Installation:
Snakemake relies on Python3, which must be installed. You can install Snakemake in various ways, but the usage of a package manager such as Anaconda is recommended. 
See https://snakemake.readthedocs.io/en/stable/getting_started/installation.html for further information. You can set up a Snakemake environment and install necessary software within the environment in order to minimize potential conflicts with other software packages installed with Anaconda.
To install snakemake in a anaconda environment, run the following commands:
```
conda create -n snakemake  #this creates the env; the name is arbitrary
conda activate snakemake   #this activates the env; 
conda install snakemake       #this installs snakemake
```

It is up to the user whether to incorporate analysis steps such as data trimming, but it is recommended.
The cutadapt wrapper script TrimGalore! can be used for that,  which is also included in the example Snakemake workflows, but any other trimming software can be used as well.

TrimGalore! can be installed in the snakemake environment via typing 
```
conda install -c bioconda trim-galore
```
or following the installation guide at https://github.com/FelixKrueger/TrimGalore

Further/other software might be required, depending on your type of data. If you intend to mine genes from databases like NCBI, you might need in addition the SRA-toolkit. For your own sequencing data in fastq format, this is not required. 

### Prerequisites:
Before starting the analyses, the user needs to provide all necessary input data. For the example workflow, you would need to provide the following files:
- Snakefile
- The configuration file in .yaml format, in this example called config.yaml
- Input raw data
- A protein reference file for MitoGeneExtractor in FASTA format.

The Snakefile and the config.yaml can be found [here](https://github.com/cmayer/MitoGeneExtractor/tree/last-reviews-before-publication/example-analysis-with-Snakemake-workflow-using-compressed-fastq-files-as-input)

The individual steps executed by snakemake are described [here](https://github.com/cmayer/MitoGeneExtractor/blob/last-reviews-before-publication/example-analysis-with-Snakemake-workflow-using-compressed-fastq-files-as-input/README-snakemake-workflow-for-fastq-files-as-input.md). There, you can als find how to download example input data from NCBI. The input data is assumed to be in the directory ./raw_data/

The protein references for the example workflow can be found [here](https://github.com/cmayer/MitoGeneExtractor/tree/last-reviews-before-publication/example-analysis-with-Snakemake-workflow-using-compressed-fastq-files-as-input). The reference is assumed to be in the directory ./protein_references/

Further references can be found [here](https://github.com/cmayer/MitoGeneExtractor/tree/main/Amino-Acid-references-for-taxonomic-groups)

Everything else will be handled by snakemake.

If you have set up everything, you can execute snakemake as described [here](https://github.com/cmayer/MitoGeneExtractor/blob/last-reviews-before-publication/example-analysis-with-Snakemake-workflow-using-compressed-fastq-files-as-input/README-snakemake-workflow-for-fastq-files-as-input.md).

For further information please refer to:  
https://snakemake.github.io  
https://snakemake.readthedocs.io/en/stable/
