# Example analysis using a Snakemake workflow and MitoGeneExtractor to search for mitochondrial genes in multiple SRA files:

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

This should result in downloading the two files: SRR12554982.sra and SRR12554985.sra.
Copy the SRA files you downloaded to the current working directory or modify the paths in the Snakefile.

- The configuration file in .yaml format. See config.yaml file in this folder.

- A protein reference for MitoGeneExtractor in FASTA format

- (Optional but recommended. Needed in this example workflow.) If raw sequencing reads shall be preprocessed with the aim to remove sequencing adaptors and low quality regions, install a program such as TrimGalore. TrimGalore, which is a wrapper to the cutadapt software is used in the present workflow for this purpose.
<!-- Expert users can use the cutadap.yaml file to install  file (Anaconda environment for cutadapt safed in .yaml file)  -->

- A program that replaces spaces with underscores in fasta sequence names.

This example is work in progress. See updates soon.