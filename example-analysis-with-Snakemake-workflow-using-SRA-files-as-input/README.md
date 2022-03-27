# Example analysis using a Snakemake workflow and MitoGeneExtractor to search for mitochondrial genes in multiple SRA files:



## Prerequisites:

- Download the following SRA files form NCBI:

https://www.ncbi.nlm.nih.gov/sra/?term=SRR12554982
https://www.ncbi.nlm.nih.gov/sra/?term=SRR12554985

The recommended way to download this data is to use the "prefetch" command from the [SRA toolkidt](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software).
This should result in downloading two files: SRR12554982.sra and SRR12554985.sra.

- Copy the SRA files you downloaded to the current working directory or modify the paths in the Snakefile.

- Create a config.yaml file listing all files that shall be processed. See example file in this folder.

- If raw sequencing reads shall be preprocessed with the aim to remove sequencing adaptors and low quality regions, install a program such as TrimGalore. TrimGalore, which is a wrapper to the cutadapt software is suded in the present workflow for this purpose.

- Install a program that replaces spaces with underscores in fasta sequence names.

