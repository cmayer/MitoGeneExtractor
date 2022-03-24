# Readme for the example analysis of the MitoGeneExtractor program (without Snakemake workflow).

## Prerequisites
- This example requires that the exonerate program is installed and that it can be found in the system path.
If you exonerate program should not be in the system path, the runner.sh file can be modified by adding the
option -e /path-to-exonerate-program/exonerate-program-name so that MitoGeneExtractor can find it.

- Mito has to be compiled. Enter the source directory mitogeneextractor and run the command "./make".

## The example data
One of the sequencing libraries analysed in Brasseur et al. 2022 was:
https://www.ncbi.nlm.nih.gov/sra/?term=SRR12554985
This file has been downloaded and extracted as described in Brasseur et al. 2022 using the NCBI SRA-tools 
version 1.10, it has been preprocessed and converted to fasta as described in Brasseur et al. 2022. 
The resulting fasta file has been analysed with the first command of the runner.sh file.

Run the analysis:
./runner.sh

This runs the analysis with two amino acid reference sequences for COI. The general vertebrate reference for COI and the specific COI reference for Passeriformes.
The results are almost identical. ....

