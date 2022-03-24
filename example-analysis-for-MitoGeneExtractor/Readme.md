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

In order not to exceed total file size limits of github and still be able to include the example file in the repository, the example input file SRR12554985_trimmed_reduced.fas.zip was created by extracting from the original file exactly those reads that match with the COI sequence.
Still, wenn running the example here, the sequences are extracted and aligned.
You can also recreate the full fasta file by downloading the the SRR12554985.sra file and following Brasseur et al. 2022 to prepare the fasta file.

## Running the analysis:
./runner.sh

This runs the analysis twice, for two amino acid reference sequences for COI. (i) with the general vertebrate reference for COI and (ii) with the specific COI reference for Passeriformes. The two commands read as follows:

../MitoGeneExtractor-v1.9.1 -d SRR12554985_trimmed_reduced.fas -p ../Amino-Acide-references-for-tanomic-groups/xxx-need-to-be-added.fasta -V vulgar-SRR12554985_PasseriformesReference.txt -o SRR12554985_align_PasseriformesReference.fas -n 0 -c SRR12554985_cons_PasseriformesReference.fas -t 0.5 -r 1 -C 2

../MitoGeneExtractor-v1.9.1 -d SRR12554985_trimmed_reduced.fas -p ../Amino-Acide-references-for-tanomic-groups/COI-vertebrata-protein-consensus-50_whole.fasta -V vulgar-SRR12554985_vertebrateReference.txt -o SRR12554985_align_vertebrateReference.fas -n 0 -c SRR12554985_cons_vertebrateReference.fas -t 0.5 -r 1 -C 2

Note that if you move the example folder, the path to the MitoGeneExtractor-v1.9.1 has to be changed accordingly.

The options used in this example are descrived in detail in the main readme or when calling the program with the option -h:
../MitoGeneExtractor-v1.9.1 -h

The the results files for the consensus sequences are almost identical. ....

