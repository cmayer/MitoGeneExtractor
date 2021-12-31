# MitoGeneExtractor

MitoGeneExtractor can be used to conveniently extract mitochondrial genes from sequencing libraries.

Mitochondrial reads are often found as byproduct in sequencing libraries, e.g. from ... .

## How MitoGeneExtractor works:
MitoGeneExtractor aligns all given nucleotide sequences against a
protein reference sequence to obtain a multiple sequence alignment. The
recommended use case is to extract mitochondrial genes from sequencing
libraries, e.g. from hybrid enrichment libraries sequenced on the
Illumina platform. The individual alignments are computed by calling the
exonerate program. 

Exonerate is a very efficient alignment program which allows to align protein and nucleotide sequences.
Nucleotide sequences which cannot be aligned to the protein reference will not be included in
the output. Exonerate should be able to align several 100k short reads in a few
minutes using a single CPU core. Therefore, this approach can be used for projects of any size.

## Supported Platforms:
All platforms are supported for which a C++ compiler is available.
It can be compiled by users without root previliges.
This program has been tested on Linux, HPC Unix platforms and MacOS.

## Installation:
MitoGeneExtractor requires either Exonerate output files as input, or the Exonerate program has to be installed, so that these output files can be generated
when running MitoGeneExtractor. As of writing this, the most recent version of the Exonerate program is 2.4, which is avaiable e.g. here:
https://github.com/nathanweeks/exonerate

On MacOS you will need to install the command line developer tools. For a manual how to do this, see 
https://osxdaily.com/2014/02/12/install-command-line-tools-mac-os-x/.


To install MitoGeneExtractor, do one of the following:
- Clone the MitoGeneExtractor project to your computer. The link can be found by klicking on the "Code" pulldown menu at the top of this page.
- Download the zipped project folder and extract the folder. The link can be found by klicking on the "Code" pulldown menu at the top of this page.  

Now enter the MitoGeneExractor-vx.x folder on the command line and run the make program by typing "make" and hitting return. The make program should be preinstalled on all Linux distributions. On MacOS it is included in the command line developer tools (see above). 
This will generate an executable called MitoGeneExtractor. Either copy this to a directory in your path, to reference it by its full path on the command line.

PUT THIS SOMEWHERE ELSE:
In rare cases, the Exonerate program  quits/crashes unexpectedly.  Interestingly, when rerunning Exonerate, there is a good chance that it will not abort. Why rerunning Exonerate can be successful is a mystery. It has been checked that successful runs always create the same expected output. MitoGeneExtractor can identify, if Exonerate aborted the run and will try upto 10 times by rerunning Exonerate.  The log output of MitoGeneExtractor will report the run count for the succesfull run.

## Get help:
Type MitoGeneExtractor -h (if MitoGeneExtractor is in your path) or Path-to-MitoGeneExtractor/MitoGeneExtractor -h to list all command line options of MitoGeneExtractor. 

In order to run MitoGeneExtractor you need your input sequences in fasta format as well as an amino acid reference sequence of your mitochondrial gene of interest. Example references are included in the reference-sequence-examples folder of this project. For the COI gene, one can specify as a reference the amino acid sequence of the barcode region, or if intended, the full COI sequence. If the full COI sequence shall be extracted, we suggest to create a reference specific for your taxonomic group, since the COI gene can differ considerably in the first and last few amino acids for specific groups with respect to references designed for larger groups. For the barcode region of COI this is normally not a problem. In principle all mitochondrial genes can be extracted with this approach.

## Example analysis:
Let us assume that you want to extract sequences from a short read file called input-reads.fas and that your amino acid reference file is called COI-whole-xxx.fas. Files with these names are found in the "example-analysis" folder.

Run the example analysis with the command (see run-analysis.sh):
xxxx





## Reference: Please cite when using MitoGeneExtractor
For more information see:
Brasseur ... 2022 ...

Please cite this paper when using MitoGeneExtractor.



