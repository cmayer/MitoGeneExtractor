# MitoGeneExtractor

MitoGeneExtractor can be used to conveniently extract mitochondrial protein-coding genes from next-generation sequencing libraries.
Mitochondrial reads are often found as byproduct in sequencing libraries obtained from whole-genome sequencing, RNA-sequencing or various kinds of reduced representation libraries (e.g. hybrid enrichment libraries).

## List of use recommended use cases:
- Extract mitochondrial protein-coding genes across a broad taxonomic range from sequencing libraries.
  Successfully tested for
  * Illumina short read libraries, whole genomic, transcriptomic, reduced representation (hybrid enrichment libraries, RAD sequencing)
  * PacBio long read libraries.

- Mine plastome protein-coding genes (matK or rbcL) from sequencing libraries
  Successfully tested for 
  * (??) Illumina short read libraries

- Mine mitochondrial protein-coding genes from transcriptome assemblies. While this sometimes works, we recommended to mine these genes from quality trimmed sequencing reads rather than from the assembly of the reads. We saw example for which a reconstruction from quality trimmed reads was successful while it failed completely when using the assembly. Results might depend on parameters passed to the assembler.

- Mine/excise protein-coding genes from whole mitochondrial genomes, which is often simpler than referring to the annotation (if one is available at all).

- Check for **contamination*+ in sequencing libraries. Contamination is a common problem. By searching for COI sequences in taxonomic group, one should obtain COI sequencing reads from the target species as well as from the contaminating species. 

- Off label usage: Mine prokaryotic genes from assemblies. (Not tested, but in principle, all genes which can be directly translated into amino acid sequences can be reconstructed with MitoGeneExtractor)

## List of input sources for which did not work in our tests:
- Long reads from MinION, Oxford Nanopore (Tested on a small number of libraries. See the supplementary materials of the publication for details.) 


## Arguments pro MitoGeneExtractor:

Several tools exist that are able to reconstruct whole or partial mitochondrial genomes from sequencing libraries. Most of them extract genes from assemblies. We found several examples in which assemblies contained strongly reduced amounts of mitochondrial sequences, in particular in the presence of conflicts sequences in the region of interest, e.g. if Numts are present or if the library contains DNA from different specimen.
If mitochondrial genes cannot be assembled, assembly based tools cannot find these genes.

For MGE this means that we recommend to extract protein coding mitochondrial genes from (quality trimmed) reads rather than assemblies if possible.
We have seen examples where the extraction from assemblies worked equally well as the extraction from unassembled reads, but we have also seen cases where the extraction from unassembled reads was successful, but failed for the assembly.

## How MitoGeneExtractor works - the algorithm:
MitoGeneExtractor aligns all given input nucleotide sequences against a protein reference sequence to obtain a multiple sequence alignment. The intended use case is to extract mitochondrial protein coding genes from sequencing libraries. The individual alignments are computed by calling the Exonerate program. 

Exonerate is a very efficient alignment program which allows to align protein and nucleotide sequences.
Nucleotide sequences which cannot be aligned to the protein reference will not be included in the output. Exonerate should be able to align several 100k short reads in a few minutes using a single CPU core. Therefore, this approach can be used for projects of any size.
Exonerate can align amino acid sequences also to long nucleotide sequences. For this reason, MitoGeneExtractor can also mine sequences from assemblies or from long read libraries. It can even be used to extract genes of interest from whole mitochondrial or nuclear genome/transcriptome assemblies. 

## Input

### Required by MitoGeneExtractor
MitoGeneExtractor requires two input files:
- The amino acid reference in fasta file format. For MitoGeneExtractor version 1.9.5 or newer this file can contain multiple protein coding reference genes and/or their variants. This allows to extract all protein coding genes of interest in one program run. 

- The nucleotide reads/assemblies/genomes in the fasta or fastq format. Since version 1.9.5 any number of fasta or fastq files (e.g. files from paired-end sequencing or multiple replicates) can be specified as program parameters. They will automatically be concatenated and analysed in a sigle run. Since the paired-end information is not exploited, paired-end libraries can be combined with single-end data.

***Recommendation:*** Since quality scores are not used in the analysis, we recommend to pass quality trimmed reads to MitoGeneExtractor.

### Optional input
The user can specify the vulgar file, i.e. the output file produced by Exonerate, that corresponds to the input sequences.
This avoids aligning all reads against the reference(s) again, if only the MGE parameters are changed.

** Caution: ** MGE can only find obvious inconsistencies between the sequence input files and the vulgar file. If the vulgar file contains only partial results, this will not be noticed and leads to incomplete results.


## Supported Platforms:
MitoGeneExtractor: All platforms are supported for which a C++ compiler is available.
It can be compiled by users without root privileges.
This program has been tested on Linux, HPC Unix platforms and MacOS. It should be possible in principle to compile it on windows. However, there is not support from the authors for this platform.

Snakemake workflow: Linux, HPC Unix platforms and MacOS. It should be possible to run this on Windows systems but the authors have no hardware to test this.

## Installation:
MitoGeneExtractor requires either an Exonerate output file as input, or the Exonerate program has to be installed, so that such an Exonerate output can be generated by MitoGeneExtractor. As of writing this, the most recent version of the Exonerate program is 2.4, which is available e.g. here:
https://github.com/nathanweeks/exonerate

On MacOS you will need to install the command line developer tools. For a manual how to do this, see 
https://osxdaily.com/2014/02/12/install-command-line-tools-mac-os-x/.


To install MitoGeneExtractor, do one of the following:
- Clone the MitoGeneExtractor project to your computer. The link can be found by clicking on the "Code" pulldown menu at the top of this page.
- Download the zipped project folder and extract the folder. The link can be found by clicking on the "Code" pulldown menu at the top of this page.  

Now enter the MitoGeneExractor-vx.x folder on the command line and run the make program by typing "make" and hitting return. The make program should be preinstalled on all Linux distributions. On MacOS it is included in the command line developer tools (see above). 

```{r, eval=TRUE}
cd/MitoGeneExractor-vx.x
make
```

The make program will generate an executable called MitoGeneExtractor_. Either copy this to a directory in your path or reference it by its full path on the command line.

<!---
In rare cases, the Exonerate program  quits/crashes unexpectedly.  Interestingly, when rerunning Exonerate, there is a good chance that it will not abort. Why rerunning Exonerate can be successful is a mystery. It has been checked that successful runs always create the same expected output. MitoGeneExtractor can identify, if Exonerate aborted the run and will try up to 10 times by rerunning Exonerate. The log output of MitoGeneExtractor will report the run count for the successful run. --->

## Get help and a full list of command line options:
Type  
```{r, eval=TRUE}
MitoGeneExtractor -h
```
if MitoGeneExtractor is in your PATH and otherwise
```{r, eval=TRUE}
Path-to-MitoGeneExtractor/MitoGeneExtractor -h 
```
to display a full list of command line options of MitoGeneExtractor. 

In order to run MitoGeneExtractor you need your input read data in fasta/fast format as well as reference fasta file with one or more amino acid reference sequence of your mitochondrial gene(s) of interest. Example references are included in the **Amino-Acid-references-for-taxonomic-groups** folder of this project [here](https://github.com/cmayer/MitoGeneExtractor/tree/last-reviews-before-publication/Amino-Acid-references-for-taxonomic-groups). For the COI gene, one can specify as a reference the amino acid sequence of the barcode region, or if intended, the full COI sequence. If the full COI sequence shall be extracted, we suggest to create a reference specific for your taxonomic group, since the COI gene can differ considerably in the first and last few amino acids for specific groups with respect to references designed for larger groups. For the barcode region of COI this is normally not a problem. 

## Example analysis:
An example analysis for the MitoGeneExtractor program can be found in the **example-analysis-for-MitoGeneExtractor** folder [here](https://github.com/cmayer/MitoGeneExtractor/tree/last-reviews-before-publication/example-analysis-for-MitoGeneExtractor). The Readme.md file in this folder provided the necessary information to run the example analysis and provides further details.

***Quickstart:***
Assume the input file (sequencing reads in fasta format, transcriptome assembly, genome assembly) are stored in the file: query-input.fas.
Furthermore assume that the amino acid reference sequence is stored in the COI-reference.fas file.
Then the following command could be used to attempt to reconstruct the COI sequence from the read data in the query-input.fas file:

```{r, eval=TRUE}
MitoGeneExtractor  -d query-input.fas -p COI-reference.fas -V vulgar.txt -o out-alignment.fas -n 0 -c out-consensus.fas -t 0.5 -r 1 -C 2
```
Specifying the name of the vulgar file is optional, but recommended as this is the most-time consuming step. If the file exists, it is used as input instead of calling Exonerate to create it. If it does not exist, the name is used to create the vulgar file. The -C 2 option specifies the genetic code (here: vertebrate mitochondrial), the -t 0.5 option specifies the consensus threshold and the -r 1 and -n 0 options are used for a stricter alignment quality (see options for details). 

If your read data is in fastq format, you could run the same analysis via this command:
```{r, eval=TRUE}
MitoGeneExtractor  -d query-input.fq -p COI-reference.fas -V vulgar.txt -o out-alignment.fas -n 0 -c out-consensus.fas -t 0.5 -r 1 -C 2
```
If you have multiple input files (e.g. paired-end data (PE) and single-end (SE) data) you cand specify this as follows:
```{r, eval=TRUE}
MitoGeneExtractor  -d PE_query-input_1.fq PE_query-input_2.fq SE_query-input.fq -p COI-reference.fas -V vulgar.txt -o out-alignment.fas -n 0 -c out-consensus.fas -t 0.5 -r 1 -C 2
```
Note, that the order of file names does not matter.

## Prepared Snakemake workflows
You can find a description how data preprocessing and MitoGeneExtractor analyses can be implemented in Snakemake [here](https://github.com/cmayer/MitoGeneExtractor/blob/last-reviews-before-publication/Use-Snakemake-manual.md)

A Snakefile which starts with .sra files as input can be found [here](https://github.com/cmayer/MitoGeneExtractor/tree/last-reviews-before-publication/example-analysis-with-Snakemake-workflow-using-SRA-files-as-input)

A Snakefile which starts with .fastq data can be found [here](https://github.com/cmayer/MitoGeneExtractor/tree/last-reviews-before-publication/example-analysis-with-Snakemake-workflow-using-compressed-fastq-files-as-input)


## Command line options:
A full list of the command line options is available when typing
MitoGeneExtractor -h

**-d <string>, --dna_sequences_file <string>:** Name (potentially including the path) of the nucleotide sequence file in the fasta format. Sequences are expected to be unaligned without gaps. Typically, these are short or long reads but could also be assembled fragments of any length. (Required parameter) 

**-p <string>, --prot_reference_file <string>:** Protein sequence file in the fasta format. This is the sequence used to align the reads against. File is expected to have exactly one reference sequence. (Required parameter) 

**-o  <string>** Name of alignment output file. (Required parameter) 

**-V <string>, --vulgar_file <string>:** Name of Exonerate vulgar file. If the specified file exists, it will be used for this analysis. If it does not exist, MitoGeneExtractor will run Exonerate in order to create the file. The created file will then be used to proceed. If no file is specified with this option, a temporary file called tmp-vulgar.txt will be created and removed after the program run. In this case a warning will be printed to the console, since the vulgar file cannot be used again. (Optional, but recommended parameter) 

**-e <string>, --exonerate_program <string>:** Name of the exonerate program in the system path OR the path to the exonerate program including the program name. Default: exonerate. (Optional parameter)

**-n <int>, --numberOfBpBeyond <int>:** Specifies the number of base pairs that are shown beyond the Exonerate alignment. A value of 0 means that the sequence is clipped at the point the Exonerate alignment ends. Values of 1 and 2 make sense, since exonerate does not consider partial matches of the DNA to the amino acid sequence, so that partial codons would always be clipped, even if the additional base pairs would match with the expected amino acid. Values >0 lead to the inclusion of sequence segments that do not align well with the amino acid sequence and have to be treated with caution. They might belong to chimera, NUMTs, or other problematic sequences. Larger values might be included e.g. if problematic sequences with a well matching seed alignment are of interest. CAUTION: Bases included with this option might not be aligned well or could even belong to stop codons! They should be considered as of lower quality compared to other bases. Bases that are added with this option are added as lower case characters to the output alignment file. A sequence coverage of bases not belonging to these extra bases can be requested with the --minSeqCoverageInAlignment_uppercase option. Default: 0. Type: integer. (Optional parameter)

**-c <string>, --consensus_file <string>:** If this option is specified, a consensus sequence of all aligned reads is written to the file with the specified name. Normally, this is the intended output. Default: No consensus is written, since no good default output file is known. (Optional parameter)

**-t <float>, --consensus_threshold <float>:** This option modifies the consensus threshold. Default: 0.5 which corresponds to 50%. Type: Decimal number. (Not required)

**-D, --includeDoubleHits:** Include reads with two alignment results found by exonerate. Default: No. (Optional parameter)

**-g, --onlyGap:** Include only reads which aligned with a gap. Useful for finding problems in the set of reads or the correspondence between the reference and the reads.

**--noGaps:** Do not include reads which aligned with a gap. Default: No. (Optional parameter)

**--report_gaps_mode <int>:** Gaps can be reported in different ways. With this
     option the reporting mode can be specified: 1: report leading and
     trailing gaps with '-' character. Report internal gaps (introduced
     with options -G or -g) with '~' character. 2: report leading and
     trainling gaps with '-' character. Report internal gaps (introduced
     with options -G or -g) with '-' characters. 3: Remove all gap
     characters in output. In this case sequences are extracted but are
     reported with respect to the reference. Default: 1. (Optional parameter)

**-f <int>, --frameshift_penalty <int>** The frameshift penalty passed to exonerate. The option value has to be a negative integer. Default: -9. More negative values lead to lower scores and by this can have the following effects: (i) hit regions are trimmed since trimming can lead to a better final alignment score, (ii) they can also lead to excluding a read as a whole if the final score is too low and trimming does lead to a higher score. The default of the exonerate program is -28. A value of -9 (or other values less negative than -28) lead to more reads in which the best alignment has a frameshift. In order to remove reads that do not align well, one can use a less negative value for the frameshift penalty and then exclude hits with a frameshift, see -F option). (Optional parameter)

**-C <int>, --genetic_code <int>:** The number of the genetic code to use in Exonerate, if this step is required. Default: 2, i.e. vertebrate mitochondrial code. (Type: integer). [See genetic code list at NCBI.](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). (Not required but recommended). ** A misspecification of the genetic code leads to unusable results. Make sure the default is correct, or specify the genetic code.** (Optional parameter)

**-r <float>, --relative_score_threshold <float>:** Specified the relative alignment score threshold for Exonerate hits to be considered. The relative score is the score reported by Exonerate divided by the alignment length. Default 1. Reasonable thresholds are between 0.7 and 2.0.  (Optional parameter)

**--minSeqCoverageInAlignment_total <int>:** Specifies the absolute value of the minimum alignment coverage for computing the consensus sequence. For the coverage, all nucleotides count, also those lower case nucleotides that have been added beyond the exonerate alignment region. Default: 1. Increasing this value increases the number of unknown nucleotides in the consensus sequence.  (Optional parameter)

**--minSeqCoverageInAlignment_uppercase <int>:** Specifies the absolute value of the minimum alignment coverage for computing the consensus sequence. As coverage, only upper case nucleotides are taken into account, i.e. no nucleotides are counted that have been added beyond the Exonerate alignment region. Bases beyond the Exonerate alignment are added with the -n or --numberOfBpBeyond option. If no bases are added beyond the Exonerate alignment (default), the effect of this option is identical to the minSeqCoverageInAlignment_total option. Default: 1. Increasing this value increases the number of unknown nucleotides in the consensus sequence.  (Optional parameter)

**-s <int>, --minExonerateScoreThreshold <int>:** The score threshold passed to exonerate to decide whether to include or not include the hit in the output. Typ: integer, optional parameter.  (Optional parameter)

**--verbosity <int>:** Specifies how much run time information is printed to the console. Values: 0: minimal output, 1: important notices, 2: more notices, 3: basic progress, 4: detailed progress, 50-100: debug output, 1000: all output. Default: 1. (Optional parameter)


## Applications:
- Extract COI and other protein coding mitochondrial genes in a sequencing library or transcriptome. 

**Strategy:** Provide amino acid references for all genes of interest in the reference file and run MitoGeneExtractor.

- Investigate the amount of contamination in a sequencing library.

**Strategy:** Provide the amino acid references for COI of your target group and potentially distantly related contamination. 
Contamination of closely related taxa will show up as multiple variants in the alignment file. Contamination of distantly related taxa can be found as hits to distantly related COI sequence.

## Project outlook:

Currently, we are exploring the utility of using HMMs, namely nhmmer as another option and alternative to exonerate.

## Authors of the publication:
- Marie Brasseur, ZFMK, Bonn, Germany
- Jonas Astrin, ZFMK, Bonn, Germany
- Matthias Geiger, ZFMK, Bonn, Germany
- Christoph Mayer, ZFMK, Bonn, Germany

## Authors of the software project:
- Christoph Mayer, ZFMK, Bonn, Germany: MitoGeneExtractor program.
- Marie Brasseur, ZFMK, Bonn, Germany: Snakemake pipeline and analyses for publication.

## Reference: Please cite when using MitoGeneExtractor
Brasseur ... 2023 ...




