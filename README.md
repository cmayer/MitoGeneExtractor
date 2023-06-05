# MitoGeneExtractor

MitoGeneExtractor can be used to conveniently extract mitochondrial protein-coding genes from next-generation sequencing libraries.
Mitochondrial reads are often found as byproduct in sequencing libraries obtained from whole-genome sequencing, RNA-sequencing or various kinds of reduced representation libraries (e.g. hybrid enrichment libraries).

## List of use recommended use cases:
- Extract mitochondrial protein-coding genes across a broad taxonomic range from sequencing libraries.
  Successfully tested for
  * Illumina short read libraries, namely whole genomic, transcriptomic, reduced representation (e.g. hybrid enrichment libraries, RAD sequencing) libraries.
  * PacBio long read libraries.

- Mine plastome protein-coding genes (matK or rbcL) from sequencing libraries
  Successfully tested for 
  * Illumina short read libraries (genomic and transcriptomic)

- Mine mitochondrial protein-coding genes from transcriptome assemblies. While this sometimes works, we recommended to mine these genes from quality trimmed sequencing reads rather than from the assembly of the reads. We saw example for which a reconstruction from quality trimmed reads was successful while it failed completely when using the assembly. Results might depend on parameters passed to the assembler.

- Mine/excise protein-coding genes from whole mitochondrial genomes, which is often simpler than referring to the annotation (if one is available at all).

- Check for **contamination** in sequencing libraries. Contamination is a common problem. By mining COI sequences from as NGS library, one should obtain COI sequencing reads from the target species as well as from the contaminating species.

- Off label usage: Mine prokaryotic genes from assemblies. (Not tested, but in principle, all genes which can be directly translated into amino acid sequences can be reconstructed with MitoGeneExtractor)

## List of input sources for which MitoGeneExtractor did not work in our tests:
- Long reads from MinION, Oxford Nanopore (Tested on a small number of libraries. See the supplementary materials of the publication for details.) 


## Arguments pro MitoGeneExtractor:

Several tools exist that are able to reconstruct whole or partial mitochondrial genomes from sequencing libraries. All of them extract sequences from assemblies. We found several examples in which assemblies contained strongly reduced amounts of mitochondrial sequence data compared to the raw reads, in particular in the presence of conflicting sequences, e.g. if NUMTs are present or if the library contains DNA from different specimen.
If mitochondrial sequences cannot be assembled, assembly based tools cannot find the genes.

For MGE this means that we recommend to extract protein coding mitochondrial genes from (quality trimmed) reads rather than assemblies if possible.
We have seen examples where the extraction from assemblies worked equally well as the extraction from unassembled reads, but we have also seen cases where the extraction from unassembled reads was successful, but failed when using the assembly.

## How MitoGeneExtractor works - the algorithm:
MitoGeneExtractor aligns all given input nucleotide sequences against a protein reference sequence to obtain a multiple sequence alignment. The intended use case is to extract mitochondrial protein coding genes from sequencing libraries. The individual alignments are computed by calling the Exonerate program. 

Exonerate is a very efficient alignment program which allows to align protein and nucleotide sequences.
Nucleotide sequences which cannot be aligned to the protein reference will not be included in the output. Exonerate should be able to align several 100k short reads in a few minutes using a single CPU core. Therefore, this approach can be used for projects of any size.
Exonerate can align amino acid sequences also to long nucleotide sequences. For this reason, MitoGeneExtractor can also mine sequences from assemblies or from long read libraries. It can even be used to extract genes of interest from whole mitochondrial or nuclear genome/transcriptome assemblies. 

## Input to MitoGeneExtractor

### Required by MitoGeneExtractor
MitoGeneExtractor requires two input files:
- The amino acid reference in fasta file format. For MitoGeneExtractor version 1.9.5 or newer this file can contain multiple protein coding reference genes and/or their variants. This allows to extract all protein coding genes of interest in one program run. 

Many example references are included in the **Amino-Acid-references-for-taxonomic-groups** folder of this project [here](https://github.com/cmayer/MitoGeneExtractor/tree/last-reviews-before-publication/Amino-Acid-references-for-taxonomic-groups). For the COI gene, one can specify as a reference the amino acid sequence of the barcode region, or if intended, the full COI sequence. If the full COI sequence shall be extracted, we suggest to create a reference specific for your taxonomic group, since the COI gene can differ considerably in the first and last few amino acids for specific groups with respect to references designed for larger groups. For the barcode region of COI this is normally not a problem.

- The nucleotide reads/assemblies/genomes in the fasta or fastq format. Since version 1.9.5 any number of fasta or fastq files (e.g. files from paired-end sequencing or multiple replicates) can be specified as program parameters. They will automatically be concatenated and analysed in a single run. Since the paired-end information is not exploited, paired-end libraries can be combined with single-end data.

**Recommendation:** Since quality scores are not used in the analysis, we recommend to pass quality trimmed reads to MitoGeneExtractor.

### Optional input
The user can specify a previously computed vulgar file, i.e. the output file produced by Exonerate.
The vulgar file has to correspond to the input sequences!
Specifying an existing file avoids aligning all reads against the reference(s) again, if only the MGE parameters are changed.

**If you specify a vulgar file name:**
- If the file exists, it will be used.
- If the file does not exist, MGE will run Exonerate to create a new vulgar file and save it using the specified filename.

**If you do not specify a vulgar file name:**
- MGE will run Exonerate to create a new vulgar file and remove it after it has been used.

**Caution:** MGE can only find obvious inconsistencies between the sequence input files and the vulgar file. If the vulgar file contains only partial results (e.g. from a previous run with less data), this will not be noticed and leads to incomplete results.


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

The make program will generate an executable called MitoGeneExtractor_vx.y.z, where x.y.z is the version number. Either copy this to a directory in your path or reference it by its full path on the command line.

<!---
In rare cases, the Exonerate program  quits/crashes unexpectedly.  Interestingly, when rerunning Exonerate, there is a good chance that it will not abort. Why rerunning Exonerate can be successful is a mystery. It has been checked that successful runs always create the same expected output. MitoGeneExtractor can identify, if Exonerate aborted the run and will try up to 10 times by rerunning Exonerate. The log output of MitoGeneExtractor will report the run count for the successful run. --->

## Get help and a full list of command line options:
Type  
```{r, eval=TRUE}
MitoGeneExtractor-vx.y.z -h
```
if MitoGeneExtractor is in your PATH and otherwise
```{r, eval=TRUE}
Path-to-MitoGeneExtractor/MitoGeneExtractor-vx.y.z -h 
```
to display a full list of command line options of MitoGeneExtractor. 

<!---
In order to run MitoGeneExtractor you need your input read data in fasta/fast format as well as reference fasta file with one or more amino acid reference sequence of your mitochondrial gene(s) of interest. 
--->


## Quickstart:
Assume the input file (sequencing reads in fasta format, transcriptome assembly, genome assembly) are stored in the file: query-input.fas.
Furthermore assume that the amino acid reference sequence is stored in the COI-reference.fas file.
Then the following command could be used to attempt to reconstruct the COI sequence from the read data in the query-input.fas file:

```{r, eval=TRUE}
MitoGeneExtractor-vx.y.z  -d query-input.fas -p COI-reference.fas -V vulgar.txt -o out-alignment.fas -n 0 -c out-consensus.fas -t 0.5 -r 1 -C 2
```
Specifying the name of the vulgar file is optional, but recommended as this is the most-time consuming step. If the file exists, it is used as input instead of calling Exonerate to create it. If it does not exist, the name is used to create the vulgar file. The -C 2 option specifies the genetic code (here: vertebrate mitochondrial), the -t 0.5 option specifies the consensus threshold and the -r 1 and -n 0 options are used for a stricter alignment quality (see options for details). 

If your read data is in fastq format, you could run the same analysis via this command:
```{r, eval=TRUE}
MitoGeneExtractor-vx.y.z  -d query-input.fq -p COI-reference.fas -V vulgar.txt -o out-alignment.fas -n 0 -c out-consensus.fas -t 0.5 -r 1 -C 2
```
If you have multiple input files (e.g. paired-end data (PE) and single-end (SE) data) you cand specify this as follows:
```{r, eval=TRUE}
MitoGeneExtractor-vx.y.z  -d PE_query-input_1.fq PE_query-input_2.fq SE_query-input.fq -p COI-reference.fas -V vulgar.txt -o out-alignment.fas -n 0 -c out-consensus.fas -t 0.5 -r 1 -C 2
```
Note, that the order of file names does not matter. It is also possible to simultaneously specify input data in fastq and fasta format.


## Example analysis:
An example analysis for the MitoGeneExtractor program can be found in the **example-analysis-for-MitoGeneExtractor** folder [here](https://github.com/cmayer/MitoGeneExtractor/tree/last-reviews-before-publication/example-analysis-for-MitoGeneExtractor). The Readme.md file in this folder provided the necessary information to run the example analysis and provides further details.


## Prepared Snakemake workflows
You can find a description how data preprocessing and MitoGeneExtractor analyses can be implemented in Snakemake [here](https://github.com/cmayer/MitoGeneExtractor/blob/last-reviews-before-publication/Use-Snakemake-manual.md)

An example snakemake analysis can be found [here](https://github.com/cmayer/MitoGeneExtractor/tree/last-reviews-before-publication/example-analysis-with-Snakemake-workflow-using-compressed-fastq-files-as-input)


## Command line options:
A full list of the command line options is available when typing
MitoGeneExtractor-vx.y.z -h

**-d <string>,  --dna_fasta_file <string>  (accepted multiple times)**   
 Specifies the input query nucleotide sequence files in the fasta format. Sequences are expected not to include gap characters. This option can be specified multiple times if multiple input files shall be analysed in one run. If sequence files contain reads, they should have been quality filtered before being used as input for this program. This option can be combined with multiple input files in the fastq format (see -q option).

**-q <string>,  --dna_fastq_file <string>  (accepted multiple times)**  
 Specifies the input query nucleotide sequence files in the fastq format. This option can be specified multiple times if multiple input files shall be analysed in one run. All input files will be converted to a fasta file without taking into account the quality scores. Sequence files should be quality filtered before being used as input for this program. This option can be combined with multiple input files in the fasta format (see -d option).

**-p <string>,  --prot_reference_file <string>**  
 Specifies the fasta file containing the amino acid reference sequences. This file can contain one or multiple reference sequences. All input nucleotide sequences are aligned against all references. Hits with a score higher than the minimum are considered. If a sequence matches multiple reference genes/variants, the sequence will be assigned to the reference for which the alignment score is higher or to both if the scores are equal. 

**-o <string>,  -- <string> (required)**  
Specifies the base name of alignment output file(s). Aligned input sequences are written to a file with the name: BaseName + sequenceNameOfRefernce + .fas for each reference sequence.

**-V <string>,  --vulgar_file <string>**  
 Specifies the name of Exonerate vulgar file. If the specified file exists, it will be used for the analysis. If it does not exist MitoGeneExtractor will run Exonerate in order to create the file with this name. The created file will then be used to proceed. If no file is specified with this option, a temporary file called tmp-vulgar.txt will be created and removed after the program run. In this case a warning will be printed to the console.

**-e <string>,  --exonerate_program <string>**  
 Specifies the name of the Exonerate program in system path OR the path to the Exonerate program including the program name. Default: Exonerate

**-V <string>, --vulgar_file <string>**   Name of Exonerate vulgar file. If the specified file exists, it will be used for this analysis. If it does not exist, MitoGeneExtractor will run Exonerate in order to create the file. The created file will then be used to proceed. If no file is specified with this option, a temporary file called tmp-vulgar.txt will be created and removed after the program run. In this case a warning will be printed to the console, since the vulgar file cannot be used again. (Optional, but recommended parameter) 

**-e <string>, --exonerate_program <string>**   Name of the Exonerate program in the system path OR the path to the Exonerate program including the program name. Default: exonerate. (Optional parameter)

**-n <int>,  --numberOfBpBeyond <int>**  
 Specifies the number of base pairs that are shown beyond the Exonerate alignment. A value of 0 means that the sequence is clipped at the point the Exonerate alignment ends. Values >0 can lead to the inclusion of sequence segments that do not align well with the amino acid sequence and have to be treated with caution. They might belong to chimera, NUMTs, or other problematic sequences. Larger values might be included e.g. if problematic sequences with a well matching seed alignment are of interest. CAUTION: Bases included with this option might not be aligned well or could even belong to stop codons! They should be considered as of lower quality compared to other bases. Bases that are added with this option are added as lower case characters to the output alignment file. A sequence coverage of bases not belonging to these extra bases can be requested with the --minSeqCoverageInAlignment_uppercase option. Default: 0.

**-c <string>,  --consensus_file <string> (required)**
Specifies the base name of the consensus sequence output file(s). A consensus sequence with the name baseName + reference-sequence-name + .fas  is written for each reference sequence.

**-t <float>,  --consensus_threshold <float>**  
 This option modifies the consensus threshold. Default: 0.5 which corresponds to 50%.

**-D,  --includeDoubleHits**  
 Include input sequences with two alignment results against the same reference.

**--noGaps**  
 Do not include reads for which the alignment with the reference contains gaps.

**-g,  --onlyGap**  
 Include only reads which aligned with a gap.

**--report_gaps_mode <int>**  
 Gaps can be reported in different ways. With this option the reporting mode can be specified: 1: report leading and trailing gaps with '-' character. Report internal gaps (introduced with options -G or -g) with '~' character. 2: report leading and trailing gaps with '-' character. Report internal gaps (introduced with options -G or -g) with '-' characters. 3: Remove all gap characters in output. In this case sequences are extracted but are reported with respect to the reference. Default: 1.

**-f <int>,  --frameshift_penalty <int>**  
 The frameshift penalty passed to Exonerate. Default: -9. Higher values lead to lower scores and by this can have the following effects: (i) hit regions are trimmed since trimming can lead to a better final alignment score, (ii) they can also lead to excluding a read as a whole if the final score is too low and trimming does lead to a higher score. The default of the Exonerate program is -28. A value of -9 (or other values lower than -28) lead to more reads in which the best alignment has a frameshift. In order to remove reads that do not align well, one can use a smaller frameshift penalty and then exclude hits with a frameshift, see -F option).

**-C <int>,  --genetic_code <int>**  
 The number of the genetic code to use in Exonerate, if this step is required. See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for details. Default: 2, i.e. vertebrate mitochondrial code.

**-s <int>,  --minExonerateScoreThreshold <int>**  
 The score threshold passed to Exonerate to decide whether to include or not include the hit in the output.

**-r <float>,  --relative_score_threshold <float>**  
 Specified the relative alignment score threshold for Exonerate hits to be considered. The relative score is the score reported by Exonerate divided by the alignment length. Default 1. Reasonable thresholds are between 0.7 and 2.0.

**--minSeqCoverageInAlignment_total <int>**  
 Specifies the absolute value of the minimum alignment coverage for computing the consensus sequence. For the coverage, all nucleotides count, also lower case nucleotides that have been added beyond the Exonerate alignment region. Default: 1. Increasing this value increases the number of unknown nucleotides in the consensus sequence.

**--minSeqCoverageInAlignment_uppercase <int>**  
 Specifies the absolute value of the minimum alignment coverage for computing the consensus sequence. As coverage, only upper case nucleotides are taken into account, i.e. no nucleotides are counted that have been added beyond the Exonerate alignment region. Bases beyond the Exonerate alignment are added with the -n or --numberOfBpBeyond option. If no bases are added beyond the Exonerate alignment (default), the effect of this option is identical to the minSeqCoverageInAlignment_total option. Default: 1. Increasing this value increases the number of unknown nucleotides in the consensus sequence.

**--temporaryDirectory <string>**  
 MGE has to create potentially large temporary files, e.g. if multiple input files are specified, or if fastq file are specified. With this option these files will not be created in the directory the program was launched, but in the specified tmp directory. 

**--treat-references-as-individual**  
Input sequences which can be aligned with different reference sequences are by default assigned only to the references for which the alignment score is equal to the best score achieved by this input sequence. This score competition is switched off if this option is specified. This treats multiple references as if they are specified in independent program runs. 

**--keep-concat-input-file**  
If multiple input files are specified, MGE first creates a concatenated file. By default this file is removed. Use this option if you want to keep this file.

**--verbosity <int>**  
 Specifies how much run time information is printed to the console. Values: 0: minimal output, 1: important notices, 2: more notices, 3: basic progress, 4: detailed progress, 50-100: debug output, 1000: all output.

**--,  --ignore_rest**  
 Ignores the rest of the labeled arguments following this flag.

**--version**  
Displays version information and exits.

**-h,  --help**  
Displays usage information and exits.


<!---
## Applications:
- Extract COI and other protein coding mitochondrial genes in a sequencing library or transcriptome. 

**Strategy**  
Provide amino acid references for all genes of interest in the reference file and run MitoGeneExtractor.

- Investigate the amount of contamination in a sequencing library.

**Strategy**  
 Provide the amino acid references for COI of your target group and potentially distantly related contamination. 
Contamination of closely related taxa will show up as multiple variants in the alignment file. Contamination of distantly related taxa can be found as hits to distantly related COI sequence.

--->

## Project outlook:

Currently, we are exploring the utility of using HMMs, namely nhmmer as another option and alternative to exonerate.

## Authors of the publication:
- Marie Brasseur, ZFMK/LIB, Bonn, Germany
- Jonas Astrin, ZFMK/LIB, Bonn, Germany
- Matthias Geiger, ZFMK/LIB, Bonn, Germany
- Christoph Mayer, ZFMK/LIB, Bonn, Germany

## Authors of the software project:
- Christoph Mayer, ZFMK/LIB, Bonn, Germany: MitoGeneExtractor program.
- Marie Brasseur, ZFMK/LIB, Bonn, Germany: Snakemake pipeline and analyses for publication.

## Reference: When using MitoGeneExtractor please cite:
Brasseur, M.V., Astrin, J.J., Geiger, M.F., Mayer, C., 2023. MitoGeneExtractor: Efficient extraction of mitochondrial genes from next-generation sequencing libraries. Methods in Ecology and Evolution.
[https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14075](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14075)

Since MitoGeneExtractor uses the Exonerate program, please also cite:

Slater, G.S.C., Birney, E., 2005. Automated generation of heuristics for biological sequence comparison. BMC Bioinformatics 6, 31. https://doi.org/10.1186/1471-2105-6-31





