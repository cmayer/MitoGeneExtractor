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

## Platforms:
All platforms which provide a C++ compiler.
This program has been tested on Linux and MacOS.

## Installation:


## Basic usage:



## Reference: Please cite when using MitoGeneExtractor
For more information see:
Brasseur ... 2022 ...

Please cite this paper when using MitoGeneExtractor.



