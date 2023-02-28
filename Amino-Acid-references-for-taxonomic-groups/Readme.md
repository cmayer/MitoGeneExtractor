# Amino acid reference sequences

## Idea: 
The idea behind MitoGeneExtractor is to align the input sequences (typically reads) with the amino acid sequences of protein coding genes of interest. All input sequences that can be aligned successfully are added to the alignment. Finally a consensus sequence is determined. Currently, the Exonerate program is used to align the input sequences with the amino acid references.

## Reason we use an amino acid sequence as a reference:

For mitochondrial genes one finds almost no length altering mutations and small amounts of sequence variation in the amino acid sequence even for large taxonomic groups. 

This makes the amino acid sequences an ideal reference for extracting these genes. 

We even found that general reference sequences for vertebrates and arthropod can be used to extract large parts of the mitochondrial genes for all taxa from these group we have tested. Only the slightly less conserved beginning and end of these genes require more specific references.

In arthropods, some taxa have short deletion in the COI sequence sequences and an adapted reference could be used to better capture these difference. We will document this here in more detail in the future.


## Example:

Comparison of general COI reference for all vertebrate species and the reference used in Brasseuer et al. specific for Passeriformes (general vertebrate reference on top, Passeriformes (song birds) specific reference blow):

```{r, eval=TRUE}
-MFIXRWLFSTNHKDIGTLYLLFGAWAGMVGTALSLLIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFMVMPIMIGGFGNWLVPLMIGAPDMAFPRMN
MTFINRWLFSTNHKDIGTLYLIFGAWAGMVGTALSLLIRAELGQPGALLGDDQVYNVVVTAHAFVMIFFMVMPIMIGGFGNWLVPLMIGAPDMAFPRMN

NMSFWLLPPSFLLLLASSXVEAGAGTGWTVYPPLAGNLAHAGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAXSQYQTPLFVWSVLITAVLLLL
NMSFWLLPPSFLLLLASSTVEAGVGTGWTVYPPLAGNLAHAGASVDLAIFSLHLAGISSILGAINFITTAINMKPPALSQYQTPLFVWSVLITAVLLLL

SLPVLAAGITMLLTDRNLNTTFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWAHHMFTVG
SLPVLAAGITMLLTDRNLNTTFFDPAGGGDPVLYQHLFWFFGHPEVYILILPGFGIISHVVAYYAGKKEPFGYMGMVWAMLSIGFLGFIVWAHHMFTVG

MDVDTRAYFTSATMIIAIPTGVKVFSWLATLHGGNIKWXPAMLWALGFIFLFTVGGLTGIVLANSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMGGFVH
MDVDTRAYFTSATMIIAIPTGIKVFSWLATLHGGTIKWDPPMLWALGFIFLFTIGGLTGIVLANSSLDIALHDTYYVVAHFHYVLSMGAVFAILAGFTH

WFPLFTGYTLHXTWAKIHFXXMFVGVNLTFFPQHFLGLSGMPRRYSDYPDAYTXWNTVSSXGSFISLTAVILMXFIIWEAFAAKREVLXVELTXTNXEW
WFPLFTGYTLHSTWAKXHFGVMFVGVNLTFFPQHFLGLAGMPRRYSDYPDAYTLWNTISSVGSLISLTAVIMLVFIIWEAFASKRKALQPELTSTNVEW

LHGCPPPYHTFEEP---------
IHGCPPPFHTFEEPAFVQVQSQE
```

We have left out the stop codon at the end of the sequence in this example. For song birds there are two main variants of the COI gene, one with a stop codon at position 1551 and one with a stop codon at position 1554.

The references differ slightly. In particular, the Passeriformes has one additional amino acid at the beginning and 9 additional amino acids at the end.
These would not be extracted with the general reference. In the example analysis we show however that the middle part that includes the barcode region is extracted with both references and that the extracted barcode regions are identical. The general vertebrate reference can be a good starting point in many cases. For many other vertebrate groups the differences are even smaller than for song birds.

If the whole COI shall be extracted from start to end, specific references should be preferred.


## The folder Amino-Acid-references-for-taxonomic-groups

Reference sequences for different mitochondrial genes and different taxonomic groups shall be added over time to this repository, so that all users cannot only profit from the MitoGeneExtractor program, but also from a list of available references.


## The amino acid reference file given to MitoGeneExtractor

The reference file, which must be a file in the fasta format, has to contain at least one, but can contain any number of amino acid reference sequences. This can be different genes we want to extract in one program run and also different variants of the same genes. If this file contains multiple sequences, all input sequences will aligned with all references with MitoGeneExtractor.

The sequence names of the reference sequences will be used as part of the output filenames. Therefore they should not contain symbols not allowed in filenames. MGE might change sequence names if they contain illegal characters.

### Using multiple and different reference genes

The reference file can contain references for different genes. All genes will be reconstructed if possible. This speeds up the analyses for multiple genes.

Multiple reference genes and multiple variants of these can be combined in one reference file.

### Multiple variants of the same gene

If multiple significantly different variants of the same gene can occur in the input sequences, all variants can be added to the reference file. Input sequences that align to multiple references will be assigned to the reference for which the alignment has the highest alignment score. If input sequences have the same alignment score to multiple reference sequences, they are assigned to all references.

Multiple reference genes and multiple variants of these can be combined in one reference file.

### To which reference sequences are the input sequences assigned if they align with multiple reference sequences?

1) If a reads aligns to the same reference twice or more often, this is usually considered to be a problem. 
   - If the ```--includeDoubleHits``` command line option is specified, the read with the highest score is included and the one with the inferior score is ignored.
   - If the ```--includeDoubleHits``` command line option is not specified, non of the reads will be used.

2) If a read aligns to different reference sequences, e.g. to multiple sightly different references of the same gene reads are added to the references as follows:
  - If the ```--treat-references-as-individual``` command line option is specified, reads are assigned to all references.
  - If the ```--treat-references-as-individual``` command line option is not specified, reads are assigned to the reference for which the alignment score is highest. If the same best alignment score is obtained when aligning the read to multiple references, the read is added to all references.

2) allows multiple references to be specified and reads being assigned to all references in regions that are identical, but they are assigned to the best fitting reference in regions the references differ.

