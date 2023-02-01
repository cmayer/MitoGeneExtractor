# Amino acid reference sequences

Amino acid reference sequences are needed to align the query sequences, e.g. sequencing reads against. 
For conserved genes, such as mitochondrial genes, the same reference
can be used for large taxonomic groups. Even a single reference for all vertebrate species and a single reference for all arthropod sequences exist for the COI gene, which can be used to obtain a large proportion of the full COI gene for these groups.
Only the slightly less conserved beginning and end of these genes require more specific references.

Example:

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

The references differ slightly. In particular, the Passeriformes has one additional amino acid at the beginning and 9 additional amino acids at the end.
These would not be extracted with the general reference. In the example analysis we show however that the barcode region extracted with both references
are identical in the given example. The general vertebrate reference can be a good starting point in many cases. For many other vertebrate groups the differences
are even smaller than for song birds.

If the whole COI shall be extracted from start to end, specific references should be preferred.


Reference sequences for different mitochondrial genes shall be added over time to this repository, so that all users can not only profit from the MitoGeneExtractor program, but also from a list of available references.

## The reference file

The reference file, which must be a file in the fasta format, has to contain one or multiple amino acid reference sequences.

In an MGE run, the input nucleotide sequences will be aligned (with Exonerate) against all reference sequences in this file and reads that align well are retained.

The sequence names of the reference sequences will be used as part of the output filenames. Therefore they should not contain symbols not allowed in filenames. MGE might change sequence names if they contain illegal characters.

### Using multiple and different reference genes

The reference file can contain references for different genes. All genes will be reconstructed if possible. This speeds up the analyses for multiple genes.

Multiple reference genes and multiple variants of these can be combined in one reference file.

### Multiple variants of the same gene

If multiple significantly different variants of the same gene can occur in the input sequences, all variants can be added to the reference file. Input sequences that align to multiple references will be assigned to the reference for which the alignment has the highest alignment score. If input sequences have the same alignment score to multiple reference sequences, they are assigned to all references.

Multiple reference genes and multiple variants of these can be combined in one reference file.

