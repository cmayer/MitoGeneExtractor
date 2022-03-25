# Amino acid reference sequences

Amino acid reference sequences are needed to align the query sequences, e.g. sequencing reads against. 
For conserved genes, such as mitochondrial genes, the same reference
can be used for large taxonomic groups. Even a single reference for all vertebrate species and a single refernce for all arthropod sequences exist for the COI gene, which can be used to obtain a large proportion of the full COI gene for these groups.
Only the slighly less conserved beginning and end of these genes require more specific references.

Example:

Comparison of general COI reference for all vertebrate species and the reference using in Brasseuer et al specific for Passeriformes (general vertebrate reference on top, Passeriformes blow):

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

The reference differ stlighly. In particular, the Passeriformes has one additional amino acid at the beginning and 9 additional amino acids at the end.
These would not be extracted with the general reference. In the example anaysis we show however that the barcode region extracted with both references
are identical in the given example. The general vertebrate reference can be a good starting point in many cases.

If the whole COI shall be extracted from start to end, specific references should be preferred.


References shall be added over time to this repository, so that all users can not only profit from the MitoGeneExtractor program, but also from a list of avaiable references.
