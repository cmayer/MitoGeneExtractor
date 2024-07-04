# Users have asked many question. The answers are likely of interest for others as well.

## List of questions:

### MitoGeneExtractor reports: "ERROR: Running exonerate failed. The generated vulgar file is incomplete and should be removed manually. Exiting."

MitoGeneExtractor calls the exonerate program to align the input sequences against the specified references.
This error occurs if the exonerate program crashes. So far all cases in which exonerate crashed, this was most likely due to insufficient RAM.
Rerunning exonerate on a machine with more memory solved the problem. Currently we are working a version of MitoGeneExtractor that chops the input file into smaller files and runs them separately. 

### Sequences found with MitoGeneExtractor have a minimum length of about 60.

This can be changes by specifying the -s x option of MGE where x is the score passed to exonerate.
The logic behind this is the following. Exonerate aligns the nucleotide input sequences to the amino acid reference.
This is done by translating the nucleotide sequences to all 3 reading frames and aligning them to the amino acid reference.
Each aligned amino acid gets a similarity score which by default is the blosum62 score. If the total score is above 100, which is the 
default score set by exonerate, the alignment is included in the vulgar file returned by exonerate.
Since the mean blosum62 score of matching amino acids is close to 5, about 60 nucleotides or 20 amoni acids are needed to exceed the score
threshold. The exact theshold depends on whether all amino acids match exactly, which they might not.

If you want that even shorter sequences are reported by MitoGeneExtractor, you need to use the -s parameter to specify a threshold score lower than
100. A rough estimate of the threshold score required to report input sequences of length Lmin is: Lmin*5/3.

Lowering the score threshold has the following downside: More bad alignments are reported. In order to avoid this effect, the -r <float>, --relative_score_threshold <float> parameter has been introduced. A typical perfectly matching alignment has a relative score of 5 per codon or 3 nucleotides, which corresponds to a relative score of 5/3=1.67. The default threshold of MitoGeneExtractor for the relative alignment score is 1.0, which allows a small number of mismatching amino acids in the exonerate alingment. With the -r <float>, --relative_score_threshold <float> parameter the allowed number of non indenctical amino acids in be controlled.

