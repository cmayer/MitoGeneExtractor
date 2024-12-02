# ChangeLog

## Version 1.9.3
First Version on Github.

## Version 1.9.5, February 5th, 2023

* Added more mitochondrial reference sequences.
* Multiple amino acid references can now be analyses in one run. Input sequences will be assigned to the reference for which it has the highest score.
If scores are equal, input sequences will be added to multiple references.
* Now multiple fasta file and multiple fastq will can be specified. They will be merged before being analysed. The old -d and new -q option are now MultiArg options.
* A new opion was added: --temporaryDirectory <string>, --keep-concat-input-file


## Version 1.9.6beta, December 1st, 2024

* refactored the code for future changes.
* Fixed bugs that could result in crashes with segmentation fault. I have so far not found changes in the output, but I cannot exclude this at the moment. 
* currently, the changes are only in the new branch https://github.com/cmayer/MitoGeneExtractor/tree/refactor-bugfix-v1.9.6
