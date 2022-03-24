## Extraction from a reduced sequencing library file with fasta reads. Reference as in Brasseur et al. 2022:
/Daten/Programmieren/C++/Programme/align-trim-DNA-against-Protein/MitoGeneExtractor-v1.9.1-dist/MitoGeneExtractor-v1.9.1 -d SRR12554985_trimmed_reduced.fas -p ../Amino-Acide-references-for-tanomic-groups/xxx-need-to-be-added.fasta -V vulgar-SRR12554985_PasseriformesReference.txt -o SRR12554985_align_PasseriformesReference.fas -n 0 -c SRR12554985_cons_PasseriformesReference.fas -t 0.5 -r 1 -C 2

/Daten/Programmieren/C++/Programme/align-trim-DNA-against-Protein/MitoGeneExtractor-v1.9.1-dist/MitoGeneExtractor-v1.9.1 -d SRR12554985_trimmed_reduced.fas -p ../Amino-Acide-references-for-tanomic-groups/COI-vertebrata-protein-consensus-50_whole.fasta -V vulgar-SRR12554985_vertebrateReference.txt -o SRR12554985_align_vertebrateReference.fas -n 0 -c SRR12554985_cons_vertebrateReference.fas -t 0.5 -r 1 -C 2

