## Extraction from a reduced sequencing library file with fasta reads. Reference as in Brasseur et al. 2022:
../MitoGeneExtractor-v1.9.3 -d SRR12554985_trimmed_reduced.fas -p ../Amino-Acid-references-for-taxonomic-groups/COI-references/COI-fulllength-Passeriformes-protein-reference.fasta -V vulgar-SRR12554985_PasseriformesReference.txt -o SRR12554985_align_PasseriformesReference.fas -n 0 -c SRR12554985_cons_PasseriformesReference.fas -t 0.5 -r 1 -C 2

../MitoGeneExtractor-v1.9.3 -d SRR12554985_trimmed_reduced.fas -p ../Amino-Acid-references-for-taxonomic-groups/COI-references/COI-fulllength-general-vertebrata-protein-reference_from-consensus-50.fasta -V vulgar-SRR12554985_vertebrateReference.txt -o SRR12554985_align_vertebrateReference.fas -n 0 -c SRR12554985_cons_vertebrateReference.fas -t 0.5 -r 1 -C 2

