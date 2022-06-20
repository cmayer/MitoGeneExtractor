# Create protein references from NCBI protein sequences
In case you do not have any protein reference sequence available, you can download protein sequences from NCBI for your gene/taxon of interest and create a new reference. Here, we describe one possible way to achieve this, using COI and Songbirds (*Passeriformes*) as example. 
It might be worth to see whether there is already a suitable protein reference available under [this link](https://github.com/cmayer/MitoGeneExtractor/tree/main/Amino-Acid-references-for-taxonomic-groups).

If not, consider sharing your reference there, so that others benefit from that after you followed this step-by-step instruction:

1.) Visit the NCBI protein DB: https://www.ncbi.nlm.nih.gov/protein/. Click on 'advanced' to make an advanced search.

![1](https://user-images.githubusercontent.com/79691910/174590203-f34ce2f2-5622-4c77-97a2-fb65a2be7b13.jpg)

2.) Enter for 'Gene Name' your mitochondrial gene of interest (COX1), for 'Organism' the name of your taxonomic level of interest (*Passeriformes*) and for 'Sequence Length' the specific length range of your amino acid sequence. 
The complete CDS of COI has 516 amino acids, so we enter 510:520. Press 'search'.

![2](https://user-images.githubusercontent.com/79691910/174590215-02636974-a6a6-4eeb-807c-f0ae3d68cad1.jpg)

3.) Download all hits by pressing 'Send to:' and select 'File' in 'FASTA' format. Press 'Create File'. 

![3](https://user-images.githubusercontent.com/79691910/174590221-5519980a-a4fa-42af-8e77-51d8d2db0bff.jpg)

You will have a fasta file that should look somehow like this:


>UQK94632.1 cytochrome c oxidase subunit I (mitochondrion) [Corvus splendens]
MTFINRWLFSTNHKDIGTLYLIFGAWAGMVGTALSLLIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFM
VMPIMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSFLLLLASSTVEAGAGTGWTVYPPLAGNMAHA
GASVDLAIFSLHLAGISSILGAINFITTAINMKPPALSQYQTPLFVWSVLITAVLLLLSLPVLAAGITML
LTDRNLNTTFFDPAGGGDPVLYQHLFWFFGHPEVYILILPGFGIISHVVAYYAGKKEPFGYMGMVWAMLS
IGFLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAIPTGIKVFSWLATLHGGTIKWDPPMLWALGFIFLFT
IGGLTGIVLANSSLDIALHDTYYVVAHFHYVLSMGAVFAILAGFTHWFPLFTGYTLHSTWAKIHFGVMFV
GVNLTFFPQHFLGLAGMPRRYSDYPDAYTLWNTISSVGSLISLTAVIMLMFIIWEAFASKRKALQPELVN
TNVEWIHGCPPPFHTFEEPAFVQVQE

>UQK94619.1 cytochrome c oxidase subunit I (mitochondrion) [Corvus splendens]
MTFINRWLFSTNHKDIGTLYLIFGAWAGMVGTALSLLIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFM
VMPIMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSFLLLLASSTVEAGAGTGWTVYPPLAGNMAHA
GASVDLAIFSLHLAGISSILGAINFITTAINMKPPALSQYQTPLFVWSVLITAVLLLLSLPVLAAGITML
LTDRNLNTTFFDPAGGGDPVLYQHLFWFFGHPEVYILILPGFGIISHVVAYYAGKKEPFGYMGMVWAMLS
IGFLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAIPTGIKVFSWLATLHGGTIKWDPPMLWALGFIFLFT
IGGLTGIVLANSSLDIALHDTYYVVAHFHYVLSMGAVFAILAGFTHWFPLFTGYTLHSTWAKIHFGVMFV
GVNLTFFPQHFLGLAGMPRRYSDYPDAYTLWNTISSVGSLISLTAVIMLMFIIWEAFASKRKALQPELVN
TNVEWIHGCPPPFHTFEEPAFVQVQE


*Optional*: Use the supplied python script to retain only one sequence per genus. If only a few sequences are deposited in NCBI for your gene/taxon of interest, you can consider omitting this step. 
If you used the script, your fasta file should have the Genus as header, followed by the amino acid sequence.

5.) Use the cleaned or uncleaned fasta file to make an alignment. The alignment can be done e.g. with AliView (https://ormbunkar.se/aliview/). Using a graphical interface has the advantage that you can visually inspect your sequences and remove sequences that corrupt the alignment (e.g. due to Indels or bad quality sequences with missing data).
In AliView, select 'File' -> 'Open file', navigate to the fasta file location on your system and open the file.
The protein sequence is pretty conserved for most of the taxa, although some genera are considerably different. This is fine, but we should filter for bad quality sequences with lots of missing data (e.g. more than 3 'X' in a row). 

![6](https://user-images.githubusercontent.com/79691910/174590297-48de81de-6a90-4875-8154-1899a5148d6d.jpg)

To delete a sequence, select the taxon in the left panel, right click on it and press 'Delete selected sequence(s)'. Further, it is advised to exclude sequences that are longer than one would expect, because the additional nucleotide positions will be included in the alignment.
Additionally, you can exclude sequences that are much shorter than you would expect, although these should not compromise the alignment if they are rare.

![7](https://user-images.githubusercontent.com/79691910/174590317-d6677ea1-cca7-467e-9f6f-4bdd10911db6.jpg)

6.) Once you are done, select 'Align' -> 'Realign everything'. You can use the default alignment algorithm of AliView (MUSCLE).

![8](https://user-images.githubusercontent.com/79691910/174590472-be4475d1-20ec-4d63-a2cc-2fa19e071e36.jpg)

7.) After the first round of alignment, you will probably have some sequences, that corrupt the alignment because they have some uncommon Indels (Columns with almost only gaps). Remove the sequences that introduce these gaps, then press again 'Realign everything'. 
If now other samples with Indels pop up, remove them again and realign the data. Repeat this, until no more columns exist where most samples have missing data (gaps).

![9](https://user-images.githubusercontent.com/79691910/174590496-56068dd3-c936-4d8c-815d-d78589bf49c5.jpg)

8.) You are done when your alignment looks like this: No major gap columns within the CDS. Gaps at the start/end are not very problematic. However, they can yield in a longer reference than the actual legnth of the CDS if they are not removed.
In the example below, the COI gene actually starts at position 4. In some cases, the sequences in database entries are a bit too long and here they have three additional amino acids at the start.
You can remove such sequences or leave them in. If they are included, keep in mind that MitoGeneExtractor will reconstruct slightly longer sequences (3 AA = 9 nucleotides more in the beginning).

![11](https://user-images.githubusercontent.com/79691910/174590684-227c483e-9608-4008-bfde-2de5007b2fa2.jpg)

9.) Save the final alignment by pressing 'File' -> 'Save as Fasta'. Create one final consensus sequence as reference for MitoGeneExtractor. You could use for that e.g. EMBOSS (https://www.ebi.ac.uk/Tools/msa/emboss_cons/).

In case you have any questions/issues, post them in Github.
