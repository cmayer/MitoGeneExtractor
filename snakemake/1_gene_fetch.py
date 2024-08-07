###Grabs closest cox1 protein reference sequence for each sample using BOLD-downloaded taxonomic rankings, 
###outputs the cox1.fasta seqs to a user specified dir, and outputs Process ID, matched term, accession number, 
###ref seq name, and absolute path to ref seq to .csv

import csv
import sys
import time
import os
from Bio import Entrez, SeqIO



#Set your email and API key
Entrez.email = "b.price@nhm.ac.uk"  #Add your email for Entrez here
Entrez.api_key = "82df5a6f5cf735302d3cf1fcf48b206cfe09"  #Add your NCBI API key here



def fetch_refseq_protein_sequences_by_taxonomy(taxonomy, gene_name, retmax=1):
    for rank in taxonomy:
        if rank:
            search_term = f"{gene_name}[Gene] AND {rank}[Organism] AND refseq[filter]"
            try:
                search_handle = Entrez.esearch(db="protein", term=search_term, retmax=retmax)
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                if search_results["IdList"]:
                    ids = search_results["IdList"]
                    print(f"Fetching sequences with IDs: {ids}")  
                    fetch_handle = Entrez.efetch(db="protein", id=ids, rettype="gb", retmode="text")
                    records = list(SeqIO.parse(fetch_handle, "genbank"))
                    fetch_handle.close()
                    
                    return records, rank
                else:
                    print(f"No sequences found for {rank}. Trying next level...")
            except Exception as e:
                print(f"Error fetching data for {rank}: {e}")
            
#Delay between searches
            time.sleep(1)
    
    print(f"No sequences found for taxonomy {taxonomy} and gene {gene_name}.")
    return [], None




def usage():
    """
    Print usage information for the script.
    """
    print("""
    Usage: python 1_mge_fetch.py <input_csv_file> <gene_name> <output_directory>

    <input_csv_file>: Path to input CSV file containing taxonomy information from BOLD download.
    <gene_name>: Name of gene to search for in NCBI RefSeq database (e.g. cox1).
    <output_directory>: Directory to save output 'protein_references' folder (will save .fasta(s) and output.csv in protein_references dir within specified output directory).
    """)



if __name__ == "__main__":
    if len(sys.argv) != 4:
        usage()
        sys.exit(1)

    input_file = sys.argv[1]
    gene_name = sys.argv[2]
    output_directory = sys.argv[3]

    if not os.path.exists(input_file):
        print(f"Error: The file '{input_file}' does not exist.")
        usage()
        sys.exit(1)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    protein_references_dir = os.path.join(output_directory, "protein_references")
    summary_output = []


#Create directory for output files
    os.makedirs(protein_references_dir, exist_ok=True)

    with open(input_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            taxonomy = [row['Species'], row['Genus'], row['Family'], row['Order'], row['Class'], row['Phylum']]
            process_id = row['Process ID']
            records, matched_rank = fetch_refseq_protein_sequences_by_taxonomy(taxonomy, gene_name, retmax=1)


#Write seqs to FASTA with process ID as filename and seq header
            if records:
                for record in records:
                    accession_number = record.id  # Save the original accession number
                    record.id = process_id  # Use process ID as the sequence ID
                    record.description = ""  # Clear the description
                    fasta_filename = f"{process_id}.fasta"
                    fasta_filepath = os.path.abspath(os.path.join(protein_references_dir, fasta_filename))
                    with open(fasta_filepath, "w") as output_handle:
                        SeqIO.write([record], output_handle, "fasta")                

#Add to summary output.csv
                    summary_output.append({
                        'process_id': process_id,
                        'matched_term': matched_rank,
                        'accession_number': accession_number,  
                        'reference_name': process_id,  
                        'reference_path': fasta_filepath
                    })


#Delay between searches
            time.sleep(1)


#Write output.csv to protein_references directory
    summary_output_file = os.path.join(protein_references_dir, f"{os.path.splitext(os.path.basename(input_file))[0]}_mge_fetch_sum_out.csv")
    with open(summary_output_file, "w", newline='') as csvfile:
        fieldnames = ['process_id', 'matched_term', 'accession_number', 'reference_name', 'reference_path']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        writer.writerows(summary_output)


    print(f"Summary output saved to {summary_output_file}.")
