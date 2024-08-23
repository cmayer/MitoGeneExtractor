###Grabs closest cox1 protein reference sequence for each sample using BOLD-downloaded taxonomic rankings, 
###outputs the cox1.fasta seqs to a user specified dir, and outputs Process ID, matched term, accession number, 
###ref seq name, absolute path to ref seq, and validates returned protein reference taxonomy to BOLD-downloaded taxonomy

import csv
import sys
import time
import os
from Bio import Entrez, SeqIO
from typing import Optional, List




#Set your email and API key
Entrez.email = "############"  # Add your email for Entrez here
Entrez.api_key = ##############  # Add your NCBI API key here




def fetch_refseq_protein_sequences_by_taxonomy(taxonomy, gene_name, retmax=1):
    best_rank = None
    best_records = []
    
#Iterate over ranks from Order to Species
    for rank in taxonomy:
        if rank:
            search_term = f"{gene_name}[Gene] AND {rank}[Organism] AND refseq[filter]"
            try:
                search_handle = Entrez.esearch(db="protein", term=search_term, retmax=retmax)
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                if search_results["IdList"]:
                    ids = search_results["IdList"]
                    print(f"Fetching sequences with IDs: {ids} for rank: {rank}")  # Print rank here
                    fetch_handle = Entrez.efetch(db="protein", id=ids, rettype="gb", retmode="text")
                    records = list(SeqIO.parse(fetch_handle, "genbank"))
                    fetch_handle.close()
                    
#Update best records and rank if seq found
                    best_records = records
                    best_rank = rank
                    break  # Stop searching once a match is found
                else:
                    print(f"No sequences found for {rank}. Trying next level...")
            except Exception as e:
                print(f"Error fetching data for {rank}: {e}")

            time.sleep(1)
    
    if not best_records:
        print(f"No sequences found for taxonomy {taxonomy} and gene {gene_name}.")
    
    return best_records, best_rank




#Fetch taxonomic ID for returned protein reference
def fetch_tax_id_by_name(organism_name):
    search_handle = Entrez.esearch(db="taxonomy", term=organism_name, retmax=1)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    if search_results["IdList"]:
        return search_results["IdList"][0]
    return None




#Compare taxonomic ranks between BOLD taxonomy and returned protein reference
def validate_taxonomy(input_taxonomy, fetched_taxonomy):
    ranks_to_check = ["Order", "Class", "Phylum"]
    for rank in ranks_to_check:
        input_rank_value = input_taxonomy.get(rank, None)
        fetched_rank_value = next((lineage['ScientificName'] for lineage in fetched_taxonomy['LineageEx'] if lineage['Rank'] == rank.lower()), None)
        
        if input_rank_value and fetched_rank_value:
            if input_rank_value.lower() != fetched_rank_value.lower():
                return False
    return True




def usage():
    print("""
    Usage: python 1_gene_fetch.py path/to/input_taxonomy_csv [gene_name] path/to/protein_references/output_dir

    <input_taxonomy_csv>: Path to input CSV file containing taxonomy information from BOLD download.
    <gene_name>: Name of gene to search for in NCBI RefSeq database (e.g. cox1).
    <output_directory>: Path to directory to save output files (will save .fasta(s) and output.csv in this directory. must specify 'protein_references' dir)).
    """)




if __name__ == "__main__":
    if len(sys.argv) != 4:
        usage()
        sys.exit(1)

    input_file = sys.argv[1]
    gene_name = sys.argv[2]
    user_output_directory = sys.argv[3]

    if not os.path.exists(input_file):
        print(f"Error: The file '{input_file}' does not exist.")
        usage()
        sys.exit(1)

#Check user-specified output directory exists
    if not os.path.exists(user_output_directory):
        os.makedirs(user_output_directory)

    summary_output = []

#Search 'down' taxonomic ranks, starting from Order in BOLD taxonomy
    with open(input_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            taxonomy = {
                'Order': row['Order'],
                'Family': row['Family'],
                'Genus': row['Genus'],
                'Species': row['Species']
            }
            process_id = row['Process ID']
            records, matched_rank = fetch_refseq_protein_sequences_by_taxonomy(taxonomy.values(), gene_name, retmax=1)


#Fetch taxonomy from NCBI for matched_rank using organism name extracted from returned protein sequence to get corresponding tax_id
#Fetch full taxonomy using tax_id and compare against input taxonomy
            if records:
                organism_name = records[0].annotations['organism']  
                tax_id = fetch_tax_id_by_name(organism_name)
                
                if tax_id:
                    fetched_taxonomy = fetch_taxonomy_by_id(tax_id) 
                    matched_validated = validate_taxonomy(taxonomy, fetched_taxonomy) 
                else:
                    matched_validated = False 
            else:
                matched_validated = False 

#Write seqs to fasta with process ID as filename and seq header
            if records:
                for record in records:
                    accession_number = record.id  
                    record.id = process_id  
                    record.description = ""  
                    fasta_filename = f"{process_id}.fasta"
                    fasta_filepath = os.path.abspath(os.path.join(user_output_directory, fasta_filename))
                    with open(fasta_filepath, "w") as output_handle:
                        SeqIO.write([record], output_handle, "fasta")                

#Add fields to summary output.csv
                    summary_output.append({
                        'process_id': process_id,
                        'matched_term': matched_rank,
                        'accession_number': accession_number,  
                        'reference_name': process_id,  
                        'reference_path': fasta_filepath,
                        'matched_validated': 'true' if matched_validated else 'false'
                    })

#Print accession number for fetched seq and its taxonomic rank
                    print(f"Fetched sequence ID: {accession_number} for rank: {matched_rank}")

            time.sleep(1)

#Write output.csv to user-specified directory
    summary_output_file = os.path.join(user_output_directory, f"{os.path.splitext(os.path.basename(input_file))[0]}_gene_fetch_sum_out.csv")
    with open(summary_output_file, "w", newline='') as csvfile:
        fieldnames = ['process_id', 'matched_term', 'accession_number', 'reference_name', 'reference_path', 'matched_validated']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_output)

    print(f"Summary output saved to {summary_output_file}.")
