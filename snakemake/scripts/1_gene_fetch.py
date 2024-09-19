###Grabs closest cox1 protein reference sequence for each sample using BOLD-downloaded taxonomic rankings, 
###outputs the cox1.fasta seqs to a user specified dir, and outputs Process ID, matched term, accession number, 
###ref seq name, absolute path to ref seq, and validates returned protein reference taxonomy to BOLD-downloaded taxonomy


import csv
import sys
import time
import os
from Bio import Entrez, SeqIO
from typing import Optional, List




#USAGE
def usage():
    print("""
    Usage: python 1_gene_fetch.py path/to/input_taxonomy_tsv [gene_name] path/to/protein_references/output_dir

    <input_taxonomy_tsv>: Path to input tsv file containing taxonomy information from BOLD download.
    <gene_name>: Name of gene to search for in NCBI RefSeq database (e.g. cox1).
    <output_directory>: Path to directory to save output files (will save .fasta(s) and output.csv in this directory. must specify 'protein_references' dir)).
    """)




#Set email and API key
Entrez.email = "####@##.ac.uk"  # Add your email for Entrez here
Entrez.api_key = "####"  # Add your NCBI API key here



#Fetch protein reference using input taxonomic ranks
def fetch_protein_seq_by_taxonomy(taxonomy, gene_name, retmax=1):
    best_rank = None
    best_records = []
    found_sequence = False
    rank_priority = ["Order", "Family", "Genus", "Species"]

    #Iterate over ranks from Order to Species
    for rank, rank_name in zip(taxonomy, rank_priority):
        if rank:
            search_term = f"{gene_name}[Gene] AND {rank}[Organism] AND refseq[filter]"
            try:
                search_handle = Entrez.esearch(db="protein", term=search_term, retmax=retmax)
                search_results = Entrez.read(search_handle)
                search_handle.close()

                if search_results["IdList"]:
                    ids = search_results["IdList"]
                    print(f"Fetching {gene_name} sequence: {ids} for {rank_name}: {rank}")
                    fetch_handle = Entrez.efetch(db="protein", id=ids, rettype="gb", retmode="text")
                    records = list(SeqIO.parse(fetch_handle, "genbank"))
                    fetch_handle.close()

                    #Update best records and rank if found a better match (lower rank)
                    best_records = records
                    best_rank = f"{rank_name}: {rank}"
                    found_sequence = True  # Set the flag indicating a sequence was found
                else:
                    print(f"No sequences found for {rank_name}: {rank}. Trying next level...")
                    if found_sequence:
                        break  # Exit the loop early if a sequence was already found

            except Exception as e:
                print(f"Error fetching data for {rank_name}: {rank}: {e}")

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



#Fetch full taxonomy of returned protein seq using tax_id
def fetch_taxonomy_by_id(tax_id):
    try:
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        if records:
            return records[0]  # Return the first (and usually only) record
        else:
            return None
    except Exception as e:
        print(f"Error fetching taxonomy for tax_id {tax_id}: {e}")
        return None


#Compare 'Order' taxonomic rank between BOLD taxonomy and returned protein reference
def validate_order(input_order, fetched_taxonomy):
    if not fetched_taxonomy:
        print("Fetched taxonomy is None, cannot validate.")
        return False

    fetched_order = next((lineage['ScientificName'] for lineage in fetched_taxonomy.get('LineageEx', []) if lineage['Rank'] == 'order'), None)

    if input_order and fetched_order:
        if input_order.lower() != fetched_order.lower():
            return False
    return True



#Specify format of output taxonomy validation 
def format_taxonomy(lineage):
    # Filter out clades and start from 'Kingdom'
    rank_order = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    taxonomy_str = "; ".join(f"{lineage_item['Rank'].capitalize()}: {lineage_item['ScientificName']}" 
                              for lineage_item in lineage 
                              if lineage_item['Rank'].lower() in rank_order and lineage_item['Rank'].lower() != 'clade')
    return taxonomy_str



#MAIN
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
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            taxonomy = {
                'Order': row['Order'],
                'Family': row['Family'],
                'Genus': row['Genus'],
                'Species': row['Species']
            }
            process_id = row['Process ID']
            records, matched_rank = fetch_protein_seq_by_taxonomy(taxonomy.values(), gene_name, retmax=1)

             #Fetch NCBI taxonomy for matched_rank using organism name extracted from returned protein sequence to get corresponding tax_id
            if records:
                organism_name = records[0].annotations['organism']
                tax_id = fetch_tax_id_by_name(organism_name)

                if tax_id:
                    fetched_taxonomy = fetch_taxonomy_by_id(tax_id)
                    matched_validated = validate_order(taxonomy['Order'], fetched_taxonomy)
                else:
                    fetched_taxonomy = None
                    matched_validated = False
            else:
                fetched_taxonomy = None
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

                    #Convert fetched taxonomy to string
                    if fetched_taxonomy:
                        lineage = fetched_taxonomy.get('LineageEx', [])
                        taxonomy_str = format_taxonomy(lineage)
                    else:
                        taxonomy_str = "Not available"

                    #Add fields to summary output.csv
                    summary_output.append({
                        'process_id': process_id,
                        'matched_term': matched_rank,
                        'accession_number': accession_number,
                        'reference_name': process_id,
                        'reference_path': fasta_filepath,
                        'matched_validated': 'true' if matched_validated else 'false',
                        'ncbi_taxonomy': taxonomy_str  # Add filtered taxonomy to output
                    })

                    # Print accession number for fetched seq and its taxonomic rank
                    print(f"Fetched sequence accession: {accession_number} for {matched_rank}")

            time.sleep(1)

    #Write output.csv to user-specified directory
    summary_output_file = os.path.join(user_output_directory, f"{os.path.splitext(os.path.basename(input_file))[0]}_gene_fetch_sum_out.csv")
    with open(summary_output_file, "w", newline='') as csvfile:
        fieldnames = ['process_id', 'matched_term', 'accession_number', 'reference_name', 'reference_path', 'matched_validated', 'ncbi_taxonomy']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_output)

    print(f"Summary output saved to {summary_output_file}.")
