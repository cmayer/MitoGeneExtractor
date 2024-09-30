###Grabs closest cox1 protein reference sequence for each sample using BOLD-downloaded taxonomic rankings, 
###outputs the cox1.fasta seqs to a user specified dir and outputs a protein_reference.csv


import csv
import sys
import time
import os
from Bio import Entrez, SeqIO
from functools import wraps



# Set your email and API key
Entrez.email = ########"  # Add your email for Entrez here
Entrez.api_key = "#####"  # Add your NCBI API key here




def usage():
    """
    Prints the usage instructions for the script.
    """
    print("""
    Usage: python 1_gene_fetch.py <input_taxonomy_file> <gene_name> <output_directory> <samples.csv>

    <input_taxonomy_file>: Path to input taxonomy file containing taxonomy information from BOLD download.
                           Supported formats: CSV or TSV (based on file extension).
    <gene_name>: Name of gene to search for in NCBI RefSeq database (e.g., cox1).
    <output_directory>: Path to directory to save output files (will save .fasta files and summary CSV in this directory).
                        The directory will be created if it does not exist.
    <samples.csv>: Path to input CSV file containing Process IDs to match (with column 'ID').

    Example:
    python 1_gene_fetch.py taxonomy.csv cox1 ./protein_references samples.csv
    """)




def fetch_protein_seq_by_taxonomy(taxonomy, gene_name, retmax=1):
    """
    Fetches protein sequences from NCBI RefSeq based on taxonomy and gene name.

    Parameters:
    - taxonomy (iterable): Taxonomic ranks (order, family, genus, species).
    - gene_name (str): Name of the gene to search for.
    - retmax (int): Maximum number of records to retrieve.

    Returns:
    - tuple: (list of SeqRecord objects, matched_rank as str)
    """
    best_rank = None
    best_records = []
    found_sequence = False
    rank_priority = ["order", "family", "genus", "species"]

    # Iterate over ranks from order to species
    for rank, rank_name in zip(taxonomy, rank_priority):
        if rank:
            search_term = f"{gene_name}[Gene] AND {rank}[Organism] AND refseq[filter]"
            try:
                search_handle = Entrez.esearch(db="protein", term=search_term, retmax=retmax)
                search_results = Entrez.read(search_handle)
                search_handle.close()

                if search_results["IdList"]:
                    ids = search_results["IdList"]
                    print(f"Fetching {gene_name} sequence: {ids} for {rank_name.capitalize()}: {rank}")
                    fetch_handle = Entrez.efetch(db="protein", id=ids, rettype="gb", retmode="text")
                    records = list(SeqIO.parse(fetch_handle, "genbank"))
                    fetch_handle.close()

                    # Update best records and rank if found a better match (lower rank)
                    best_records = records
                    best_rank = f"{rank_name.capitalize()}: {rank}"
                    found_sequence = True  #Set  flag indicating sequence found
                else:
                    print(f"No sequences found for {rank_name.capitalize()}: {rank}. Trying next level...")
                    if found_sequence:
                        break  #Exit loop if sequence already found

            except Exception as e:
                print(f"Error fetching data for {rank_name.capitalize()}: {rank}: {e}")

            time.sleep(1)

    if not best_records:
        print(f"No sequences found for taxonomy {taxonomy} and gene {gene_name}.")

    return best_records, best_rank




def fetch_tax_id_by_name(organism_name):
    """
    Fetches the taxonomic ID for a given organism name from NCBI.

    Parameters:
    - organism_name (str): Scientific name of the organism.

    Returns:
    - str or None: Taxonomic ID if found, else None.
    """
    search_handle = Entrez.esearch(db="taxonomy", term=organism_name, retmax=1)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    if search_results["IdList"]:
        return search_results["IdList"][0]
    return None




def fetch_taxonomy_by_id(tax_id):
    """
    Retrieves full taxonomy information for a given taxonomic ID.

    Parameters:
    - tax_id (str): Taxonomic ID.

    Returns:
    - dict or None: Taxonomy information if found, else None.
    """
    try:
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        if records:
            return records[0]
        else:
            return None
    except Exception as e:
        print(f"Error fetching taxonomy for tax_id {tax_id}: {e}")
        return None




def validate_order(input_order, fetched_taxonomy):
    """
    Validates if the 'Order' rank in fetched taxonomy matches the input taxonomy.

    Parameters:
    - input_order (str): The 'order' rank from input taxonomy.
    - fetched_taxonomy (dict): The fetched taxonomy information from NCBI.

    Returns:
    - bool: True if matched or input_order is not provided, False otherwise.
    """
    if not fetched_taxonomy:
        print("Fetched taxonomy is None, cannot validate.")
        return False

    fetched_order = next(
        (lineage['ScientificName'] for lineage in fetched_taxonomy.get('LineageEx', []) if lineage['Rank'].lower() == 'order'),
        None
    )

    if input_order and fetched_order:
        if input_order.lower() != fetched_order.lower():
            return False
    return True




def format_taxonomy(lineage):
    """
    Formats the taxonomy lineage into a readable string.

    Parameters:
    - lineage (list): List of lineage dictionaries from NCBI taxonomy.

    Returns:
    - str: Formatted taxonomy string.
    """
    # Filter out clades and start from 'Kingdom'
    rank_order = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    taxonomy_str = "; ".join(
        f"{lineage_item['Rank'].capitalize()}: {lineage_item['ScientificName']}"
        for lineage_item in lineage
        if lineage_item['Rank'].lower() in rank_order and lineage_item['Rank'].lower() != 'clade'
    )
    return taxonomy_str




def determine_delimiter(file_path):
    """
    Determines the delimiter of the input taxonomy file based on its extension.

    Parameters:
    - file_path (str): Path to the taxonomy file.

    Returns:
    - str: Delimiter character.
    """
    _, ext = os.path.splitext(file_path)
    if ext.lower() == '.csv':
        # Further check for common delimiters
        with open(file_path, 'r') as f:
            first_line = f.readline()
            if ';' in first_line:
                return ';'
            else:
                return ','
    elif ext.lower() == '.tsv':
        return '\t'
    else:
        print("Unsupported file extension. Please provide a .csv or .tsv file for taxonomy.")
        sys.exit(1)




if __name__ == "__main__":
    if len(sys.argv) != 5:
        usage()
        sys.exit(1)

    input_file = sys.argv[1]
    gene_name = sys.argv[2]
    user_output_directory = sys.argv[3]
    samples_file = sys.argv[4]


    #Check input taxonomy exists
    if not os.path.exists(input_file):
        print(f"Error: The file '{input_file}' does not exist.")
        usage()
        sys.exit(1)


    #Check samples.csv exists
    if not os.path.exists(samples_file):
        print(f"Error: The file '{samples_file}' does not exist.")
        usage()
        sys.exit(1)


    #Determine delimiter based on file extension
    delimiter = determine_delimiter(input_file)
    print(f"Using delimiter '{delimiter}' for taxonomy file.")


    #Check user-specified output directory exists, create if not
    if not os.path.exists(user_output_directory):
        os.makedirs(user_output_directory)
        print(f"Created output directory: {user_output_directory}")


    #Load Process IDs from samples.csv
    valid_process_ids = set()
    try:
        with open(samples_file, newline='') as samples_csv:
            reader = csv.DictReader(samples_csv)
            for row in reader:
                valid_process_ids.add(row['ID'].strip())
        print(f"Loaded {len(valid_process_ids)} Process IDs from '{samples_file}'.")
    except Exception as e:
        print(f"Error reading samples file '{samples_file}': {e}")
        sys.exit(1)

    summary_output = []


    #Process taxonomy file
    try:
        with open(input_file, newline='') as taxonomy_file:
            reader = csv.DictReader(taxonomy_file, delimiter=delimiter)

            fieldnames_lower = [field.lower() for field in reader.fieldnames]

            for row in reader:

                row_lower = {k.lower(): v for k, v in row.items()}

                process_id = row_lower.get('process id', '').strip()

                if not process_id:
                    print("Warning: Missing 'Process ID' in a row. Skipping.")
                    continue

                print(f"\nProcessing Process ID: {process_id}")



                #Skip rows where Process ID doesn't match the samples.csv Process IDs
                if process_id not in valid_process_ids:
                    print(f"Skipping Process ID: {process_id} (not in samples.csv)")
                    continue


                #Extract taxonomy ranks with case-insensitive keys
                taxonomy = {
                    'order': row_lower.get('order', '').strip(),
                    'family': row_lower.get('family', '').strip(),
                    'genus': row_lower.get('genus', '').strip(),
                    'species': row_lower.get('species', '').strip()
                }

                print(f"Taxonomy: order='{taxonomy['order']}', family='{taxonomy['family']}', genus='{taxonomy['genus']}', species='{taxonomy['species']}'")

                records, matched_rank = fetch_protein_seq_by_taxonomy(taxonomy.values(), gene_name, retmax=1)


                #Fetch NCBI taxonomy for matched_rank using organism name extracted from the returned protein sequence to get corresponding tax_id
                if records:
                    organism_name = records[0].annotations.get('organism', '')
                    tax_id = fetch_tax_id_by_name(organism_name)

                    if tax_id:
                        fetched_taxonomy = fetch_taxonomy_by_id(tax_id)
                        matched_validated = validate_order(taxonomy['order'], fetched_taxonomy)
                    else:
                        fetched_taxonomy = None
                        matched_validated = False
                else:
                    fetched_taxonomy = None
                    matched_validated = False


                #Write sequences to fasta with process ID as filename and seq header
                if records:
                    for record in records:
                        accession_number = record.id
                        record.id = process_id
                        record.description = ""
                        fasta_filename = f"{process_id}.fasta"
                        fasta_filepath = os.path.abspath(os.path.join(user_output_directory, fasta_filename))
                        try:
                            with open(fasta_filepath, "w") as output_handle:
                                SeqIO.write([record], output_handle, "fasta")
                            print(f"Written FASTA file: {fasta_filepath}")
                        except Exception as e:
                            print(f"Error writing FASTA file '{fasta_filepath}': {e}")
                            continue

                        if fetched_taxonomy:
                            lineage = fetched_taxonomy.get('LineageEx', [])
                            taxonomy_str = format_taxonomy(lineage)
                        else:
                            taxonomy_str = "Not available"


                        #Add fields to summary output.csv
                        summary_output.append({
                            'process_id': process_id,
                            'matched_term': matched_rank if matched_rank else "None",
                            'accession_number': accession_number,
                            'reference_name': process_id,
                            'reference_path': fasta_filepath,
                            'matched_validated': 'true' if matched_validated else 'false',
                            'ncbi_taxonomy': taxonomy_str  # Add filtered taxonomy to output
                        })


                        #Print accession number for fetched seq and its taxonomic rank
                        print(f"Fetched sequence accession: {accession_number} for {matched_rank}")

                time.sleep(1)

    except Exception as e:
        print(f"Error processing taxonomy file '{input_file}': {e}")
        sys.exit(1)


    #Write summary_output to CSV
    summary_output_file = os.path.join(
        user_output_directory,
        f"{os.path.splitext(os.path.basename(input_file))[0]}_gene_fetch_summary.csv"
    )
    try:
        with open(summary_output_file, "w", newline='') as summary_csv:
            fieldnames = [
                'process_id',
                'matched_term',
                'accession_number',
                'reference_name',
                'reference_path',
                'matched_validated',
                'ncbi_taxonomy'
            ]
            writer = csv.DictWriter(summary_csv, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(summary_output)
        print(f"\nSummary output saved to {summary_output_file}.")
    except Exception as e:
        print(f"Error writing summary CSV '{summary_output_file}': {e}")
        sys.exit(1)

        print(f"Error writing summary CSV '{summary_output_file}': {e}")
        sys.exit(1)
