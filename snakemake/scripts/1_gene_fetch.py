import csv
import sys
import time
import os
import logging
from Bio import Entrez, SeqIO
from functools import wraps
from ratelimit import limits, sleep_and_retry
from requests.exceptions import RequestException




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
    python 1_gene_fetch.py path/to/taxonomy.csv cox1 path/to/protein_references path/to/samples.csv
    """)




#Set email and API key
Entrez.email = "#####@#####"  
Entrez.api_key = "#########"  


#Consta

nts for rate limiting (Entrez API limits with key = 10/sec)
MAX_CALLS_PER_SECOND = 10
ONE_SECOND = 1




#Configure detailed logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("gene_fetch-up_rank.log")
    ]
)
logger = logging.getLogger()





#Retry decorator for Entrez search functions
def retry(ExceptionToCheck, tries=3, delay=1, backoff=2):
    def deco_retry(f):
        @wraps(f)
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay
            while mtries > 1:
                try:
                    return f(*args, **kwargs)
                except ExceptionToCheck as e:
                    msg = f"{e}, Retrying in {mdelay} seconds..."
                    logger.warning(msg)
                    time.sleep(mdelay)
                    mtries -= 1
                    mdelay *= backoff
            return f(*args, **kwargs)
        return f_retry
    return deco_retry





@sleep_and_retry
@limits(calls=MAX_CALLS_PER_SECOND, period=ONE_SECOND)
@retry((RequestException, RuntimeError), tries=4, delay=1, backoff=2)
def entrez_esearch(**kwargs):
    """
    Wrapper for Entrez.esearch with rate limiting and retry.
    """
    return Entrez.esearch(**kwargs)





@sleep_and_retry
@limits(calls=MAX_CALLS_PER_SECOND, period=ONE_SECOND)
@retry((RequestException, RuntimeError), tries=4, delay=1, backoff=2)
def entrez_efetch(**kwargs):
    """
    Wrapper for Entrez.efetch with rate limiting and retry.
    """
    return Entrez.efetch(**kwargs)






def fetch_protein_seq_by_taxonomy(taxonomy, gene_name, retmax=1000):
    """
    Fetches the longest protein sequence from NCBI based on taxonomy and gene name.
    Iterates from species to order, updating the best sequence found at each level.

    Parameters:
    - taxonomy (iterable): Taxonomic ranks (order, family, genus, species).
    - gene_name (str): Name of the gene to search for.
    - retmax (int): Maximum number of records to retrieve at each level.

    Returns:
    - tuple: (list of SeqRecord objects, matched_rank as str, ncbi_taxonomy as list)
    """
    best_rank = None
    best_record = None
    ncbi_taxonomy = []

    rank_priority = ["species", "genus", "family", "order"]  
    length_threshold = 500

    #Iterate over ranks from species to order
    for rank_name in rank_priority:
        rank = taxonomy.get(rank_name)

        if rank:  #Check if the current rank is not empty
            search_term = f"{gene_name}[Gene] AND {rank}[Organism] AND 500:1000[SLEN]"
            logger.debug(f"Search term for {rank_name}: {search_term}")

            try:
                total_ids = []
                search_handle = entrez_esearch(db="protein", term=search_term, retmax=retmax)  
                search_results = Entrez.read(search_handle)
                search_handle.close()

                logger.info(f"Total records found for {rank_name} ({rank}): {search_results['Count']}")

                if not search_results["IdList"]:
                    logger.warning(f"No sequences found for {rank_name.capitalize()}: {rank}. Trying next level...")
                    continue

                total_ids.extend(search_results["IdList"])
                logger.info(f"Fetched {len(total_ids)} IDs for {rank_name}: {rank}")

                #Fetch records in batches
                for start in range(0, len(total_ids), retmax):
                    batch_ids = total_ids[start:start + retmax]
                    logger.debug(f"Fetching IDs: {batch_ids}")
                    fetch_handle = entrez_efetch(db="protein", id=batch_ids, rettype="gp", retmode="text")
                    records = list(SeqIO.parse(fetch_handle, "genbank"))
                    fetch_handle.close()

                    if records:
                        #Check each record, update best_record if necessary and fetch full NCBI taxonomy from record
                        for record in records:
                            record_length = len(record.seq)
                            fetched_taxonomy = record.annotations.get("taxonomy", [])  
                            ncbi_taxonomy_str = ", ".join(fetched_taxonomy)  

                            logger.info(f"Fetched NCBI taxonomy for {record.id}: {ncbi_taxonomy_str}")

                            if record_length >= length_threshold:
                                if best_record is None:
                                    best_record = record
                                    best_rank = f"{rank_name.capitalize()}: {rank}"
                                    ncbi_taxonomy = fetched_taxonomy  
                                    logger.info(f"First sequence found at {best_rank}: {record.id} with length {record_length}")
                                elif record_length > len(best_record.seq):
                                    best_record = record
                                    best_rank = f"{rank_name.capitalize()}: {rank}"
                                    ncbi_taxonomy = fetched_taxonomy  #Update NCBI taxonomy if a better record is found
                                    logger.info(f"Replacing sequence with longer sequence at {best_rank}: {record.id}, length: {record_length}")
                            else:
                                logger.info(f"Sequence at {rank_name}: {record.id} is too short (length {record_length}). Skipping...")

                    #Only stop if a valid sequence has been found (length > 500)
                    if best_record is not None and len(best_record.seq) >= length_threshold:
                        logger.info(f"Stopping search as sequence found at {best_rank}")
                        break

            except Exception as e:
                logger.error(f"Error fetching data for {rank_name.capitalize()} ({rank}): {e}")

        #If a valid sequence is found, break out of the loop
        if best_record is not None:
            break  #Stop searching after the first valid sequence

    if best_record is None:
        logger.warning(f"No sequences found for taxonomy levels: {taxonomy}.")
        return [], None, []

    logger.info(f"Final best sequence found at {best_rank}: {best_record.id}, length: {len(best_record.seq)}")
    return [best_record], best_rank, ncbi_taxonomy  #Return fetched NCBI taxonomy





def validate_taxonomy(ncbi_taxonomy):
    """
    Validate the NCBI taxonomy at the phylum level.
    
    Parameters:
    - ncbi_taxonomy (list): List of taxonomic ranks from NCBI.
    
    Returns:
    - bool: True if the taxonomy belongs to Arthropoda, Mollusca, or Annelida, else False.
    """
    valid_phyla = {"Arthropoda", "Mollusca", "Annelida"}

    if not ncbi_taxonomy:
        logger.warning("NCBI taxonomy is empty.")
        return False  #Empty taxonomy

    #Split taxonomy string into list
    taxonomy_list = [rank.strip() for rank in ncbi_taxonomy]

    #Check for at least 4 ranks to access phylum (index 3 = phylum)
    if len(taxonomy_list) >= 4:
        phylum = taxonomy_list[3]  
        if phylum in valid_phyla:
            logger.info(f"Validated taxonomy with Phylum: {phylum} as valid.")
            return True
        else:
            logger.warning(f"Taxonomy Phylum: {phylum} is not valid.")
            return False
    else:
        logger.warning("Taxonomy list does not have enough ranks to determine phylum.")
        return False




def determine_delimiter(file_path):
    """
    Determine the delimiter for the CSV file based on its extension.

    Parameters:
    - file_path (str): Path to the input file.

    Returns:
    - str: The delimiter used in the CSV file.
    """
    if file_path.endswith('.tsv'):
        return '\t'
    return ','  #Default to comma




def read_existing_summary(summary_output_file):
    """
    Read existing summary entries from the output summary CSV file and determine if the file is empty.

    Parameters:
    - summary_output_file (str): Path to the summary output file.

    Returns:
    - tuple: (list of existing summary entries, is_file_empty)
             existing_entries: list of entries from the file.
             is_file_empty: boolean indicating if the file is empty.
    """
    existing_entries = []
    
    #Check if file exists
    if os.path.exists(summary_output_file):
        #Check if file is empty
        if os.stat(summary_output_file).st_size == 0:
            return existing_entries, True  
        else:
            #If file not empty, read existing entries
            with open(summary_output_file, 'r', newline='') as summary_csv:
                reader = csv.DictReader(summary_csv)
                existing_entries = list(reader)
            return existing_entries, False  
    
    return existing_entries, True  





def main():
    if len(sys.argv) != 5:
        usage()
        sys.exit(1)

    input_file = sys.argv[1]
    gene_name = sys.argv[2]
    user_output_directory = sys.argv[3]
    samples_file = sys.argv[4]

    #Set retmax global variable (i.e. how many records to search at each rank)
    retmax = 1000

    #Check inputs exist
    if not os.path.exists(input_file):
        logger.error(f"The file '{input_file}' does not exist.")
        usage()
        sys.exit(1)

    if not os.path.exists(samples_file):
        logger.error(f"The file '{samples_file}' does not exist.")
        usage()
        sys.exit(1)

    #Determine input delimiter based on file extension
    delimiter = determine_delimiter(input_file)
    logger.info(f"Using delimiter '{delimiter}' for taxonomy file.")

    #Check/Create output directory
    if not os.path.exists(user_output_directory):
        try:
            os.makedirs(user_output_directory)
            logger.info(f"Created output directory: {user_output_directory}")
        except Exception as e:
            logger.error(f"Failed to create output directory '{user_output_directory}': {e}")
            sys.exit(1)

    #Load Process IDs from samples.csv
    valid_process_ids = set()
    try:
        with open(samples_file, newline='') as samples_csv:
            reader = csv.DictReader(samples_csv)
            for row in reader:
                valid_process_ids.add(row['ID'].strip())
        logger.info(f"Loaded {len(valid_process_ids)} Process IDs from '{samples_file}'.")
    except Exception as e:
        logger.error(f"Error reading samples file '{samples_file}': {e}")
        sys.exit(1)

    summary_output = []

    #Process input taxonomy file
    try:
        with open(input_file, newline='') as taxonomy_file:
            reader = csv.DictReader(taxonomy_file, delimiter=delimiter)

            #Create a list of field names in lowercase for consistency
            fieldnames_lower = [field.lower() for field in reader.fieldnames]

            #Check for existing summary entries
            summary_output_file = os.path.join(
                user_output_directory,
                f"{os.path.splitext(os.path.basename(input_file))[0]}_gene_fetch_summary.csv"
            )
            existing_entries = read_existing_summary(summary_output_file)

            for row in reader:
                row_lower = {k.lower(): v for k, v in row.items()}

                process_id = row_lower.get('process id', '').strip()  

                if not process_id:
                    logger.warning("Missing 'Process ID' in a row. Skipping.")
                    continue

                logger.info(f"\nProcessing Process ID: {process_id}")

                #Skip unmatched Process IDs
                if process_id not in valid_process_ids:
                    logger.info(f"Skipping Process ID: {process_id} (not in samples.csv)")
                    continue

                #Extract taxonomic ranks from input taxonomy file
                taxonomy = {
                    'order': row_lower.get('order', '').strip(),
                    'family': row_lower.get('family', '').strip(),
                    'genus': row_lower.get('genus', '').strip(),
                    'species': row_lower.get('species', '').strip()
                }

                logger.info(f"Taxonomy: order='{taxonomy['order']}', family='{taxonomy['family']}', genus='{taxonomy['genus']}', species='{taxonomy['species']}'")

                #Fetch records using fetch_protein_seq_by_taxonomy function
                records, matched_rank, ncbi_taxonomy = fetch_protein_seq_by_taxonomy(taxonomy, gene_name, retmax=retmax)

                #Validate taxonomy
                tax_validated = validate_taxonomy(ncbi_taxonomy)

                #Write Sequences to fasta file for each Process ID
                if records:
                    for record in records:
                        accession_number = record.id
                        record.id = process_id  #Using process_id here
                        record.description = ""
                        sequence_length = len(record.seq)  #Get the length of the sequence

                        fasta_filename = f"{process_id}.fasta"
                        fasta_filepath = os.path.abspath(os.path.join(user_output_directory, fasta_filename))

                        #Append the new sequence to the FASTA file (or create if it doesn't exist)
                        try:
                            with open(fasta_filepath, "a") as output_handle: 
                                SeqIO.write([record], output_handle, "fasta")
                            logger.info(f"Written FASTA file: {fasta_filepath}")
                        except Exception as e:
                            logger.error(f"Error writing FASTA file '{fasta_filepath}': {e}")
                            continue

                        #Add info to summary output
                        summary_output.append({
                            'process_id': process_id,
                    	   'matched_term': matched_rank if matched_rank else "None",
                    	   'accession_number': accession_number,
                            'length': sequence_length,
                            'reference_name': process_id,
                            'reference_path': fasta_filepath,
                            'tax_validated': 'true' if tax_validated else 'false',
                            'ncbi_taxonomy': ", ".join(ncbi_taxonomy)  #Include NCBI taxonomy in CSV output
                        })
                        
                        logger.info(f"Fetched sequence accession: {accession_number}, length: {sequence_length} for {matched_rank}")
                else:
                    summary_output.append({
                        'process_id': process_id,
                        'matched_term': "None",
                        'accession_number': "None",
                        'length': "None",
                        'reference_name': process_id,
                        'reference_path': "None",
                        'tax_validated': "false",
                        'ncbi_taxonomy': "None"
                    })

                    logger.warning(f"No protein sequence found for Process ID {process_id}.")


            #Write summary output csv
            try:

                #Read existing summary entries and check if the file is empty
                existing_entries, file_empty = read_existing_summary(summary_output_file)

                with open(summary_output_file, "a" if not file_empty else "w", newline='') as summary_csv:
                    fieldnames = [
                        'process_id',
                        'matched_term',
                        'accession_number',
                        'length',
                        'reference_name',
                        'reference_path',
                        'tax_validated',
                        'ncbi_taxonomy'
                    ]
                    writer = csv.DictWriter(summary_csv, fieldnames=fieldnames)

                    #If the file is empty, write the header
                    if file_empty:
                        writer.writeheader()

                    #Append new rows to the file
                    writer.writerows(summary_output)

                logger.info(f"\nSummary output saved to {summary_output_file}.")

            except Exception as e:
                logger.error(f"Error writing summary CSV '{summary_output_file}': {e}")
                sys.exit(1)

    except Exception as e:
        logger.error(f"Error processing taxonomy file '{input_file}': {e}")
        sys.exit(1)



if __name__ == "__main__":
    main()
