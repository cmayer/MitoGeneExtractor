###Grabs closest cox1 protein reference sequence for each sample using BOLD-downloaded taxonomic rankings, 
###outputs the cox1.fasta seqs to a user specified dir and outputs a protein_reference.csv


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
    python 1_gene_fetch.py taxonomy.csv cox1 ./protein_references samples.csv
    """)


#Set your email and API key
Entrez.email = "d.parsons@nhm.ac.uk"  
Entrez.api_key = "1866f9734a06f26bc5895a84387542ac9308"  



#Constants for rate limiting (Entrez API limits with key = 10/sec)
MAX_CALLS_PER_SECOND = 10
ONE_SECOND = 1



#Configure detailed logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("gene_fetch.log")
    ]
)
logger = logging.getLogger()








#Retry decorator for entrez search functions
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

    #Iterate over ranks from order to species
    for rank, rank_name in zip(taxonomy, rank_priority):
        if rank:
            search_term = f"{gene_name}[Gene] AND {rank}[Organism] AND refseq[filter]"
            logger.debug(f"Search term for {rank_name}: {search_term}")
            try:
                search_handle = entrez_esearch(db="protein", term=search_term, retmax=retmax)
                search_results = Entrez.read(search_handle)
                search_handle.close()

                if search_results["IdList"]:
                    ids = search_results["IdList"]
                    logger.info(f"Fetching {gene_name} sequence: {ids} for {rank_name.capitalize()}: {rank}")
                    fetch_handle = entrez_efetch(db="protein", id=ids, rettype="gb", retmode="text")
                    records = list(SeqIO.parse(fetch_handle, "genbank"))
                    fetch_handle.close()

                    #Update best records and rank if found a better match (lower rank)
                    best_records = records
                    best_rank = f"{rank_name.capitalize()}: {rank}"
                    found_sequence = True 
                else:
                    logger.warning(f"No sequences found for {rank_name.capitalize()}: {rank}. Trying next level...")
                    if found_sequence:
                        break  

            except Exception as e:
                logger.error(f"Error fetching data for {rank_name.capitalize()}: {rank}: {e}")


    if not best_records:
        logger.warning(f"No sequences found for taxonomy {taxonomy} and gene {gene_name}.")

    return best_records, best_rank




@retry((RequestException, RuntimeError), tries=4, delay=1, backoff=2)
def fetch_tax_id_by_name(organism_name):
    """
    Fetches the taxonomic ID for a given organism name from NCBI.

    Parameters:
    - organism_name (str): Scientific name of the organism.

    Returns:
    - str or None: Taxonomic ID if found, else None.
    """


    try:
        search_handle = entrez_esearch(db="taxonomy", term=organism_name, retmax=1)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        if search_results["IdList"]:
            tax_id = search_results["IdList"][0]
            logger.debug(f"Tax ID for {organism_name}: {tax_id}")
            return tax_id
        else:
            logger.warning(f"No Tax ID found for organism: {organism_name}")
            return None
    except Exception as e:
        logger.error(f"Error fetching Tax ID for {organism_name}: {e}")
        return None




@retry((RequestException, RuntimeError), tries=4, delay=1, backoff=2)
def fetch_taxonomy_by_id(tax_id):
    """
    Retrieves full taxonomy information for a given taxonomic ID.

    Parameters:
    - tax_id (str): Taxonomic ID.

    Returns:
    - dict or None: Taxonomy information if found, else None.
    """


    try:
        handle = entrez_efetch(db="taxonomy", id=tax_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        if records:
            logger.debug(f"Fetched taxonomy for Tax ID {tax_id}")
            return records[0]
        else:
            logger.warning(f"No taxonomy records found for Tax ID {tax_id}")
            return None
    except Exception as e:
        logger.error(f"Error fetching taxonomy for tax_id {tax_id}: {e}")
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
        logger.warning("Fetched taxonomy is None, cannot validate.")
        return False

    fetched_order = next(
        (lineage['ScientificName'] for lineage in fetched_taxonomy.get('LineageEx', []) if lineage['Rank'].lower() == 'order'),
        None
    )

    if input_order and fetched_order:
        if input_order.lower() != fetched_order.lower():
            logger.warning(f"Order mismatch: input '{input_order}' vs fetched '{fetched_order}'")
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


    #Filter out clades and start from 'Kingdom'
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
                logger.debug("Detected ';' as delimiter for CSV file.")
                return ';'
            else:
                logger.debug("Detected ',' as delimiter for CSV file.")
                return ','
    elif ext.lower() == '.tsv':
        logger.debug("Detected tab delimiter for TSV file.")
        return '\t'
    else:
        logger.error("Unsupported file extension. Please provide a .csv or .tsv file for taxonomy.")
        sys.exit(1)




def main():
    if len(sys.argv) != 5:
        usage()
        sys.exit(1)

    input_file = sys.argv[1]
    gene_name = sys.argv[2]
    user_output_directory = sys.argv[3]
    samples_file = sys.argv[4]



    #Check Inputs exist
    if not os.path.exists(input_file):
        logger.error(f"The file '{input_file}' does not exist.")
        usage()
        sys.exit(1)

    if not os.path.exists(samples_file):
        logger.error(f"The file '{samples_file}' does not exist.")
        usage()
        sys.exit(1)


    #Determine input delimitera based on file extension
    delimiter = determine_delimiter(input_file)
    logger.info(f"Using delimiter '{delimiter}' for taxonomy file.")


    #Check/Create output dire
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

            fieldnames_lower = [field.lower() for field in reader.fieldnames]

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


                #Fetch protein sequences using function
                records, matched_rank = fetch_protein_seq_by_taxonomy(taxonomy.values(), gene_name, retmax=1)


                #Fetch NCBI taxonomy and validate order to mitigate homonym matches (experimental)
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


                #Write Sequences to fasta file for each Process ID
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
                            logger.info(f"Written FASTA file: {fasta_filepath}")
                        except Exception as e:
                            logger.error(f"Error writing FASTA file '{fasta_filepath}': {e}")
                            continue

                        if fetched_taxonomy:
                            lineage = fetched_taxonomy.get('LineageEx', [])
                            taxonomy_str = format_taxonomy(lineage)
                        else:
                            taxonomy_str = "Not available"


                        #Add info to summary output
                        summary_output.append({
                            'process_id': process_id,
                            'matched_term': matched_rank if matched_rank else "None",
                            'accession_number': accession_number,
                            'reference_name': process_id,
                            'reference_path': fasta_filepath,
                            'matched_validated': 'true' if matched_validated else 'false',
                            'ncbi_taxonomy': taxonomy_str  # Add filtered taxonomy to output
                        })

                        logger.info(f"Fetched sequence accession: {accession_number} for {matched_rank}")

    except Exception as e:
        logger.error(f"Error processing taxonomy file '{input_file}': {e}")
        sys.exit(1)


    #Write summary output csv
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
        logger.info(f"\nSummary output saved to {summary_output_file}.")
    except Exception as e:
        logger.error(f"Error writing summary CSV '{summary_output_file}': {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
