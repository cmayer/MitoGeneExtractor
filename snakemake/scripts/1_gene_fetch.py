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
    print("""\
    Usage: python 1_gene_fetch.py <gene_name> <output_directory> <samples.csv>

    <gene_name>: Name of gene to search for in NCBI RefSeq database (e.g., cox1).
    <output_directory>: Path to directory to save output files (will save .fasta files and summary CSV in this directory).
                        The directory will be created if it does not exist.
    <samples.csv>: Path to input CSV file containing TaxIDs to search (with column 'taxid').

    Example:
    python 1_gene_fetch.py cox1 ./protein_references/output_dir path/to/samples.csv
    """)




# Set email and API key
Entrez.email = "b.price@nhm.ac.uk"  
Entrez.api_key = "82df5a6f5cf735302d3cf1fcf48b206cfe09"  




# Constants for rate limiting (Entrez API limits with key = 10/sec)
MAX_CALLS_PER_SECOND = 10
ONE_SECOND = 1




# Configure detailed logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("gene_fetch.log")
    ]
)
logger = logging.getLogger()




# Retry decorator for Entrez search functions
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




def fetch_taxonomy_by_taxid(taxid):
    """
    Fetches the full taxonomy and rank for a given TaxID from the NCBI database.

    Parameters:
    - taxid (str): NCBI taxonomy ID.

    Returns:
    - tuple: Full taxonomy (list), rank of the TaxID (str)
    """
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        if records:
            record = records[0]
            
            # Extract NCBI lineage
            taxonomy = record.get("Lineage", "").split("; ")
            
            # Extract rank of given TaxID & include at end of NCBI lineage
            rank = record.get("Rank", "No rank available")
            scientific_name = record.get("ScientificName", "")
            if scientific_name and scientific_name not in taxonomy:
                taxonomy.append(scientific_name)

            # Return taxonomy and rank
            return taxonomy, rank
        else:
            logger.warning(f"No taxonomy found for TaxID: {taxid}")
            return [], "No rank available"
    except Exception as e:
        logger.error(f"Error fetching taxonomy for TaxID {taxid}: {e}")
        return [], "No rank available"



def fetch_protein_seq_by_taxid(taxid, gene_name, retmax=1000):
    """
    Fetches the longest protein sequence from NCBI based on taxid and gene name,
    moving up the taxonomic hierarchy if necessary.

    Parameters:
    - taxid (str): NCBI taxonomy ID.
    - gene_name (str): Name of the gene to search for.
    - retmax (int): Maximum number of records to retrieve.

    Returns:
    - tuple: (best_record, best_taxonomy as list, matched_rank as str)
    """
    best_record = None
    best_taxonomy = []
    length_threshold = 500

    # Get taxonomy for initial input taxid
    full_taxonomy, matched_rank = fetch_taxonomy_by_taxid(taxid)

    # Construct search_term based on taxid, gene name, and length thresholds
    search_term = f"{gene_name}[Gene] AND txid{taxid}[Organism] AND 500:1000[SLEN]"
    logger.info(f"Constructed search term for TaxID {taxid}: {search_term}")

    try:
        total_ids = []
        search_handle = entrez_esearch(db="protein", term=search_term, retmax=retmax)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        logger.info(f"Total records found for TaxID {taxid}: {search_results['Count']}")

        if not search_results["IdList"]:
            logger.warning(f"No sequences found for TaxID: {taxid}. Attempting to traverse up the ranks...")
            # If no records found, traverse up the lineage
            for rank in reversed(full_taxonomy):
                logger.info(f"Searching for rank: {rank}")
                search_term = f"{gene_name}[Gene] AND {rank}[Organism] AND 500:1000[SLEN]"
                search_handle = entrez_esearch(db="protein", term=search_term, retmax=retmax)
                search_results = Entrez.read(search_handle)
                search_handle.close()

                if search_results["IdList"]:
                    total_ids = search_results["IdList"]
                    logger.info(f"Total records found for {rank}: {len(total_ids)}")

                    # Fetch records in batches
                    for start in range(0, len(total_ids), retmax):
                        batch_ids = total_ids[start:start + retmax]
                        fetch_handle = entrez_efetch(db="protein", id=batch_ids, rettype="gp", retmode="text")
                        records = list(SeqIO.parse(fetch_handle, "genbank"))
                        fetch_handle.close()

                        for record in records:
                            record_length = len(record.seq)
                            fetched_taxonomy = record.annotations.get("taxonomy", [])
                            ncbi_taxonomy_str = ", ".join(fetched_taxonomy)

                            if record_length >= length_threshold:
                                if best_record is None or record_length > len(best_record.seq):
                                    best_record = record
                                    best_taxonomy = fetched_taxonomy
                                    logger.info(f"Longest sequence found so far: {record.id} (length {record_length}) with taxonomy: {ncbi_taxonomy_str}")
                            else:
                                logger.info(f"Sequence at {rank}: {record.id} is too short (length {record_length}). Skipping...")

                if best_record is not None:
                    logger.info(f"Stopping search as sequence found for {rank}")
                    break

        else:
            total_ids.extend(search_results["IdList"])
            logger.info(f"Fetched {len(total_ids)} IDs for TaxID: {taxid}")

            # Fetch records in batches
            for start in range(0, len(total_ids), retmax):
                batch_ids = total_ids[start:start + retmax]
                fetch_handle = entrez_efetch(db="protein", id=batch_ids, rettype="gp", retmode="text")
                records = list(SeqIO.parse(fetch_handle, "genbank"))
                fetch_handle.close()

                if records:
                    for record in records:
                        record_length = len(record.seq)
                        fetched_taxonomy = record.annotations.get("taxonomy", [])
                        ncbi_taxonomy_str = ", ".join(fetched_taxonomy)

                        if record_length >= length_threshold:
                            if best_record is None or record_length > len(best_record.seq):
                                best_record = record
                                best_taxonomy = fetched_taxonomy
                                logger.info(f"Longest sequence found so far: {record.id} (length {record_length}) with taxonomy: {ncbi_taxonomy_str}")
                        else:
                            logger.info(f"Sequence at TaxID {taxid}: {record.id} is too short (length {record_length}). Skipping...")

    except Exception as e:
        logger.error(f"Error fetching data for TaxID {taxid}: {e}")

    if best_record is None:
        logger.warning(f"No sequences found for TaxID {taxid} after searching all ranks.")
    
    return best_record, best_taxonomy, matched_rank  



def main():
    if len(sys.argv) != 4:
        usage()
        sys.exit(1)

    gene_name = sys.argv[1]
    output_directory = sys.argv[2]
    samples_csv = sys.argv[3]

    # Create output_directory if it does not exist
    os.makedirs(output_directory, exist_ok=True)


    # Check if protein_references.csv already exists    
    output_csv_path = os.path.join(output_directory, "protein_references.csv")

    if os.path.exists(output_csv_path):
        logger.info(f"Output CSV file '{output_csv_path}' already exists. Skipping processing.")
        return

    # Prepare the output CSV file with headers
    with open(output_csv_path, mode='w', newline='') as summary_file:
        summary_writer = csv.writer(summary_file)
        summary_writer.writerow(['process_id', 'taxid', 'accession_number', 'sequence_length', 
                                 'matched_rank', 'ncbi_taxonomy', 'reference_name', 'reference_path'])
        summary_file.flush()

        # Read input samples.csv file and process each TaxID
        with open(samples_csv, newline='') as samples_file:
            reader = csv.DictReader(samples_file)
            process_ids = [row['ID'] for row in reader]  
            num_process_ids = len(process_ids) 

            logger.info(f"Detected {num_process_ids} process IDs in input samples.csv.")

            # Reset reader to iterate over again
            samples_file.seek(0)
            next(reader)  

            for row in reader:
                taxid = row['taxid']
                process_id = row['ID'] 
                logger.info(f"Processing TaxID: {taxid} for Process ID: {process_id}")

                # Fetch full taxonomy and rank for logging
                full_taxonomy, matched_rank = fetch_taxonomy_by_taxid(taxid)
                logger.info(f"Full taxonomy for TaxID {taxid}: {', '.join(full_taxonomy)}")

                # Define output sequence file path & check if it already exists in output_directory
                output_seq_path = os.path.join(output_directory, f"{process_id}.fasta")

                if os.path.exists(output_seq_path):
                    logger.info(f"Sequence file '{output_seq_path}' already exists. Skipping fetching sequence.")
                    continue

                # Fetch protein sequence by TaxID and gene name and return metadata
                best_record, best_taxonomy, matched_rank = fetch_protein_seq_by_taxid(taxid, gene_name)

                # If a valid record was found, write it to the output files
                if best_record:
                    # Change header to the process_id
                    best_record.id = process_id  
                    best_record.name = ""  
                    best_record.description = ""  

                    SeqIO.write(best_record, output_seq_path, "fasta")
                    logger.info(f"Written sequence to '{output_seq_path}'.")

                    # Extract reference_name and reference_path
                    reference_name = process_id  
                    reference_path = os.path.abspath(output_seq_path)  

                    # Write fields to protein_references.csv
                    summary_writer.writerow([
                        process_id,
                        taxid,
                        best_record.id,
                        len(best_record.seq),
                        matched_rank,
                        "; ".join(best_taxonomy) if best_taxonomy else "",
                        reference_name,
                        reference_path
                    ])
                    summary_file.flush()  
                else:
                    logger.warning(f"No valid record found for TaxID {taxid}.")

if __name__ == "__main__":
    main()
