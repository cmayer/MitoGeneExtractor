import csv
import sys
import os
import logging
import time
from Bio import Entrez, SeqIO
from functools import wraps
from ratelimit import limits, sleep_and_retry
from requests.exceptions import RequestException
from time import sleep
from random import uniform
from urllib.error import HTTPError





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
Entrez.email = "#######################"  
Entrez.api_key = "#####################"  




# Constants for rate limiting (Entrez API limits with key = 10/sec)
MAX_CALLS_PER_SECOND = 10
ONE_SECOND = 1




# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger()



def retry(ExceptionToCheck, tries=4, initial_delay=5, backoff=2, max_delay=120):
    """
    Retry decorator with exponential backoff and longer delays.
    
    Parameters:
    - tries: number of times to try (not retry) before giving up
    - initial_delay: initial delay between retries in seconds
    - backoff: multiplier applied to delay between retries
    - max_delay: maximum delay between retries in seconds
    """
    def deco_retry(f):
        @wraps(f)
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, initial_delay
            while mtries > 1:
                try:
                    return f(*args, **kwargs)
                except ExceptionToCheck as e:
                    # Add specific handling for common NCBI errors
                    if isinstance(e, HTTPError):
                        if e.code == 400:  # Bad Request
                            logger.error(f"Bad request error: {e}. Skipping...")
                            return None
                        elif e.code == 429:  # Too Many Requests
                            mdelay = min(mdelay * 3, max_delay)  # Triple delay for rate limit errors
                            
                    msg = f"{str(e)}, Retrying in {mdelay} seconds..."
                    logger.warning(msg)
                    
                    # Add jitter to avoid synchronized retries
                    actual_delay = mdelay + uniform(-0.1 * mdelay, 0.1 * mdelay)
                    sleep(actual_delay)
                    
                    mtries -= 1
                    mdelay = min(mdelay * backoff, max_delay)
            try:
                return f(*args, **kwargs)
            except ExceptionToCheck as e:
                logger.error(f"Final attempt failed: {str(e)}")
                return None
        return f_retry
    return deco_retry


@retry((HTTPError, RuntimeError, IOError), tries=4, initial_delay=5, backoff=2, max_delay=120)
@sleep_and_retry
@limits(calls=10, period=1.1)  # Slightly longer period to be safe
def entrez_esearch(**kwargs):
    """
    Wrapper for Entrez.esearch with enhanced error handling.
    """
    try:
        handle = Entrez.esearch(**kwargs)
        result = Entrez.read(handle)
        handle.close()
        return result
    except RuntimeError as e:
        if "Search Backend failed" in str(e):
            sleep(uniform(2, 4))  # Additional delay for backend failures
        raise
    except Exception as e:
        logger.error(f"Unexpected error in entrez_esearch: {str(e)}")
        raise

@retry((HTTPError, RuntimeError, IOError), tries=4, initial_delay=5, backoff=2, max_delay=120)
@sleep_and_retry
@limits(calls=10, period=1.1)
def entrez_efetch(**kwargs):
    """
    Wrapper for Entrez.efetch with enhanced error handling.
    """
    try:
        handle = Entrez.efetch(**kwargs)
        return handle
    except HTTPError as e:
        if e.code == 400:  # Bad Request
            logger.error(f"Bad request error for query parameters: {kwargs}")
            return None
        if e.code == 429:  # Too Many Requests
            logger.warning("Rate limit exceeded, backing off...")
            sleep(uniform(2, 4))
        raise
    except RuntimeError as e:
        if "Search Backend failed" in str(e):
            sleep(uniform(2, 4))
        raise
    except Exception as e:
        logger.error(f"Unexpected error in entrez_efetch: {str(e)}")
        raise




def fetch_taxonomy_by_taxid(taxid):
    """
    Fetches taxonomy with rank information for a given TaxID.
    Maintains logging while handling errors gracefully.
    """
    try:
        taxid = str(taxid).strip()
        
        # Add delay before taxonomy fetch to help prevent rate limiting
        sleep(uniform(0.2, 0.5))
        
        try:
            handle = entrez_efetch(db="taxonomy", id=taxid, retmode="xml")
            if handle is None:  # Handle was None due to Bad Request
                logger.warning(f"Could not fetch taxonomy data for TaxID {taxid}")
                return [], {}, "no rank", {}
                
            records = Entrez.read(handle)
            handle.close()

            if records and len(records) > 0:
                record = records[0]
                logger.info(f"Successfully fetched taxonomy record for TaxID {taxid}")
                
                # Get the regular taxonomy list
                taxonomy = record.get("Lineage", "").split("; ")
                scientific_name = record.get("ScientificName")
                if scientific_name and scientific_name not in taxonomy:
                    taxonomy.append(scientific_name)
                
                # Initialize dictionaries for ranks and taxids
                taxon_ranks = {}
                taxon_ids = {}
                
                # Process LineageEx for detailed rank information
                lineage_ex = record.get("LineageEx", [])
                for node in lineage_ex:
                    if "ScientificName" in node:
                        name = node["ScientificName"]
                        if "Rank" in node:
                            taxon_ranks[name] = node["Rank"]
                        if "TaxId" in node:
                            taxon_ids[name] = node["TaxId"]
                
                # Add the current taxon's information
                if scientific_name:
                    if "Rank" in record:
                        taxon_ranks[scientific_name] = record["Rank"]
                    taxon_ids[scientific_name] = taxid
                
                return taxonomy, taxon_ranks, record.get("Rank", "no rank"), taxon_ids
            else:
                logger.warning(f"No taxonomy records found for TaxID: {taxid}")
                return [], {}, "no rank", {}
                
        except Exception as parse_error:
            if handle:
                try:
                    handle.close()
                except:
                    pass
            logger.error(f"Error parsing taxonomy response for TaxID {taxid}: {parse_error}")
            return [], {}, "no rank", {}
            
    except Exception as e:
        logger.error(f"Error fetching taxonomy for TaxID {taxid}: {e}")
        return [], {}, "no rank", {}

# And in fetch_protein_seq_by_taxid, modify the logging duplicate check:
    if full_taxonomy:
        if not hasattr(fetch_protein_seq_by_taxid, '_last_logged_taxid') or \
           fetch_protein_seq_by_taxid._last_logged_taxid != taxid:
            logger.info(f"Full taxonomy: {', '.join(full_taxonomy)}")
            fetch_protein_seq_by_taxid._last_logged_taxid = taxid


@retry((HTTPError, RuntimeError, IOError), tries=5, initial_delay=10, backoff=2, max_delay=120)
@sleep_and_retry
@limits(calls=8, period=1.1)  # Slightly more conservative rate limiting
def entrez_efetch(**kwargs):
    """
    Wrapper for Entrez.efetch with enhanced error handling.
    """
    try:
        handle = Entrez.efetch(**kwargs)
        return handle
    except HTTPError as e:
        if e.code == 400:  # Bad Request
            logger.warning(f"Bad request error for parameters: {kwargs}")
            return None
        if e.code == 429:  # Too Many Requests
            logger.warning("Rate limit exceeded, backing off...")
            sleep(uniform(5, 10))
        raise
    except RuntimeError as e:
        if "Search Backend failed" in str(e):
            sleep(uniform(5, 10))  # Longer delay for backend failures
        raise
    except Exception as e:
        logger.error(f"Unexpected error in entrez_efetch: {str(e)}")
        raise






def get_taxid_for_name(name, context_taxonomy=None):
    """
    Helper function to get taxid for a scientific name, using taxonomic context to ensure correct taxid
    """
    try:
        # Start with base search for the name
        search_term = f'"{name}"[Scientific Name]'
        if context_taxonomy:
            context_terms = [f'"{term}"[Organism]' for term in context_taxonomy if term != name]
            if context_terms:
                search_term = f'{search_term} AND ({" AND ".join(context_terms)})'
        
        logger.info(f"Searching taxonomy database with term: {search_term}")
        
        try:
            result = entrez_esearch(db="taxonomy", term=search_term)
            
            if result and "IdList" in result and result["IdList"]:
                found_taxid = result["IdList"][0]
                
                # Verify this is the correct taxid by checking its full taxonomy
                verify_handle = entrez_efetch(db="taxonomy", id=found_taxid, retmode="xml")
                verify_records = Entrez.read(verify_handle)
                verify_handle.close()
                
                if verify_records:
                    lineage = verify_records[0].get("Lineage", "")
                    logger.info(f"Found taxid {found_taxid} for {name} with lineage: {lineage}")
                    
                    # Check if this is the correct taxonomic context
                    if context_taxonomy and not any(term in lineage for term in context_taxonomy):
                        logger.warning(f"Found taxid {found_taxid} for {name} but it doesn't match taxonomic context")
                        return None
                    return found_taxid
            else:
                logger.warning(f"No taxid found for name: {name}")
                return None
                
        except RuntimeError as e:
            if "Search Backend failed" in str(e):
                logger.warning(f"Search backend failed for {name}, waiting before retry...")
                sleep(uniform(2, 4))  # Increased delay for backend failures
                raise
            else:
                logger.error(f"Runtime error searching for {name}: {e}")
                return None
                
    except Exception as e:
        logger.error(f"Error getting taxid for {name}: {e}")
        return None





def validate_taxonomy(fetched_taxonomy, validation_ranks):
    """
    Helper function to validate taxonomy matches.
    Only validates that the sequence belongs to the correct order.
    
    Parameters:
    - fetched_taxonomy: taxonomy of the sequence being validated
    - validation_ranks: dictionary containing order and family information
    """
    if 'order' not in validation_ranks:
        return True  # If we don't have order information, accept the sequence
        
    order = validation_ranks['order']
    if order not in fetched_taxonomy:
        logger.debug(f"Taxonomy validation failed - Missing order {order}: {fetched_taxonomy}")
        return False
        
    return True




def fetch_protein_seq_by_taxid(taxid, gene_name, retmax=1000, existing_taxonomy=None, existing_rank=None, existing_taxon_ranks=None, existing_taxon_ids=None):
    """
    Fetches the longest protein sequence from NCBI based on taxid and gene name.
    Uses pre-fetched taxonomy information to avoid redundant API calls.
    Only validates sequences against the order rank.
    
    Parameters:
    - taxid (str): NCBI taxonomy ID
    - gene_name (str): Name of gene to search for
    - retmax (int): Maximum number of records to retrieve per batch
    - existing_taxonomy (list): Already fetched taxonomy
    - existing_rank (str): Already fetched rank information
    - existing_taxon_ranks (dict): Already fetched rank mappings
    - existing_taxon_ids (dict): Already fetched taxon IDs
    """
    best_record = None
    best_taxonomy = []
    matched_rank = None
    length_threshold = 500

    # Use existing taxonomy data or fetch if not provided
    if not all([existing_taxonomy, existing_rank, existing_taxon_ranks, existing_taxon_ids]):
        full_taxonomy, taxon_ranks, initial_rank, taxon_ids = fetch_taxonomy_by_taxid(taxid)
        if not full_taxonomy:
            logger.error(f"Could not fetch taxonomy for TaxID {taxid}")
            return None, [], "No taxonomy found"
    else:
        full_taxonomy = existing_taxonomy
        initial_rank = existing_rank
        taxon_ranks = existing_taxon_ranks
        taxon_ids = existing_taxon_ids
        logger.debug(f"Using pre-fetched taxonomy data for TaxID {taxid}")

    # Find order from the input taxonomy using rank information
    validation_ranks = {}
    if full_taxonomy:
        for taxon, rank in taxon_ranks.items():
            if rank == "order":
                validation_ranks['order'] = taxon
                break
                
        logger.info(f"Using order ({validation_ranks.get('order', 'None')}) for validation")

    def process_sequence_batch(batch_ids, context=""):
        """
        Helper function to process a batch of sequence IDs with detailed logging.
        Returns the best record found in this batch.
        """
        local_best_record = None
        try:
            fetch_handle = entrez_efetch(db="protein", id=batch_ids, rettype="gp", retmode="text")
            records = list(SeqIO.parse(fetch_handle, "genbank"))
            fetch_handle.close()
            
            if not records:
                logger.warning(f"{context} - No sequences could be parsed from the response for IDs: {batch_ids}")
                return None
                
            logger.info(f"{context} - Successfully parsed {len(records)} sequences")
            
            for record in records:
                try:
                    record_length = len(record.seq)
                    fetched_taxonomy = record.annotations.get("taxonomy", [])
                    ncbi_taxonomy_str = ", ".join(fetched_taxonomy)
                    
                    logger.debug(f"Processing sequence {record.id}:")
                    logger.debug(f"  Length: {record_length}")
                    logger.debug(f"  Taxonomy: {ncbi_taxonomy_str}")
                    
                    if record_length >= length_threshold:
                        if validate_taxonomy(fetched_taxonomy, validation_ranks):
                            if local_best_record is None or record_length > len(local_best_record.seq):
                                local_best_record = record
                                logger.info(f"{context} - Found valid sequence: {record.id} (length {record_length})")
                                logger.info(f"  Taxonomy: {ncbi_taxonomy_str}")
                        else:
                            logger.debug(f"{context} - Sequence {record.id} failed taxonomy validation")
                    else:
                        logger.debug(f"{context} - Sequence {record.id} too short ({record_length} < {length_threshold})")
                except Exception as e:
                    logger.error(f"{context} - Error processing individual sequence {record.id}: {e}")
                    continue
                    
            return local_best_record
            
        except Exception as e:
            logger.error(f"{context} - Error fetching/processing batch: {e}")
            return None

    try:
        # Construct search term based on gene name
        cox1_synonyms = ['cox1', 'COX1', 'COI', 'COXI']
        if any(syn.lower() == gene_name.lower() for syn in cox1_synonyms):
            gene_search_term = (
                '(cox1[Gene] OR COI[Gene] OR "cytochrome c oxidase subunit 1"[Protein Name] '
                'OR "cytochrome oxidase subunit 1"[Protein Name] OR "cytochrome c oxidase subunit I"[Protein Name] '
                'OR "COX1"[Protein Name] OR "COXI"[Protein Name])'
            )
        else:
            gene_search_term = f'{gene_name}[Gene] OR "{gene_name}"[Protein Name]'
        
        # Try initial search with input taxid
        initial_search_term = f"{gene_search_term} AND txid{taxid}[Organism] AND 500:1000[SLEN]"
        logger.info(f"Initial search term with input taxid {taxid}: {initial_search_term}")

        search_results = entrez_esearch(db="protein", term=initial_search_term, retmax=retmax)
        if not search_results:
            logger.error(f"Search failed for TaxID {taxid}")
            return None, [], "Search failed"
            
        total_count = int(search_results["Count"])
        logger.info(f"Total records found for TaxID {taxid}: {total_count}")

        if total_count > 0 and search_results["IdList"]:
            matched_rank = initial_rank
            logger.info(f"Found sequences at rank: {matched_rank}")

            # Process initial results in batches
            total_ids = search_results["IdList"]
            sequences_processed = 0
            records_found = 0
            
            sleep(uniform(0.5, 1))
            
            for start in range(0, len(total_ids), retmax):
                batch_ids = total_ids[start:start + retmax]
                logger.info(f"Processing batch {start//retmax + 1} of {-(-len(total_ids)//retmax)}")
                
                batch_record = process_sequence_batch(
                    batch_ids, 
                    context=f"Initial search batch {start//retmax + 1}"
                )
                
                if batch_record:
                    records_found += 1
                    if best_record is None or len(batch_record.seq) > len(best_record.seq):
                        best_record = batch_record
                        best_taxonomy = batch_record.annotations.get("taxonomy", [])
                
                sequences_processed += len(batch_ids)
                sleep(uniform(0.5, 1))
            
            logger.info(f"Processed {sequences_processed} sequences, found {records_found} valid records")

        if not best_record:
            logger.warning(f"No valid sequences found for TaxID: {taxid}. Traversing up stored taxonomy...")
            
            # Use the stored taxonomy for traversal, skipping the last entry (current taxon)
            for rank in reversed(full_taxonomy[:-1]):
                current_rank = taxon_ranks.get(rank, 'unknown')
                
                # Stop if we've gone past order level
                if current_rank in ['class', 'subphylum', 'phylum', 'kingdom', 'superkingdom'] or rank == 'cellular organisms':
                    logger.info(f"Reached {current_rank} rank, stopping traversal as we're above order level")
                    break
                
                sleep(uniform(1, 2))
                
                rank_taxid = taxon_ids.get(rank)
                if not rank_taxid:
                    logger.debug(f"No taxid found for rank {rank}, skipping")
                    continue
                    
                logger.info(f"Searching at rank: {current_rank} (taxon: {rank})")
                rank_search_term = f"{gene_search_term} AND txid{rank_taxid}[Organism:exp] AND 500:1000[SLEN]"
                logger.info(f"Searching rank {rank} using its taxid: {rank_taxid}")
                logger.info(f"Search term: {rank_search_term}")
                
                search_results = entrez_esearch(db="protein", term=rank_search_term, retmax=retmax)
                if not search_results:
                    logger.error(f"Search failed for rank {rank}")
                    continue

                found_count = int(search_results["Count"])
                logger.info(f"Found {found_count} sequences for rank {rank}")

                if search_results["IdList"]:
                    total_ids = search_results["IdList"]
                    sequences_processed = 0
                    records_found = 0

                    for start in range(0, len(total_ids), retmax):
                        batch_ids = total_ids[start:start + retmax]
                        logger.info(f"Processing batch {start//retmax + 1} of {-(-len(total_ids)//retmax)} for rank {rank}")
                        
                        batch_record = process_sequence_batch(
                            batch_ids,
                            context=f"Rank {rank} batch {start//retmax + 1}"
                        )
                        
                        if batch_record:
                            records_found += 1
                            if best_record is None or len(batch_record.seq) > len(best_record.seq):
                                best_record = batch_record
                                best_taxonomy = batch_record.annotations.get("taxonomy", [])
                                matched_rank = current_rank
                        
                        sequences_processed += len(batch_ids)
                        sleep(uniform(0.5, 1))
                    
                    logger.info(f"Processed {sequences_processed} sequences at rank {rank}, found {records_found} valid records")

                    if best_record:
                        logger.info(f"Found valid sequence at rank {rank}, stopping traversal")
                        break
                else:
                    logger.info(f"No sequences found for rank {rank}")

    except Exception as e:
        logger.error(f"Error fetching data for TaxID {taxid}: {e}")

    if best_record is None:
        logger.warning(f"No sequences found for TaxID {taxid} after searching all ranks")
        matched_rank = "No match"
    
    return best_record, best_taxonomy, matched_rank




  
def get_process_id_column(header):
    """
    Identifies the process ID column name from possible variations.
    
    Parameters:
    - header (list): List of column names from CSV
    
    Returns:
    - str: Matching column name or None if not found
    """
    valid_names = ['ID', 'process_id', 'Process ID', 'process id', 'Process id', 'PROCESS ID', 'sample', 'SAMPLE', 'Sample']
    for name in valid_names:
        if name in header:
            return name
    return None




def main():
    if len(sys.argv) != 4:
        usage()
        sys.exit(1)

    gene_name = sys.argv[1]
    output_directory = sys.argv[2]
    samples_csv = sys.argv[3]

    # Create output_directory if it does not exist
    os.makedirs(output_directory, exist_ok=True)

    # Configure logging to write to output directory
    log_file = os.path.join(output_directory, "gene_fetch.log")
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
    logger.addHandler(file_handler)
    
    # Create failed_searches.csv in output directory
    failed_searches_path = os.path.join(output_directory, "failed_searches.csv")
    if not os.path.exists(failed_searches_path):
        with open(failed_searches_path, "w", newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['process_id', 'taxid', 'error_type', 'timestamp'])

    # Check if protein_references.csv already exists
    output_csv_path = os.path.join(output_directory, "protein_references.csv")

    # Dictionary to keep track of existing process_ids
    existing_entries = {}

    if os.path.exists(output_csv_path):
        logger.info(f"Output CSV file '{output_csv_path}' already exists. Checking existing entries...")
        with open(output_csv_path, mode='r') as summary_file:
            reader = csv.DictReader(summary_file)
            for row in reader:
                existing_entries[row['process_id']] = row

    # Prepare the output CSV file with append mode
    file_exists = os.path.exists(output_csv_path)
    with open(output_csv_path, mode='a', newline='') as summary_file:
        summary_writer = csv.writer(summary_file)
        if not file_exists:
            summary_writer.writerow(['process_id', 'taxid', 'accession_number', 'sequence_length', 
                                   'matched_rank', 'ncbi_taxonomy', 'reference_name', 'reference_path'])
            summary_file.flush()

        # Read input samples.csv file and process each TaxID
        with open(samples_csv, newline='') as samples_file:
            reader = csv.DictReader(samples_file)
            
            # Identify the process ID column
            process_id_col = get_process_id_column(reader.fieldnames)
            if not process_id_col:
                logger.error("Could not find process ID column. Please use one of: 'ID', 'process_id', 'Process ID', 'process id', 'Process id', 'PROCESS ID', 'sample', 'SAMPLE', 'Sample'")
                sys.exit(1)
                
            # Get list of process IDs
            samples_file.seek(0)
            next(reader)  # Skip header
            process_ids = [row[process_id_col] for row in reader]
            num_process_ids = len(process_ids)

            logger.info(f"Detected {num_process_ids} process IDs in input samples.csv using column '{process_id_col}'")

            # Reset reader to iterate over again
            samples_file.seek(0)
            next(reader)  # Skip header

            for row in reader:
                try:
                    taxid = row['taxid']
                    process_id = row[process_id_col]

                    # Skip if already processed
                    if process_id in existing_entries:
                        logger.info(f"Process ID {process_id} already has an entry in protein_references.csv. Skipping.")
                        continue

                    logger.info(f"===== Processing TaxID: {taxid} for Process ID: {process_id} =====")

                    # Define output fasta file path & check if it already exists
                    output_seq_path = os.path.join(output_directory, f"{process_id}.fasta")
                    if os.path.exists(output_seq_path):
                        logger.info(f"Sequence file '{output_seq_path}' already exists. Skipping.")
                        continue

                    # Get full taxonomy first for logging
                    full_taxonomy, taxon_ranks, initial_rank, taxon_ids = fetch_taxonomy_by_taxid(taxid)
                    if full_taxonomy:
                        logger.info(f"Full taxonomy for TaxID {taxid}: {', '.join(full_taxonomy)}")
                        
                        best_record, best_taxonomy, matched_rank = fetch_protein_seq_by_taxid(
                            taxid, 
                            gene_name, 
                            existing_taxonomy=full_taxonomy,
                            existing_rank=initial_rank,
                            existing_taxon_ranks=taxon_ranks,
                            existing_taxon_ids=taxon_ids
                        )

                        if best_record:
                            # Process successful result
                            original_accession = best_record.id
                            best_record.id = process_id
                            best_record.name = ""
                            best_record.description = ""

                            SeqIO.write(best_record, output_seq_path, "fasta")
                            logger.info(f"Written sequence to '{output_seq_path}'")

                            # Write to summary CSV
                            summary_writer.writerow([
                                process_id,
                                taxid,
                                original_accession,
                                len(best_record.seq),
                                matched_rank,
                                "; ".join(best_taxonomy) if best_taxonomy else "",
                                process_id,
                                os.path.abspath(output_seq_path)
                            ])
                            summary_file.flush()
                        else:
                            # Log failed search
                            with open(failed_searches_path, "a", newline='') as f:
                                writer = csv.writer(f)
                                writer.writerow([process_id, taxid, "No sequence found", time.strftime("%Y-%m-%d %H:%M:%S")])
                            logger.warning(f"No valid record found for TaxID {taxid}")
                    else:
                        # Log taxonomy fetch failure
                        with open(failed_searches_path, "a", newline='') as f:
                            writer = csv.writer(f)
                            writer.writerow([process_id, taxid, "Taxonomy fetch failed", time.strftime("%Y-%m-%d %H:%M:%S")])
                        logger.warning(f"Could not fetch taxonomy for TaxID {taxid}")

                except Exception as e:
                    logger.error(f"Error processing row for Process ID {process_id}: {e}")
                    with open(failed_searches_path, "a", newline='') as f:
                        writer = csv.writer(f)
                        writer.writerow([process_id, taxid, f"Unexpected error: {str(e)}", time.strftime("%Y-%m-%d %H:%M:%S")])
                    continue

if __name__ == "__main__":
    main()
