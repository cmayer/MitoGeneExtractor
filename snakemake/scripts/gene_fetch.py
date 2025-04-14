"""
Gene Fetch - NCBI Sequence Retrieval Tool

This script fetches gene sequences from NCBI databases based on taxonomy IDs (taxids).
It can retrieve both protein and nucleotide sequences, with support for various genes
including protein-coding genes (e.g., cox1, cox2, cytb, rbcl, matk) and rRNA genes (e.g., 16S, 18S).

Key Features:
- Taxonomic traversal: If sequences are not found at the input taxonomic level (e.g. species), 
  searches up higher taxonomic ranks (genus, family, etc.)
- Selects the longest high-quality sequence when multiple matches exist
- Customisable length filtering thresholds
- CDS extraction for protein-coding genes when '--type both' is specified
- Detailed logging of operations and progress tracking
- Efficient processing of repeated taxonomy queries through caching
- Optimised sequence retrieval by employing prefiltering step when >50 sequences found
- Single-taxid mode (-s/--single) for retrieving all available sequences for a specific taxon (-i not required)
- Rate limiting to comply with NCBI API guidelines (10 requests/second with API key)
- Robust error handling with automatic retries for NCBI API calls, including exponential backoff

Input:
- NCBI account email address and API key (see: https://support.nlm.nih.gov/kbArticle/?pn=KA-05317)
- CSV file containing taxonomy IDs (must have 'taxid' and 'ID' column)
- Gene name (e.g., 'cox1', '16s', 'rbcl', 'matk')
- Output directory path (will create new directories)
- Sequence type ('protein', 'nucleotide', or 'both')
- Optional: Minimum sequence length thresholds for filtering results
- Optional: Single-taxid mode

Output:
- FASTA files containing retrieved sequences (named by ID)
- Log file
- CSV file tracking successful retrievals with taxonomy, accession number, length, path to FASTA files
- CSV file documenting failed retrievals with error types
- Single-taxid mode: CSV file(s) containing accession number, length, and sequence description

Dependencies:
- Python>=3.9
- Biopython>=1.84
- ratelimit>=2.2.1

Usage:
    python gene_fetch.py -g/--gene <gene_name> -o/--out <output_directory> --type <sequence_type> 
                        [-i/--in <samples.csv>] [-s/--single <taxid>] 
                        [--protein_size <min_size>] [--nucleotide_size <min_size>]
                        [-e/--email <email>] [-k/--api-key <api_key>]
    
    # Required arguments:
    -e/--email             Email to use for NCBI API requests
    -k/--api-key           API key to use for NCBI API requests
    -g/--gene              Gene name to search for (e.g., cox1, 16s, rbcl)
    -o/--out               Directory to save output files
    --type                 Sequence type to fetch (protein, nucleotide, or both)
    -i/--in                Input CSV file with taxonomy IDs (required unless using --single)

    # Optional arguments:
    -s/--single            Single TaxID to fetch all available sequences for
    --protein_size         Minimum protein sequence length (default: 500)
    --nucleotide_size      Minimum nucleotide sequence length (default: 1500)

Examples:
    # Single taxid mode for retrieving all available sequences
    python gene_fetch.py -e your.email@domain.com -k your_api_key -g <gene_name> -o <output_directory> -s <taxid> --type <sequence_type>

    # Get all available rbcL sequences for a specific plant family (taxid: 3700)
    python gene_fetch.py -e your.email@domain.com -k your_api_key -g rbcl -o ./orchid_rbcl -s 3700 --type both

    # Standard usage retrieving both protein and nucleotide sequences for cox1
    python gene_fetch.py -e your.email@domain.com -k your_api_key -g cox1 -o ./output_dir -i ./samples.csv --type both --protein_size 500 --nucleotide_size 1500
    
    # Fetch only protein sequences for cytochrome b
    python gene_fetch.py -e your.email@domain.com -k your_api_key -g cytb -o ./cytb_proteins -i arthropods.csv --type protein --protein_size 300
    
    # Retrieve 16S rRNA nucleotide sequences
    python gene_fetch.py -e your.email@domain.com -k your_api_key -g 16s -o ./16s_sequences -i bacteria.csv --type nucleotide --nucleotide_size 1000
    
    # Using custom email and API key
    python gene_fetch.py -e your.email@domain.com -k your_api_key -g cox1 -o ./output_dir -i ./samples.csv --type both 
    
Notes:
- For protein-coding genes, --type both first fetches protein and then the corresponding nucleotide sequence
- Progress updates are logged every 10 samples by default
- NCBI API key is recommended and will increase rate limits from 3 to 10 requests/second
- The tool avoids "unverified" sequences by default
- CDS extraction includes fallback mechanisms for atypical annotation formats
- When more than 50 matching sequences are found for a sample, the tool employs an efficient two-step process:
  1. Fetches summary information for all matches to determine sequence lengths (using NCBI esummary API)
  2. Processes only the top 10 longest sequences, significantly reducing full record API calls

Single Mode Operation:
- When using `-s/--single` flag with a taxid, the tool switches to exhaustive retrieval mode
- In single mode, default length thresholds are reduced (protein: 100aa, nucleotide: 200bp)
- All matching sequences are retrieved
- Output files are named by their accession IDs
- Useful for creating comprehensive reference databases or exploring variation within a taxon

Future Development:
- Implement more sophisticated sequence quality filtering based on metadata
- Add optional alignment of retrieved sequences
- Add support for direct GenBank submission format output
- Enhance LRU caching for taxonomy lookups to reduce API calls
- Improve efficiency of record searching and selecting the longest sequence
- Add support for additional genetic markers beyond the currently supported set

Author: D. Parsons
Version: 1.0.4
License: MIT
"""

import csv
import sys
import os
import time
import argparse
import logging
from Bio import Entrez, SeqIO
from functools import wraps, lru_cache
from ratelimit import limits, sleep_and_retry
from requests.exceptions import RequestException
from time import sleep
from random import uniform
from urllib.error import HTTPError
from dataclasses import dataclass
from typing import Optional, Tuple, List, Dict, Any
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from Bio.Seq import translate, Seq
from Bio.SeqFeature import FeatureLocation
from http.client import IncompleteRead
import re



# Initialise logger at module level
logger = logging.getLogger("gene_fetch")

def log_progress(current: int, total: int, interval: int = 10) -> None:
    """Log progress at specified intervals."""
    if current == 0:
        logger.info(f"Starting processing: 0/{total} samples processed (0%)")
    elif current == total:
        logger.info(f"Completed processing: {total}/{total} samples processed (100%)")
    elif current % interval == 0:
        percentage = (current / total) * 100
        logger.info(f"Progress: {current}/{total} samples processed ({percentage:.2f}%)")


def setup_argument_parser():
    parser = argparse.ArgumentParser(description='Fetch gene sequences from NCBI databases.')
    
    parser.add_argument('-g', '--gene', required=True,
                      help='Name of gene to search for in NCBI RefSeq database (e.g., cox1)')
    
    parser.add_argument('-o', '--out', required=True,
                      help='Path to directory to save output files')
    
    parser.add_argument('-i', '--in', required=True, dest='input_csv',
                      help='Path to input CSV file containing TaxIDs')

    parser.add_argument('-s', '--single', type=str,
                      help='Single TaxID to fetch all available sequences for')
    
    parser.add_argument('--type', required=True, choices=['protein', 'nucleotide', 'both'],
                      help='Specify sequence type to fetch')
    
    parser.add_argument('--protein_size', type=int, default=500,
                      help='Minimum protein sequence length (default: 500)')
    
    parser.add_argument('--nucleotide_size', type=int, default=1500,
                      help='Minimum nucleotide sequence length (default: 1500)')
    
    parser.add_argument('-e', '--email', type=str, required=True,
                      help='Email to use for NCBI API requests (required)')
    
    parser.add_argument('-k', '--api-key', type=str, required=True,
                      help='API key to use for NCBI API requests (required)')
    
    return parser

@dataclass
class Config:
        def __init__(self, email, api_key):
            # Email and API key are now required
            if not email:
                raise ValueError("Email address is required for NCBI API requests. Use -e/--email to provide your email.")
            if not api_key:
                raise ValueError("API key is required for NCBI API requests. Use -k/--api-key to provide your API key.")
            
            self.email = email
            self.api_key = api_key

            # With API key, we can make up to 10 requests per second
            self.max_calls_per_second = 10

            # Default batch size for fetching sequences
            self.fetch_batch_size = 100

            # Delay between batches (seconds)
            self.batch_delay = (1, 2)  # uniform random delay between 1-2 seconds

            # Set search 'type'
            self.valid_sequence_types = frozenset({'protein', 'nucleotide', 'both'})

            # Minimum nucleotide and protein lengths for 'normal' mode
            self.protein_length_threshold = 500
            self.nucleotide_length_threshold = 1500
            
            # Minimum nucleotide and protein lengths for 'single' mode
            self.min_nucleotide_size_single_mode = 200
            self.min_protein_size_single_mode = 100
            
            self.gene_search_term = ""
            
            # Define gene type categories
            self._rRNA_genes = {
                '16s': [
                    '16S ribosomal RNA[Title]',
                    '16S rRNA[Title]',
                    '16S[Title]',
                ],
                '18s': [
                    '18S ribosomal RNA[Title]',
                    '18S rRNA[Title]',
                    '18S[Title]',
                ],
                '28s': [
                    '28S ribosomal RNA[Title]',
                    '28S rRNA[Title]',
                    '28S[Title]',
                ],
                '12s': [
                    '12S ribosomal RNA[Title]',
                    '12S rRNA[Title]',
                    '12S[Title]',
                ]
            }
            
            self._protein_coding_genes = {
                'cox1': [
                    'cox1[Gene]',
                    'COI[Gene]',
                    '"cytochrome c oxidase subunit 1"[Protein Name]',
                    '"cytochrome oxidase subunit 1"[Protein Name]',
                    '"cytochrome c oxidase subunit I"[Protein Name]',
                    '"COX1"[Protein Name]',
                    '"COXI"[Protein Name]'
                ],
                'cox2': [
                    'cox2[Gene]',
                    'COII[Gene]',
                    '"cytochrome c oxidase subunit 2"[Protein Name]',
                    '"cytochrome oxidase subunit 2"[Protein Name]',
                    '"cytochrome c oxidase subunit II"[Protein Name]',
                    '"COX2"[Protein Name]',
                    '"COXII"[Protein Name]'
                ],
                'cox3': [
                    'cox3[Gene]',
                    'COIII[Gene]',
                    '"cytochrome c oxidase subunit 3"[Protein Name]',
                    '"cytochrome oxidase subunit 3"[Protein Name]',
                    '"cytochrome c oxidase subunit III"[Protein Name]',
                    '"COX3"[Protein Name]',
                    '"COXIII"[Protein Name]'
                ],
                'cytb': [
                    'cytb[Gene]',
                    'cob[Gene]',
                    '"cytochrome b"[Protein Name]',
                    '"cytochrome b"[Gene]',
                    '"CYTB"[Protein Name]'
                ],
                'nd1': [
                    'nd1[Gene]',
                    'NAD1[Gene]',
                    '"NADH dehydrogenase subunit 1"[Protein Name]',
                    '"ND1"[Protein Name]'
                ],
                'rbcl': [
                    'rbcL[Gene]',
                    'RBCL[Gene]',
                    '"ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit"[Protein Name]',
                    '"ribulose 1,5-bisphosphate carboxylase/oxygenase large subunit"[Protein Name]',
                    '"ribulose bisphosphate carboxylase large chain"[Protein Name]',
                    '"RuBisCO large subunit"[Protein Name]',
                    '"ribulose-1,5-bisphosphate carboxylase/oxygenase small subunit"[Protein Name]',
                    '"ribulose 1,5-bisphosphate carboxylase/oxygenase small subunit"[Protein Name]',
                    '"ribulose-1,5-bisphosphate carboxylase/oxygenase small chain"[Protein Name]',
                    '"RuBisCO small subunit"[Protein Name]',
                    '"Ribulose bisphosphate carboxylase/oxygenase activase"[Gene]',
                    '"rbcL gene"[Gene]',
                    '"RBCL gene"[Gene]',
                    'rbcL[Title]',
                    'RBCL[Title]',
                    '"ribulose-1,5-bisphosphate carboxylase"[All Fields]',
                    '"ribulose 1,5-bisphosphate carboxylase"[All Fields]', 
                    '"ribulose bisphosphate carboxylase"[All Fields]',
                    'RuBisCO[All Fields]'
                ],
                'matk': [
                    'matK[Gene]',
                    'MATK[Gene]',
                    '"maturase K"[Protein Name]',
                    '"maturase K"[Gene]',
                    '"maturase-K"[Protein Name]',
                    '"maturase-K"[Gene]',
                    '"Maturase K"[Protein Name]',
                    '"Maturase K"[Gene]',
                    '"matK gene"[Gene]',
                    '"MATK gene"[Gene]',
                    '"trnK-matK"[Gene]',
                    '"maturase type II intron splicing factor"[Protein Name]',
                    '"chloroplast group II intron splicing factor maturase K"[Protein Name]',
                    '"type II intron maturase K"[Protein Name]',
                    '"tRNA-lysine maturase K"[Protein Name]'
                ]
            }

        def update_thresholds(self, protein_size: int, nucleotide_size: int):
            """Update sequence length thresholds."""
            self.protein_length_threshold = protein_size
            self.nucleotide_length_threshold = nucleotide_size

        def set_gene_search_term(self, gene_name: str):
            """Set search term based on gene name and type."""
            gene_name = gene_name.lower()
            
            # Check if it's an rRNA gene
            if gene_name in self._rRNA_genes:
                self.gene_search_term = '(' + ' OR '.join(self._rRNA_genes[gene_name]) + ')'
                search_type = "rRNA"
                
            # Check if it's a protein-coding gene
            elif gene_name in self._protein_coding_genes:
                self.gene_search_term = '(' + ' OR '.join(self._protein_coding_genes[gene_name]) + ')'
                search_type = "protein-coding"
                
            else:
                # Generic search term for unknown genes
                self.gene_search_term = (
                    f'({gene_name}[Title] OR {gene_name}[Gene] OR "{gene_name}"[Text Word])'
                )
                search_type = "generic"
                
            return search_type  # Return the type for logging purposes



def ensure_directory(path: Path) -> None:
    """Ensure directory exists, create if it doesn't."""
    path.mkdir(parents=True, exist_ok=True)






def setup_logging(output_dir: Path) -> logging.Logger:
    """Setup logging with both file and console handlers."""
    ensure_directory(output_dir)
    
    # Clear existing handlers
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
    
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    log_path = output_dir / "gene_fetch.log"
    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    logger.info(f"Logging initialised. Log file: {log_path}")
    return logger




def check_ncbi_status():
    """Check NCBI service status using einfo endpoint."""
    try:
        handle = Entrez.einfo()
        result = Entrez.read(handle)
        handle.close()
        # If we can successfully query einfo, services are up
        return True
    except Exception as e:
        logger.warning(f"NCBI service check failed: {str(e)}")
        return False




def enhanced_retry(exceptions: tuple, tries: int = 4, initial_delay: int = 10, 
                  backoff: int = 2, max_delay: int = 240):
    """Enhanced retry decorator with NCBI status checking and longer delays."""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            mtries, mdelay = tries, initial_delay
            for i in range(tries):
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    if i == tries - 1:
                        logger.error(f"Final attempt failed: {str(e)}")
                        return None
                    
                    # Check if it's a server error (500)
                    if isinstance(e, HTTPError) and e.code == 500:
                        # Perform NCBI service check
                        service_status = check_ncbi_status()
                        if not service_status:
                            logger.warning("NCBI services appear to be experiencing issues")
                            # Use a longer delay for service issues
                            mdelay = min(mdelay * 4, max_delay)
                        
                    # Adjust delay based on error type
                    if isinstance(e, HTTPError):
                        if e.code == 429:  # Too Many Requests
                            mdelay = min(mdelay * 3, max_delay)
                        elif e.code >= 500:  # Server errors
                            mdelay = min(mdelay * 4, max_delay)
                    
                    # Add jitter to avoid thundering herd
                    delay_with_jitter = mdelay + uniform(-0.1 * mdelay, 0.1 * mdelay)
                    logger.warning(f"{str(e)}, Retrying in {delay_with_jitter:.2f} seconds...")
                    sleep(delay_with_jitter)
                    
                    # Progressive backoff
                    mdelay = min(mdelay * backoff, max_delay)
                    
                    # Additional delay for IncompleteRead errors
                    if isinstance(e, IncompleteRead):
                        logger.warning("Incomplete read detected, adding additional delay")
                        sleep(uniform(5, 10))  # Add 5-10 seconds extra delay
            
            return None
        return wrapper
    return decorator





class EntrezHandler:
    """Handles all Entrez API interactions with enhanced error handling."""
    def __init__(self, config: Config):
        self.config = config
        Entrez.email = config.email
        Entrez.api_key = config.api_key
        self.consecutive_errors = 0
        self.last_service_check = 0
        self.service_check_interval = 60  # Seconds between service checks
        
    def should_check_service_status(self) -> bool:
        """Determine if we should perform a service status check."""
        current_time = time.time()
        if current_time - self.last_service_check > self.service_check_interval:
            return True
        return False
        
    def handle_request_error(self, error: Exception) -> None:
        """Handle request errors and track consecutive failures."""
        self.consecutive_errors += 1
        if self.consecutive_errors >= 3 and self.should_check_service_status():
            service_status = check_ncbi_status()
            self.last_service_check = time.time()
            if not service_status:
                logger.warning("Multiple consecutive errors and NCBI service appears to be down")
                sleep(uniform(30, 60))  # Longer delay when service is down
        
    def handle_request_success(self) -> None:
        """Reset error counter on successful request."""
        self.consecutive_errors = 0
        
    @enhanced_retry((HTTPError, RuntimeError, IOError, IncompleteRead))
    @sleep_and_retry
    @limits(calls=10, period=1.1)
    def fetch(self, **kwargs) -> Optional[Any]:
        """Execute an Entrez efetch query with enhanced error handling."""
        try:
            result = Entrez.efetch(**kwargs)
            self.handle_request_success()
            return result
        except Exception as e:
            self.handle_request_error(e)
            raise

    def search(self, **kwargs) -> Optional[Dict]:
        """Execute an Entrez esearch query with enhanced error handling and pagination."""
        @enhanced_retry((HTTPError, RuntimeError, IOError, IncompleteRead))
        @sleep_and_retry
        @limits(calls=10, period=1.1)
        def _do_search(**kwargs):
            try:
                handle = Entrez.esearch(**kwargs)
                result = Entrez.read(handle)
                handle.close()
                self.handle_request_success()
                return result
            except Exception as e:
                self.handle_request_error(e)
                raise

        # First search to get total count
        initial_result = _do_search(**kwargs)
        if not initial_result:
            return None
            
        total_count = int(initial_result['Count'])
        all_ids = initial_result['IdList']
        
        # If there are more results, fetch them
        batch_size = 100  # Increased from default 20
        if total_count > len(all_ids):
            for start in range(len(all_ids), total_count, batch_size):
                kwargs['retstart'] = start
                kwargs['retmax'] = batch_size
                result = _do_search(**kwargs)
                if result and result.get('IdList'):
                    all_ids.extend(result['IdList'])
                
                # Add delay between batches
                sleep(uniform(1, 2))
        
        # Return modified result with all IDs
        initial_result['IdList'] = all_ids
        return initial_result




class SequenceProcessor:
    """Handles sequence processing and validation with WGS support."""
    
    def __init__(self, config: Config, entrez: EntrezHandler):
        self.config = config
        self.entrez = entrez
        
    def parse_contig_line(self, contig_line: str) -> Optional[Tuple[str, int, int]]:
        """
        Parse a GenBank CONTIG line to extract WGS contig information.
        
        Args:
            contig_line: The CONTIG line from a GenBank record
            
        Returns:
            Optional[Tuple[str, int, int]]: (contig_id, start, end) if found
        """
        try:
            # Remove 'join(' and ')' if present
            cleaned = contig_line.strip().replace('join(', '').replace(')', '')
            
            # Parse the contig reference
            # Format is typically: WVEN01000006.2:1..16118
            if ':' in cleaned:
                contig_id, coords = cleaned.split(':')
                if '..' in coords:
                    start, end = coords.split('..')
                    return contig_id, int(start), int(end)
        except Exception as e:
            logger.error(f"Error parsing CONTIG line '{contig_line}': {e}")
        return None
        
    def fetch_wgs_sequence(self, record: SeqRecord) -> Optional[SeqRecord]:
        try:
            # Check for CONTIG line in annotations
            contig_line = record.annotations.get('contig', None)
            if not contig_line:
                logger.warning(f"No CONTIG line found in WGS record {record.id}")
                return None
                
            # Parse the CONTIG line
            contig_info = self.parse_contig_line(contig_line)
            if not contig_info:
                logger.error(f"Could not parse CONTIG line: {contig_line}")
                return None
                
            contig_id, start, end = contig_info
            logger.info(f"Found WGS contig reference: {contig_id} positions {start}..{end}")
            
            # Fetch the actual contig sequence
            try:
                handle = self.entrez.fetch(
                    db="nucleotide",
                    id=contig_id,
                    rettype="fasta",
                    retmode="text"
                )
                contig_record = next(SeqIO.parse(handle, "fasta"))
                handle.close()
                
                if not contig_record or not contig_record.seq:
                    logger.error(f"Failed to fetch sequence for contig {contig_id}")
                    return None
                    
                # Extract the relevant portion
                start_idx = start - 1  # Convert to 0-based indexing
                sequence = contig_record.seq[start_idx:end]
                
                if not sequence:
                    logger.error(f"Extracted sequence is empty for positions {start}..{end}")
                    return None
                    
                # Create new record with the sequence
                new_record = record[:]
                new_record.seq = sequence
                logger.info(f"Successfully extracted {len(sequence)} bp from WGS contig")
                
                return new_record
                
            except Exception as e:
                logger.error(f"Error fetching WGS contig {contig_id}: {e}")
                return None
                
        except Exception as e:
            logger.error(f"Error processing WGS record: {e}")
            return None

    def fetch_nucleotide_record(self, record_id: str) -> Optional[SeqRecord]:
        """
        Fetch nucleotide sequence, including WGS records if they have an available fasta.
        """
        try:
            # Fetch the record
            handle = self.entrez.fetch(db="nucleotide", id=record_id, rettype="gb", retmode="text")
            record = next(SeqIO.parse(handle, "genbank"))
            handle.close()
            
            # Check if it's a WGS record
            is_wgs = False
            if hasattr(record, 'annotations'):
                keywords = record.annotations.get('keywords', [])
                if 'WGS' in keywords:
                    is_wgs = True
                    logger.info(f"WGS record detected for {record_id}")
            
            # For WGS records, check if there's a complete sequence available directly
            if is_wgs:
                if record.seq is not None and len(record.seq) > 0:
                    try:
                        # Verify we can access the sequence
                        seq_str = str(record.seq)
                        if seq_str and not seq_str.startswith('?'):
                            logger.info(f"WGS record {record_id} has a complete sequence of length {len(record.seq)}")
                            return record
                        else:
                            logger.info(f"WGS record {record_id} has a placeholder/incomplete sequence")
                    except Exception as e:
                        logger.error(f"Error accessing sequence content for WGS record {record_id}: {e}")
                
                # If we don't have a complete sequence, check for CONTIG line
                if 'contig' in record.annotations:
                    logger.info(f"WGS record {record_id} has a CONTIG line, attempting to fetch underlying sequence")
                    wgs_record = self.fetch_wgs_sequence(record)
                    if wgs_record:
                        return wgs_record
                    else:
                        logger.warning(f"Failed to fetch sequence from WGS CONTIG for {record_id}")
                
                # If we still don't have a sequence, log and return None
                logger.info(f"WGS record {record_id} does not have a usable sequence - skipping")
                return None
            
            # Skip unverified sequences
            if 'unverified' in record.description.lower() or 'UNVERIFIED' in record.description:
                logger.info(f"Unverified sequence detected for {record_id} - skipping")
                return None
                
            # For non-WGS and verified records, verify sequence content
            if record.seq is not None and len(record.seq) > 0:
                try:
                    # Verify we can access the sequence
                    _ = str(record.seq)
                    return record
                except Exception as e:
                    logger.error(f"Undefined sequence content for {record_id}: {e}")
                    return None
                    
            return None
            
        except Exception as e:
            logger.error(f"Error fetching nucleotide sequence for {record_id}: {e}")
            return None

    def extract_nucleotide(self, record: SeqRecord, gene_name: str, single_mode: bool = False) -> Optional[SeqRecord]:
            """
            Extract CDS region from a sequence record with robust fallbacks for various sequence types.
            """
            logger.info(f"Attempting to extract CDS for gene {gene_name} (accession {record.id})")
            
            # Prepare gene name variations for matching
            gene_variations = set()
            pattern_variations = []

            # Get minimum size threshold for single mode
            min_size = self.config.min_nucleotide_size_single_mode if single_mode else 100
            
            if gene_name.lower() in self.config._protein_coding_genes:
                # Get the gene variations from config
                variations = self.config._protein_coding_genes[gene_name.lower()]
                gene_variations = {v.split('[')[0].strip('"').lower() for v in variations}
                
                # Add common pattern variations for different writing styles
                base_gene = gene_name.lower()
                if base_gene == 'rbcl':
                    pattern_variations = [
                        'rbcl', 'rbc-l', 'rbc l', 'rubisco', 'ribulose-1,5-bisphosphate', 'ribulose bisphosphate'
                    ]
                elif base_gene.startswith('cox'):
                    number = base_gene[-1]
                    pattern_variations = [
                        f'cytochrome c oxidase'
                    ]
                elif base_gene == 'cytb':
                    pattern_variations = ['cytb', 'cyt b', 'cyt-b', 'cytochrome b', 'cytochrome-b']
                elif base_gene == 'matk':
                    pattern_variations = [
                        'matk', 'mat-k', 'mat k', 'maturase k', 'maturase-k', 'maturase', 
                        'chloroplast maturase k', 'trnk-matk'
                    ]
                elif base_gene == 'nd1':
                    pattern_variations = [
                        'nd1', 'nd-1', 'nd 1', 'nadh1', 'nadh-1', 'nadh 1',
                        'nadh dehydrogenase 1', 'nadh dehydrogenase subunit 1',
                        'nadh-dehydrogenase 1', 'nad1', 'nad-1'
                    ]
                elif base_gene == 'nd2':
                    pattern_variations = [
                        'nd2', 'nd-2', 'nd 2', 'nadh2', 'nadh-2', 'nadh 2',
                        'nadh dehydrogenase 2', 'nadh dehydrogenase subunit 2',
                        'nadh-dehydrogenase 2', 'nad2', 'nad-2'
                    ]
                elif base_gene == 'nd4':
                    pattern_variations = [
                        'nd4', 'nd-4', 'nd 4', 'nadh4', 'nadh-4', 'nadh 4',
                        'nadh dehydrogenase 4', 'nadh dehydrogenase subunit 4',
                        'nadh-dehydrogenase 4', 'nad4', 'nad-4'
                    ]
                elif base_gene == 'nd5':
                    pattern_variations = [
                        'nd5', 'nd-5', 'nd 5', 'nadh5', 'nadh-5', 'nadh 5',
                        'nadh dehydrogenase 5', 'nadh dehydrogenase subunit 5',
                        'nadh-dehydrogenase 5', 'nad5', 'nad-5'
                    ]
                elif base_gene == 'atp6':
                    pattern_variations = [
                        'atp6', 'atp-6', 'atp 6', 'atpase6', 'atpase-6', 'atpase 6',
                        'atp synthase 6', 'atp synthase subunit 6', 'atp synthase f0 subunit 6',
                        'atpase subunit 6', 'atpase subunit a'
                    ]
                elif base_gene == 'atp8':
                    pattern_variations = [
                        'atp8', 'atp-8', 'atp 8', 'atpase8', 'atpase-8', 'atpase 8',
                        'atp synthase 8', 'atp synthase subunit 8', 'atp synthase f0 subunit 8',
                        'atpase subunit 8'
                    ]
            elif base_gene == '16s' or base_gene == '16s rrna' or base_gene == 'rrn16':
                pattern_variations = [
                    '16s', '16s rrna', '16s ribosomal rna', '16s ribosomal', 
                    '16 s rrna', '16 s', 'rrn16', 'rrn 16', 
                    'small subunit ribosomal rna', 'ssu rrna', 'ssu'
                ]
            elif base_gene == '18s' or base_gene == '18s rrna' or base_gene == 'rrn18':
                pattern_variations = [
                    '18s', '18s rrna', '18s ribosomal rna', '18s ribosomal', 
                    '18 s rrna', '18 s', 'rrn18', 'rrn 18', 
                    'small subunit ribosomal rna', 'ssu rrna', 'ssu'
                ]
            elif base_gene == '28s' or base_gene == '28s rrna' or base_gene == 'rrn28':
                pattern_variations = [
                    '28s', '28s rrna', '28s ribosomal rna', '28s ribosomal', 
                    '28 s rrna', '28 s', 'rrn28', 'rrn 28', 
                    'large subunit ribosomal rna', 'lsu rrna', 'lsu'
                ]
            elif base_gene == 'its' or base_gene == 'its1' or base_gene == 'its2':
                pattern_variations = [
                    'internal transcribed spacer', 'its region', 'its1-5.8s-its2', 
                    'its 1', 'its 2', 'its1', 'its2', 'its 1-5.8s-its 2',
                    'ribosomal its', 'rrna its'
                ]
            elif base_gene == 'trnh-psba' or base_gene == 'psba-trnh':
                pattern_variations = [
                    'trnh-psba', 'psba-trnh', 'trnh psba', 'psba trnh',
                    'trnh-psba spacer', 'psba-trnh spacer', 'trnh-psba intergenic spacer',
                    'trnh psba intergenic', 'psba trnh intergenic'
                ]
            else:
                logger.warning(f"No defined variations for gene {gene_name}")
                gene_variations = {gene_name.lower()}
                
            # If we have pattern variations, add them to regular variations
            if pattern_variations:
                gene_variations.update(pattern_variations)
                
            logger.info(f"Using gene variations for matching: {gene_variations}")
            
            # STEP 1: Try to find a CDS feature with EXACT match to target gene
            found_cds = None
            
            # First, look for exact matches to our target gene name
            target_gene = gene_name.lower()
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                    
                qualifiers = []
                for field in ['gene', 'product', 'note']:
                    qualifiers.extend(feature.qualifiers.get(field, []))
                    
                # Log qualifiers for debugging
                logger.info(f"Found CDS qualifiers: {qualifiers}")
                
                # Check for exact match first 
                for qualifier in qualifiers:
                    qualifier_lower = qualifier.lower()
                    
                    # Exact match to target gene (e.g., 'cox1', 'coi')
                    if target_gene in qualifier_lower.split() or f"{target_gene}" == qualifier_lower:
                        logger.info(f"Found exact match for {target_gene} in qualifier: {qualifier}")
                        found_cds = feature
                        break
                        
                    # For cox genes, check for the specific number match
                    if target_gene.startswith('cox'):
                        if "cox" in qualifier_lower and target_gene[-1] in qualifier_lower:
                            if f"cox{target_gene[-1]}" in qualifier_lower or f"cox {target_gene[-1]}" in qualifier_lower:
                                logger.info(f"Found cox{target_gene[-1]} match in qualifier: {qualifier}")
                                found_cds = feature
                                break
                    
                    # For coi/cox1 specific matching
                    if target_gene == 'cox1':
                        if 'coi' in qualifier_lower.split() or 'co1' in qualifier_lower.split():
                            logger.info(f"Found COI/CO1 match for cox1 in qualifier: {qualifier}")
                            found_cds = feature
                            break
                
                # If we found an exact match, break the loop
                if found_cds:
                    break
                    
            # If no exact match found, try the more general matching for any variant
            if not found_cds:
                for feature in record.features:
                    if feature.type != "CDS":
                        continue
                        
                    qualifiers = []
                    for field in ['gene', 'product', 'note']:
                        qualifiers.extend(feature.qualifiers.get(field, []))
                        
                    # Check if any variation matches in the qualifiers
                    if any(var in qualifier.lower() for qualifier in qualifiers 
                          for var in gene_variations):
                        logger.info(f"Found match using variations for {gene_name}")
                        found_cds = feature
                        break
            
            # If we found a matching CDS, extract it  
            if found_cds:
                try:
                    cds_record = record[:]
                    cds_record.seq = found_cds.extract(record.seq)
                    
                    if len(cds_record.seq) >= 100:  # Sanity check for minimum length
                        logger.info(f"Successfully extracted CDS of length {len(cds_record.seq)} (accession {record.id})")
                        return cds_record
                    else:
                        logger.warning(f"Extracted CDS too short ({len(cds_record.seq)} bp) (accession {record.id})")
                except Exception as e:
                    logger.error(f"CDS extraction error for {record.id}: {e}")
                    logger.error("Full error details:", exc_info=True)
            
            # If we're not in single mode, don't use fallbacks
            if not single_mode:
                logger.debug(f"No valid CDS found for gene {gene_name} (accession {record.id})")
                return None
            
            logger.info(f"No CDS feature found, trying fallbacks for single mode (accession {record.id})")
            
            # Define reasonable size limits for different genes - only used in single mode
            max_gene_sizes = {
                'rbcl': 2000,   # Typical rbcL is ~1400bp
                'cox1': 2000,   # Typical cox1 is ~1500bp
                'cox2': 2000,   # Typical cox2 is ~1500bp
                'cox3': 2000,   # Typical cox3 is ~1500bp
                'cytb': 1800,   # Typical cytb is ~1100bp
                'nd1': 1800,    # Typical nd1 is ~1000bp
                'nd2': 1800,    # Typical nd2 is ~1000bp
                'nd4': 1800,    # Typical nd4 is ~1300bp
                'nd5': 2000,    # Typical nd5 is ~1700bp
                'matk': 2000,   # Typical matK is ~1500bp
                'atp6': 1200,   # Typical atp6 is ~800bp
                'atp8': 600,    # Typical atp8 is ~400bp
                '16s': 2000,    # Typical 16S is ~1600bp
                '18s': 2500,    # Typical 18S is ~1800bp
                '28s': 3500,    # Typical 28S can be ~3000bp
                '12s': 1500,    # Typical 12S is ~1000bp
                'its': 1000,    # Typical ITS region is ~700bp
                'its1': 500,    # Typical ITS1 is ~300bp
                'its2': 500,    # Typical ITS2 is ~350bp
                'trnh-psba': 1000  # Typical trnH-psbA is ~500-700bp
            }
            
            # Get maximum acceptable size for this gene
            max_size = max_gene_sizes.get(gene_name.lower(), 3000)  # Default to 3000 for unknown genes
            
        # FALLBACK 1: Check for gene feature with matching name but no CDS
        # Move this up in priority as it's more specific than mRNA/EST check
            for feature in record.features:
                if feature.type == "gene":
                    gene_qualifiers = feature.qualifiers.get('gene', [])
                    gene_notes = feature.qualifiers.get('note', [])
                    all_qualifiers = gene_qualifiers + gene_notes
                    
                    # Check for exact match first
                    if target_gene in [q.lower() for q in all_qualifiers]:
                        logger.info(f"Found exact gene match: {target_gene}")
                        try:
                            gene_record = record[:]
                            gene_record.seq = feature.extract(record.seq)
                            
                            if len(gene_record.seq) > max_size:
                                logger.warning(f"Extracted gene region too large ({len(gene_record.seq)} bp > {max_size} bp limit) - skipping")
                                continue
                            
                            if len(gene_record.seq) >= min_size:
                                logger.info(f"Successfully extracted gene region of length {len(gene_record.seq)}")
                                return gene_record
                            else:
                                logger.warning(f"Extracted gene region too short ({len(gene_record.seq)} bp < {min_size} bp)")
                        except Exception as e:
                            logger.error(f"Gene region extraction error: {e}")
                    
                    # More general matching
                    qualifier_text = " ".join(all_qualifiers).lower()
                    
                    # Check if any variation matches in the qualifiers
                    if any(var in qualifier_text for var in gene_variations):
                        try:
                            logger.info(f"Found matching gene feature, using gene region (accession {record.id})")
                            gene_record = record[:]
                            gene_record.seq = feature.extract(record.seq)
                            
                            # Check if the extracted sequence is too large
                            if len(gene_record.seq) > max_size:
                                logger.warning(f"Extracted gene region too large ({len(gene_record.seq)} bp > {max_size} bp limit) - skipping (accession {record.id})")
                                continue
                            
                            if len(gene_record.seq) >= min_size:  # Use min_size instead of hardcoded 100
                                logger.info(f"Successfully extracted gene region of length {len(gene_record.seq)} (accession {record.id})")
                                return gene_record
                            else:
                                logger.warning(f"Extracted gene region too short ({len(gene_record.seq)} bp < {min_size} bp) (accession {record.id})")
                        except Exception as e:
                            logger.error(f"Gene region extraction error for {record.id}: {e}")
        
        # FALLBACK 2: Check if this is an mRNA sequence with no CDS feature
            mol_type = ""
            for feature in record.features:
                if feature.type == "source" and "mol_type" in feature.qualifiers:
                    mol_type = feature.qualifiers["mol_type"][0].lower()
                    break
        
            if mol_type in ["mrna", "est"]:
                logger.info(f"Sequence is mRNA/EST, checking if description matches target gene (accession {record.id})")
                
                description_lower = record.description.lower()
                
                # Check if any of our variations appear in the description
                matching_vars = [var for var in gene_variations if var in description_lower]
                if matching_vars:
                    logger.info(f"mRNA/EST description matches gene variations: {matching_vars}")
                    
                    # Check if the sequence is too large
                    if len(record.seq) > max_size:
                        logger.warning(f"mRNA/EST sequence too large ({len(record.seq)} bp > {max_size} bp limit) - skipping (accession {record.id})")
                        return None
                        
                    if len(record.seq) >= min_size:  # Use min_size instead of hardcoded 100
                        logger.info(f"Using complete mRNA/EST of length {len(record.seq)} as it matches gene variations (accession {record.id})")
                        return record
                    else:
                        logger.warning(f"mRNA/EST sequence too short ({len(record.seq)} bp < {min_size} bp) (accession {record.id})")
                else:
                    logger.info(f"Description doesn't match any target gene variation")
        
        # FALLBACK 3: If it's a partial sequence entry, check if gene name appears in description
            description_lower = record.description.lower()
            if "partial" in description_lower:
                matching_vars = [var for var in gene_variations if var in description_lower]
                if matching_vars:
                    logger.info(f"Found partial sequence matching gene variations: {matching_vars}")
                    
                    # Check if the sequence is too large
                    if len(record.seq) > max_size:
                        logger.warning(f"Partial sequence too large ({len(record.seq)} bp > {max_size} bp limit) - skipping (accession {record.id})")
                        return None
                        
                    if len(record.seq) >= min_size:  # Use min_size instead of hardcoded 100
                        logger.info(f"Using entire partial sequence of length {len(record.seq)} (accession {record.id})")
                        return record
                    else:
                        logger.warning(f"Partial sequence too short ({len(record.seq)} bp < {min_size} bp) (accession {record.id})")
        
        # FALLBACK 4: For all records, check if the target gene is in the organism name or sequence ID
        # This is a last resort when we're in single mode and desperate for more sequences
            org_name = ""
            for feature in record.features:
                if feature.type == "source" and "organism" in feature.qualifiers:
                    org_name = feature.qualifiers["organism"][0].lower()
                    break
        
            if gene_name.lower() in record.id.lower() or gene_name.lower() in org_name:
                logger.info(f"Last resort: Gene name {gene_name} found in sequence ID or organism name (accession {record.id})")
                
                # Check if the sequence is too large
                if len(record.seq) > max_size:
                    logger.warning(f"Sequence too large ({len(record.seq)} bp > {max_size} bp limit) - skipping as last resort (accession {record.id})")
                    return None
                    
                if len(record.seq) >= min_size:  # Use min_size instead of hardcoded 100
                    logger.info(f"Using entire sequence of length {len(record.seq)} as a last resort (accession {record.id})")
                    return record
                else:
                    logger.warning(f"Sequence too short ({len(record.seq)} bp < {min_size} bp) - skipping as last resort (accession {record.id})")
        
            logger.debug(f"No valid CDS or fallback found for gene {gene_name} (accession {record.id})")
            return None

    def parse_coded_by(self, coded_by: str) -> Tuple[Optional[List[Tuple[str, Optional[Tuple[int, int]]]]], bool]:
        """
        Parse complex coded_by expressions including complement() and join() statements.
        Returns a tuple of (segments, is_complement) where segments is a list of (accession, coordinates) 
        tuples to handle split sequences, and is_complement indicates if sequence should be reverse complemented.
        """
        logger.info(f"Parsing coded_by qualifier: {coded_by}")
        try:
            # Determine if complement first
            is_complement = coded_by.startswith('complement(')
            
            # Remove outer wrapper (complement or join)
            if is_complement:
                coded_by = coded_by[10:-1]  # Remove 'complement(' and final ')'
            elif coded_by.startswith('join('):
                coded_by = coded_by[5:-1]   # Remove 'join(' and final ')'
                
            # Split by comma while preserving full coordinates
            segments_raw = []
            current_segment = ""
            in_parentheses = 0
            
            # Careful splitting to preserve full coordinates
            for char in coded_by:
                if char == '(':
                    in_parentheses += 1
                elif char == ')':
                    in_parentheses -= 1
                elif char == ',' and in_parentheses == 0:  # Only split at top level
                    segments_raw.append(current_segment.strip())
                    current_segment = ""
                    continue
                current_segment += char
            segments_raw.append(current_segment.strip())
            
            # Clean up the segments
            cleaned_segments = []
            for seg in segments_raw:
                seg = seg.strip().strip('()')
                cleaned_segments.append(seg)
                
            logger.debug(f"Cleaned segments: {cleaned_segments}")
            
            # Process all segments
            result = []
            all_coordinates_valid = True
            
            for segment in cleaned_segments:
                if not segment:  # Skip empty segments
                    continue
                    
                logger.debug(f"Processing segment: '{segment}'")
                
                # Extract accession and coordinates
                if ':' in segment:
                    accession, coords = segment.split(':')
                    accession = accession.strip()
                    coords = coords.strip()
                    logger.debug(f"Split into accession: '{accession}', coords: '{coords}'")
                    
                    if '..' in coords:
                        coord_parts = coords.split('..')
                        if len(coord_parts) != 2:
                            logger.error(f"Invalid coordinate format: {coords}")
                            all_coordinates_valid = False
                            break
                            
                        start_str, end_str = coord_parts
                        # Remove any non-digit characters
                        start_str = ''.join(c for c in start_str if c.isdigit())
                        end_str = ''.join(c for c in end_str if c.isdigit())
                        
                        logger.debug(f"Cleaned coordinate strings - Start: '{start_str}', End: '{end_str}'")
                        
                        try:
                            start = int(start_str)
                            end = int(end_str)
                            logger.debug(f"Parsed coordinates - Start: {start}, End: {end}")
                            
                            if start <= 0 or end <= 0:
                                logger.error(f"Invalid coordinates: must be positive numbers")
                                all_coordinates_valid = False
                                break
                                
                            if start > end:
                                logger.error(f"Invalid coordinate range: {start}..{end}")
                                all_coordinates_valid = False
                                break
                                
                            result.append((accession, (start, end)))
                            logger.info(f"Successfully parsed coordinates: {start}-{end} for {accession}")
                            
                        except ValueError as ve:
                            logger.error(f"Failed to parse coordinates '{coords}': {ve}")
                            all_coordinates_valid = False
                            break
                    else:
                        logger.error(f"Missing coordinate separator '..' in {coords}")
                        all_coordinates_valid = False
                        break
                else:
                    logger.error(f"Missing accession separator ':' in {segment}")
                    all_coordinates_valid = False
                    break
            
            if not all_coordinates_valid or not result:
                logger.error("Failed to parse one or more segments")
                return None, False
            
            logger.info(f"Successfully parsed {len(result)} segments")
            return result, is_complement
                
        except Exception as e:
            logger.error(f"Error parsing coded_by: {coded_by}, error: {e}")
            logger.error("Full error details:", exc_info=True)
            return None, False
    
    def fetch_nucleotide_from_protein(self, protein_record: SeqRecord, gene_name: str) -> Optional[SeqRecord]:
       """
       Fetch nucleotide sequence corresponding to a protein record.
       Handles both RefSeq coded_by qualifiers and UniProt xrefs.
       """
       try:
           logger.info(f"Attempting to fetch nucleotide sequence for protein record {protein_record.id}")
           
           # Try coded_by qualifier for RefSeq records
           cds_feature = next((f for f in protein_record.features if f.type == "CDS"), None)
           if cds_feature and 'coded_by' in cds_feature.qualifiers:
               coded_by = cds_feature.qualifiers['coded_by'][0]
               
               parsed_result = self.parse_coded_by(coded_by)
               if parsed_result:
                   segments, is_complement = parsed_result
                   
                   if not segments:
                       logger.error(f"No valid segments found in coded_by qualifier for {protein_record.id}")
                       return None
                   
                   # Fetch full sequence for first accession
                   first_accession = segments[0][0]
                   logger.info(f"Fetching nucleotide sequence for accession: {first_accession}")
                   
                   # Use enhanced nucleotide fetching
                   nucleotide_record = self.fetch_nucleotide_record(first_accession)
                   if not nucleotide_record:
                       logger.error(f"Failed to fetch nucleotide sequence for accession: {first_accession}")
                       return None
                   
                   # Extract and join all segments
                   complete_sequence = ""
                   for accession, coordinates in segments:
                       if accession != first_accession:
                           nucleotide_record = self.fetch_nucleotide_record(accession)
                           if not nucleotide_record:
                               logger.error(f"Failed to fetch additional sequence: {accession}")
                               continue
                       
                       if coordinates:
                           start, end = coordinates
                           segment_seq = str(nucleotide_record.seq[start-1:end])
                           if len(segment_seq) == 0:
                               logger.error(f"Zero-length sequence extracted using coordinates {start}..{end} (accession {accession})")
                               return None
                       else:
                           segment_seq = str(nucleotide_record.seq)
                       
                       complete_sequence += segment_seq
                   
                   # Handle complement if needed
                   if is_complement:
                       complete_sequence = str(Seq(complete_sequence).reverse_complement())
                   
                   # Create new record with complete sequence
                   new_record = nucleotide_record[:]
                   new_record.seq = Seq(complete_sequence)
                   logger.info(f"Successfully extracted nucleotide sequence of length {len(complete_sequence)} (from protein {protein_record.id})")
                   
                   return new_record
                   
           logger.warning(f"No valid nucleotide reference found in protein record {protein_record.id}")
           return None
           
       except Exception as e:
           logger.error(f"Error in fetch_nucleotide_from_protein for {protein_record.id}: {e}")
           logger.error("Full error details:", exc_info=True)
           return None

    @enhanced_retry((HTTPError, RuntimeError, IOError, IncompleteRead), tries=5, initial_delay=15)
    def fetch_taxonomy(self, taxid: str) -> Tuple[List[str], Dict[str, str], str, Dict[str, str]]:
        """
        Fetch taxonomy information for a given TaxID.
        """
        logger.info(f"Fetching taxonomy for TaxID: {taxid}")
        
        # First verify taxid format
        taxid = taxid.strip()
        if not taxid.isdigit():
            logger.error(f"Invalid TaxID format: {taxid} (must be numerical)")
            return [], {}, "", {}

        # Add initial delay to help avoid rate limiting
        sleep(uniform(0.5, 1.0))

        try:
            # Set up parameters for the request
            params = {
                "db": "taxonomy",
                "id": taxid,
                "email": self.config.email,
                "api_key": self.config.api_key,
                "tool": "gene_fetch"
            }

            max_retries = 3
            for attempt in range(max_retries):
                try:
                    handle = self.entrez.fetch(**params)
                    records = Entrez.read(handle)
                    handle.close()
                    break  # If successful, exit retry loop
                except HTTPError as e:
                    if e.code == 400:
                        if attempt < max_retries - 1:  # If not the last attempt
                            delay = (attempt + 1) * 2  # Progressive delay
                            logger.warning(f"HTTP 400 error for TaxID {taxid}, attempt {attempt + 1}/{max_retries}. Retrying in {delay} seconds...")
                            sleep(delay)
                            continue
                        else:
                            logger.error(f"Failed to fetch taxonomy after {max_retries} attempts for TaxID {taxid}")
                            return [], {}, "", {}
                    else:
                        raise  # Re-raise other HTTP errors
            else:  # If we exhaust all retries
                return [], {}, "", {}

            if not records:
                logger.error(f"No taxonomy records found for TaxID {taxid}")
                return [], {}, "", {}

            # Get the first record
            record = records[0]

            # Initialise rank information dictionaries
            rank_info = {}
            taxid_info = {}

            lineage_nodes = record.get("LineageEx", [])
            lineage = []
            
            # Process lineage nodes
            for node in lineage_nodes:
                name = node.get("ScientificName", "")
                rank = node.get("Rank", "no rank")
                node_taxid = str(node.get("TaxId", ""))
                
                # Add to lineage
                lineage.append(name)
                
                # Add rank and taxid info if valid
                if name and rank != "no rank":
                    rank_info[name] = rank
                    taxid_info[name] = node_taxid

            # Get current taxon information
            current_name = record.get("ScientificName", "")
            current_rank = record.get("Rank", "no rank")
            
            # Add current taxon to complete lineage
            complete_lineage = lineage + [current_name]            
            
            if current_rank != "no rank":
                rank_info[current_name] = current_rank
                taxid_info[current_name] = taxid
                
            logger.info(f"Successfully retrieved taxonomy information for {taxid}")
            logger.debug(f"Lineage: {complete_lineage}")
            logger.debug(f"Rank info: {rank_info}")
            
            return complete_lineage, rank_info, current_rank, taxid_info
                
        except IncompleteRead as e:
            logger.warning(f"IncompleteRead error for TaxID {taxid}: {e}")
            # Re-raise to allow the retry decorator to handle it
            raise
        except Exception as e:
            if isinstance(e, HTTPError) and e.code == 400:
                logger.error(f"HTTP 400 error for TaxID {taxid}, skipping for now")
            else:
                logger.error(f"Error fetching taxonomy for TaxID {taxid}: {e}")
                logger.error("Full error details:", exc_info=True)
            return [], {}, "", {}

    def try_fetch_at_taxid(self, current_taxid: str, rank_name: str, taxon_name: str,                          
                        sequence_type: str, gene_name: str,
                        protein_records: List[SeqRecord],      
                        nucleotide_records: List[SeqRecord],
                        best_taxonomy: List[str],
                        best_matched_rank: Optional[str],
                        fetch_all: bool = False) -> Tuple[bool, bool, List[str], Optional[str], List[SeqRecord], List[SeqRecord]]:
                protein_found = False
                nucleotide_found = False
                
                # Set minimum protein size for single mode
                min_protein_size = self.config.min_protein_size_single_mode if fetch_all else 0
                
                try:
                    # Handle protein search for 'protein' or 'both' types
                    if sequence_type in ['protein', 'both'] and (not protein_records or fetch_all):
                        # Modify search string based on fetch_all mode
                        if fetch_all:
                            protein_search = f"{self.config.gene_search_term} AND txid{current_taxid}[Organism:exp]"
                        else:
                            protein_search = (f"{self.config.gene_search_term} AND txid{current_taxid}[Organism:exp] "
                                            f"AND {self.config.protein_length_threshold}:10000[SLEN]")
                                                             
                        logger.info(f"Searching protein database at rank {rank_name} ({taxon_name}) with term: {protein_search}")
                                        
                        try:
                            protein_results = self.entrez.search(db="protein", term=protein_search)
                            if protein_results and protein_results.get("IdList"):
                                id_list = protein_results.get("IdList")
                                logger.info(f"Found {len(id_list)} protein IDs")
                                if len(id_list) > 5:  # Only log IDs if there are not too many
                                    logger.debug(f"Protein IDs: {id_list}")
                                
                                # For non-fetch_all mode, apply prefiltering if there are many IDs
                                processed_ids = id_list
                                if not fetch_all and len(id_list) > 10:
                                    logger.info(f"Prefiltering {len(id_list)} proteins based on length information")
                                    
                                    # Get summaries and sort by length
                                    try:
                                        sorted_summaries = []
                                        batch_size = 200
                                        
                                        for i in range(0, len(id_list), batch_size):
                                            batch_ids = id_list[i:i+batch_size]
                                            id_string = ','.join(batch_ids)
                                            
                                            logger.debug(f"Fetching summary for batch of {len(batch_ids)} IDs")
                                            try:
                                                handle = Entrez.esummary(db="protein", id=id_string)
                                                batch_summaries = Entrez.read(handle)
                                                handle.close()
                                                
                                                # Extract sequence lengths from summaries
                                                for summary in batch_summaries:
                                                    seq_id = summary.get('Id', '')
                                                    seq_length = int(summary.get('Length', 0))
                                                    sorted_summaries.append((seq_id, seq_length))
                                                
                                                # Add delay between batches
                                                if i + batch_size < len(id_list):
                                                    sleep(uniform(0.5, 1.0))
                                            except Exception as batch_e:
                                                logger.error(f"Error in batch summary fetch: {batch_e}")
                                                continue
                                        
                                        # Check if we got any summaries
                                        if not sorted_summaries:
                                            logger.error("Failed to fetch any sequence summaries, using all IDs")
                                        else:
                                            # Sort by length (descending)
                                            sorted_summaries.sort(key=lambda x: x[1], reverse=True)
                                            
                                            # Take only top 50 IDs by sequence length
                                            processed_ids = [item[0] for item in sorted_summaries[:10]]
                                            logger.info(f"Successfully filtered to top 10 proteins by length (longest: {sorted_summaries[0][1]} aa)")
                                        
                                    except Exception as e:
                                        logger.error(f"Error in prefiltering: {e}")
                                        logger.error("Full error details:", exc_info=True)
                                        logger.warning("Using all IDs without length filtering")
                                
                                # Log how many IDs we're processing
                                logger.info(f"Processing {len(processed_ids)} protein IDs")
                                
                                # Process the filtered or complete ID list
                                for protein_id in processed_ids:
                                    # Add logging for protein fetch attempt
                                    logger.info(f"Attempting to fetch protein sequence for {gene_name} (accession {protein_id})")
                                    
                                    handle = self.entrez.fetch(db="protein", id=protein_id,
                                                             rettype="gb", retmode="text")
                                    if handle:
                                        temp_record = next(SeqIO.parse(handle, "genbank"))
                                        handle.close()
                                        
                                        # Add logging for successful protein fetch
                                        logger.info(f"Successfully fetched protein sequence of length {len(temp_record.seq)} (accession {temp_record.id})")

                                        # Only skip UniProt records in non-single mode
                                        if not fetch_all:
                                            # Skip problematic UniProt/Swiss-Prot protein accession numbers
                                            if re.match(r'^[A-Z]\d+', temp_record.id) and not re.match(r'^[A-Z]{2,}', temp_record.id):
                                                logger.info(f"Skipping UniProtKB/Swiss-Prot protein record {temp_record.id}")
                                                continue
                                        
                                        # Check minimum protein size in single mode
                                        if fetch_all and len(temp_record.seq) < min_protein_size:
                                            logger.warning(f"Protein sequence too short ({len(temp_record.seq)} aa < {min_protein_size} aa) - skipping (accession {temp_record.id})")
                                            continue
                                        
                                        if fetch_all:
                                            protein_records.append(temp_record)
                                            protein_found = True
                                            if not best_taxonomy:
                                                best_taxonomy = temp_record.annotations.get("taxonomy", [])
                                        else:
                                            # Keep only longest sequence
                                            if not protein_records or len(temp_record.seq) > len(protein_records[0].seq):
                                                protein_records.clear()
                                                protein_records.append(temp_record)
                                                protein_found = True
                                                best_taxonomy = temp_record.annotations.get("taxonomy", [])

                                    # For normal mode (--type both), try to fetch corresponding nucleotide
                                    if protein_found and not fetch_all and sequence_type == 'both':
                                        nucleotide_record = self.fetch_nucleotide_from_protein(protein_records[0], gene_name)
                                        if nucleotide_record:
                                            nucleotide_records.clear()
                                            nucleotide_records.append(nucleotide_record)
                                            nucleotide_found = True
                                            logger.info(f"Successfully fetched corresponding nucleotide sequence")
                                        else:
                                            logger.warning("Failed to fetch corresponding nucleotide sequence")
                                            protein_records.clear()
                                            protein_found = False
                                                
                        except Exception as e:
                            logger.error(f"Error searching protein database: {e}")

                    # Handle nucleotide search
                    if ((sequence_type == 'nucleotide') or 
                        (sequence_type == 'both' and fetch_all) or  # Single taxid mode
                        (sequence_type == 'both' and not nucleotide_found)):  # Fallback for normal mode
                        
                        # Modify search string based on fetch_all mode
                        if fetch_all:
                            nucleotide_search = f"{self.config.gene_search_term} AND txid{current_taxid}[Organism:exp]"
                        else:
                            nucleotide_search = (f"{self.config.gene_search_term} AND txid{current_taxid}[Organism:exp] "
                                                f"AND {self.config.nucleotide_length_threshold}:30000[SLEN]")
                                                                 
                        logger.info(f"Searching nucleotide database at rank {rank_name} ({taxon_name}) with term: {nucleotide_search}")
                        
                        try:
                            nucleotide_results = self.entrez.search(db="nucleotide", term=nucleotide_search)
                            if nucleotide_results and nucleotide_results.get("IdList"):
                                id_list = nucleotide_results.get("IdList")
                                logger.info(f"Found {len(id_list)} nucleotide sequence IDs")
                                if len(id_list) > 5:  # Only log IDs if there are not too many
                                    logger.debug(f"Nucleotide IDs: {id_list}")

                                # Apply the same prefiltering optimization for nucleotide sequences
                                processed_ids = id_list
                                if not fetch_all and len(id_list) > 10:
                                    logger.info(f"Prefiltering {len(id_list)} nucleotide sequences based on length information")
                                    
                                    # Get summaries and sort by length
                                    try:
                                        sorted_summaries = []
                                        batch_size = 200
                                        
                                        for i in range(0, len(id_list), batch_size):
                                            batch_ids = id_list[i:i+batch_size]
                                            id_string = ','.join(batch_ids)
                                            
                                            logger.debug(f"Fetching summary for batch of {len(batch_ids)} IDs")
                                            try:
                                                handle = Entrez.esummary(db="nucleotide", id=id_string)
                                                batch_summaries = Entrez.read(handle)
                                                handle.close()
                                                
                                                # Extract sequence lengths from summaries
                                                for summary in batch_summaries:
                                                    seq_id = summary.get('Id', '')
                                                    seq_length = int(summary.get('Length', 0))
                                                    sorted_summaries.append((seq_id, seq_length))
                                                
                                                # Add delay between batches
                                                if i + batch_size < len(id_list):
                                                    sleep(uniform(0.5, 1.0))
                                            except Exception as batch_e:
                                                logger.error(f"Error in batch summary fetch: {batch_e}")
                                                continue
                                        
                                        # Check if we got any summaries
                                        if not sorted_summaries:
                                            logger.error("Failed to fetch any sequence summaries, using all IDs")
                                        else:
                                            # Sort by length (descending)
                                            sorted_summaries.sort(key=lambda x: x[1], reverse=True)
                                            
                                            # Take only top 50 IDs by sequence length
                                            processed_ids = [item[0] for item in sorted_summaries[:10]]
                                            logger.info(f"Successfully filtered to top 10 nucleotide sequences by length (longest: {sorted_summaries[0][1]} bp)")
                                        
                                    except Exception as e:
                                        logger.error(f"Error in nucleotide prefiltering: {e}")
                                        logger.error("Full error details:", exc_info=True)
                                        logger.warning("Using all IDs without length filtering")
                                
                                # Log how many IDs we're processing
                                logger.info(f"Processing {len(processed_ids)} nucleotide IDs")

                                for seq_id in processed_ids:
                                    try:
                                        logger.info(f"Attempting to fetch nucleotide sequence (accession {seq_id})")
                                        temp_record = self.fetch_nucleotide_record(seq_id)
                                        
                                        if temp_record:
                                            logger.info(f"Successfully fetched nucleotide sequence of length {len(temp_record.seq)} (accession {temp_record.id})")
                                            
                                            if gene_name not in self.config._protein_coding_genes:
                                                # For rRNA genes, use full sequence
                                                if fetch_all:
                                                    nucleotide_records.append(temp_record)
                                                    nucleotide_found = True
                                                    if not best_taxonomy:
                                                        best_taxonomy = temp_record.annotations.get("taxonomy", [])
                                                else:
                                                    # Keep only longest sequence
                                                    if not nucleotide_records or len(temp_record.seq) > len(nucleotide_records[0].seq):
                                                        nucleotide_records.clear()
                                                        nucleotide_records.append(temp_record)
                                                        nucleotide_found = True
                                                        best_taxonomy = temp_record.annotations.get("taxonomy", [])
                                            else:
                                                # For protein-coding genes, extract CDS
                                                logger.info(f"Attempting to extract CDS from nucleotide sequence (accession {temp_record.id})")
                                                cds_record = self.extract_nucleotide(temp_record, gene_name, fetch_all)
                                                if cds_record:
                                                    logger.info(f"Successfully extracted CDS of length {len(cds_record.seq)} (accession {temp_record.id})")
                                                    if fetch_all:
                                                        nucleotide_records.append(cds_record)
                                                        nucleotide_found = True
                                                        if not best_taxonomy:
                                                            best_taxonomy = temp_record.annotations.get("taxonomy", [])
                                                    else:
                                                        # Keep only longest CDS
                                                        if not nucleotide_records or len(cds_record.seq) > len(nucleotide_records[0].seq):
                                                            nucleotide_records.clear()
                                                            nucleotide_records.append(cds_record)
                                                            nucleotide_found = True
                                                            best_taxonomy = temp_record.annotations.get("taxonomy", [])
                                                else:
                                                    logger.warning(f"Failed to extract CDS from nucleotide sequence (accession {temp_record.id})")

                                    except Exception as e:
                                        logger.error(f"Error processing sequence {seq_id}: {e}")
                                        continue
                        except Exception as e:
                            logger.error(f"Error searching nucleotide database: {e}")
                            nucleotide_results = None

                    if protein_found or nucleotide_found:
                        current_match = f"{rank_name}:{taxon_name}" if rank_name else f"exact match:{taxon_name}"
                        if not best_matched_rank or (rank_name and not best_matched_rank.startswith("exact")):
                            best_matched_rank = current_match

                except Exception as e:
                    logger.error(f"Error in try_fetch_at_taxid for taxid {current_taxid}: {e}")
                    logger.error("Full error details:", exc_info=True)

                return protein_found, nucleotide_found, best_taxonomy, best_matched_rank, protein_records, nucleotide_records

    def search_and_fetch_sequences(self, taxid: str, gene_name: str, sequence_type: str, fetch_all: bool = False) -> Tuple[List[SeqRecord], List[SeqRecord], List[str], str]:  
        # Initialise empty lists for records
        protein_records = []
        nucleotide_records = []
        best_taxonomy = []
        best_matched_rank = None

        # Fetch taxonomy first
        taxonomy, taxon_ranks, initial_rank, taxon_ids = self.fetch_taxonomy(taxid)
        if not taxonomy:
            logger.error(f"Could not fetch taxonomy for TaxID {taxid}")
            return [], [], [], "No taxonomy found"

        # Get ordered list of ranks to traverse
        current_taxonomy = taxonomy[:]
        current_taxon = current_taxonomy.pop()  # Start with species
        current_rank = taxon_ranks.get(current_taxon, 'unknown')
        current_taxid = taxid

        # Traverse taxonomy from species up
        while True:
            logger.info(f"Attempting search at {current_rank} level: {current_taxon} (taxid: {current_taxid})")
            
            protein_found, nucleotide_found, best_taxonomy, best_matched_rank, protein_records, nucleotide_records = \
                self.try_fetch_at_taxid(
                    current_taxid, current_rank, current_taxon,
                    sequence_type, gene_name,
                    protein_records, nucleotide_records,
                    best_taxonomy, best_matched_rank,
                    fetch_all
                )

            # For single taxid mode with fetch_all, we only search at the exact taxid level
            if fetch_all:
                break

            # For normal mode, continue searching up taxonomy if needed
            if sequence_type == 'both':
                if protein_records and nucleotide_records:  # Need both in normal mode
                    break
            elif sequence_type == 'protein' and protein_records:
                break
            elif sequence_type == 'nucleotide' and nucleotide_records:
                break

            # Stop if we've gone too high in taxonomy or no more levels
            if (current_rank in ['class', 'subphylum', 'phylum', 'kingdom', 'superkingdom'] or 
                current_taxon == 'cellular organisms' or 
                not current_taxonomy):
                logger.info(f"Reached {current_rank} rank, stopping traversal")
                break

            # Get next level in taxonomy
            current_taxon = current_taxonomy.pop()
            current_rank = taxon_ranks.get(current_taxon, 'unknown')
            current_taxid = taxon_ids.get(current_taxon)
            if not current_taxid:
                continue

            # Add delay between attempts
            sleep(uniform(1, 2))

        # Set final matched rank
        matched_rank = best_matched_rank if best_matched_rank else "No match"
        
        # Different return logic based on mode
        if fetch_all:
            # Single taxid mode: return whatever we found
            if not protein_records and not nucleotide_records:
                logger.warning("No sequences found")
                return [], [], [], "No match"
            logger.info(f"Single taxid mode: Found {len(protein_records)} protein and {len(nucleotide_records)} nucleotide sequences")
            return protein_records, nucleotide_records, best_taxonomy, matched_rank
        else:
            # Normal mode: require both for 'both' type
            if sequence_type == 'both' and (not protein_records or not nucleotide_records):
                logger.warning("Failed to find both protein and corresponding nucleotide sequence")
                return [], [], [], "No match"
            elif sequence_type == 'protein' and not protein_records:
                logger.warning("No protein sequence found")
                return [], [], [], "No match"
            elif sequence_type == 'nucleotide' and not nucleotide_records:
                logger.warning("No nucleotide sequence found")
                return [], [], [], "No match"

        logger.info(f"Search completed. Matched at rank: {matched_rank}")
        return protein_records, nucleotide_records, best_taxonomy, matched_rank

def process_single_taxid(taxid: str, gene_name: str, sequence_type: str,
                        processor: SequenceProcessor, output_dir: Path) -> None:
        """Process a single taxid, fetching all available sequences."""
        try:
            # Fetch all sequences
            protein_records, nucleotide_records, taxonomy, matched_rank = processor.search_and_fetch_sequences(
                taxid, gene_name, sequence_type, fetch_all=True)
            if not protein_records and not nucleotide_records:
                logger.warning(f"No sequences found for taxid {taxid}")
                return
            # Create output manager
            output_manager = OutputManager(output_dir)
            # Save protein sequences
            if sequence_type in ['protein', 'both'] and protein_records:
                for i, record in enumerate(protein_records):
                    filename = f"{record.id}.fasta"
                    output_path = output_dir / filename
                    SeqIO.write(record, output_path, "fasta")
                    logger.info(f"Written protein sequence {i+1}/{len(protein_records)} to '{output_path}'")
                
                # Use the output manager's method to save summary
                output_manager.save_sequence_summary(protein_records, "protein")
                logger.info(f"###############Saved summary of {len(protein_records)} protein sequences###############")
            # Save nucleotide sequences
            if sequence_type in ['nucleotide', 'both'] and nucleotide_records:
                nucleotide_dir = output_dir / 'nucleotide'
                ensure_directory(nucleotide_dir)
                for i, record in enumerate(nucleotide_records):
                    filename = f"{record.id}.fasta"
                    output_path = nucleotide_dir / filename
                    SeqIO.write(record, output_path, "fasta")
                    logger.info(f"Written nucleotide sequence {i+1}/{len(nucleotide_records)} to '{output_path}'")
                
                # Use the output manager's method to save summary
                output_manager.save_sequence_summary(nucleotide_records, "nucleotide")
                logger.info(f"###############Saved summary of {len(nucleotide_records)} nucleotide sequences###############")
        except Exception as e:
            logger.error(f"Error processing taxid {taxid}: {e}")
            logger.error("Full error details:", exc_info=True)

class OutputManager:
    """Manages output files and directories."""
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.nucleotide_dir = output_dir / 'nucleotide'
        self.failed_searches_path = output_dir / "failed_searches.csv"
        self.sequence_refs_path = output_dir / "sequence_references.csv"
        
        self._setup_directories()
        self._setup_files()
        
    def _setup_directories(self):
        """Ensure all required directories exist."""
        ensure_directory(self.output_dir)
        ensure_directory(self.nucleotide_dir)
        
    def _setup_files(self):
        """Initialise required files if they don't exist."""
        if not self.failed_searches_path.exists():
            with open(self.failed_searches_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['process_id', 'taxid', 'error_type', 'timestamp'])
                
        if not self.sequence_refs_path.exists():
            with open(self.sequence_refs_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([
                    'process_id', 'taxid', 'protein_accession', 'protein_length',
                    'nucleotide_accession', 'nucleotide_length', 'matched_rank',
                    'ncbi_taxonomy', 'reference_name', 'protein_reference_path',
                    'nucleotide_reference_path'
                ])

    def log_failure(self, process_id: str, taxid: str, error_type: str):
        """Log failed searches."""
        with open(self.failed_searches_path, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([process_id, taxid, error_type, 
                           time.strftime("%Y-%m-%d %H:%M:%S")])

    def write_sequence_reference(self, data: Dict[str, Any]):
        """Write sequence reference information."""
        with open(self.sequence_refs_path, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                data['process_id'], data['taxid'],
                data.get('protein_accession', ''),
                data.get('protein_length', ''),
                data.get('nucleotide_accession', ''),
                data.get('nucleotide_length', ''),
                data.get('matched_rank', 'unknown'),
                data.get('taxonomy', ''),
                data['process_id'],
                data.get('protein_path', ''),
                data.get('nucleotide_path', '')
            ])

    def save_sequence_summary(self, sequences: List[SeqRecord], file_type: str):
        """
        Save summary of sequences to CSV file.
        
        Args:
            sequences: List of sequence records
            file_type: Either 'protein' or 'nucleotide'
        """
        if not sequences:
            logger.info(f"No {file_type} sequences to summarize")
            return
            
        file_path = self.output_dir / f"fetched_{file_type}_sequences.csv"
        
        with open(file_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Accession', 'Length', 'Description'])
            
            for record in sequences:
                # Get the full description/name
                description = record.description
                accession = record.id
                length = len(record.seq)
                
                writer.writerow([accession, length, description])
                
        logger.info(f"Wrote {file_type} sequence summary to {file_path}")

def get_process_id_column(header):
    """Identify the process ID column from possible variations."""
    # Debug: print what we actually see in the header
    logger.info(f"CSV header detected: {header}")
    
    valid_names = ['ID', 'process_id', 'Process ID', 'process id', 'Process id', 
                  'PROCESS ID', 'sample', 'SAMPLE', 'Sample']
    
    # Debug: print repr of each header item to see invisible characters
    for i, col in enumerate(header):
        logger.info(f"Header column {i}: {repr(col)}")
    
    # Try direct comparison first
    for col in header:
        if col in valid_names:
            logger.info(f"Found matching column: {col}")
            return col
    
    # Try trimming whitespace (in case of spaces)
    for col in header:
        trimmed = col.strip()
        if trimmed in valid_names:
            logger.info(f"Found matching column after trimming: {trimmed}")
            return col
    
    # Last resort: try case-insensitive comparison
    for col in header:
        if col.upper() in [name.upper() for name in valid_names]:
            logger.info(f"Found matching column case-insensitive: {col}")
            return col
    
    logger.error(f"No matching column found in {header}")
    return None

def process_sample(process_id: str, taxid: str, sequence_type: str, 
                  processor: SequenceProcessor, output_manager: OutputManager,
                  gene_name: str) -> None:
    """Process a single sample."""
    try:
        # Define output paths
        protein_path = output_manager.output_dir / f"{process_id}.fasta"
        nucleotide_path = output_manager.nucleotide_dir / f"{process_id}_dna.fasta"
        
        # Check if files already exist
        if ((sequence_type in ['protein', 'both'] and protein_path.exists()) or
            (sequence_type in ['nucleotide', 'both'] and nucleotide_path.exists())):
            logger.info(f"Sequence file(s) already exist for {process_id}. Skipping.")
            return
            
        # Fetch sequences (returns lists)
        protein_records, nucleotide_records, taxonomy, matched_rank = processor.search_and_fetch_sequences(
            taxid, gene_name, sequence_type)
            
        # Extract single records from lists for normal mode
        protein_record = protein_records[0] if protein_records else None
        nucleotide_record = nucleotide_records[0] if nucleotide_records else None
        
        sequences_found = False
        result_data = {
            'process_id': process_id,
            'taxid': taxid,
            'matched_rank': matched_rank,
            'taxonomy': "; ".join(taxonomy) if taxonomy else ""
        }

        # Process protein sequence
        if protein_record and sequence_type in ['protein', 'both']:
            try:
                result_data['protein_accession'] = protein_record.id
                result_data['protein_length'] = len(protein_record.seq)
                result_data['protein_path'] = str(protein_path.absolute())
                
                protein_record.id = process_id
                protein_record.description = ""
                SeqIO.write(protein_record, protein_path, "fasta")
                logger.info(f"Written protein sequence to '{protein_path}'")
                sequences_found = True
            except Exception as e:
                logger.error(f"Error writing protein sequence: {e}")
                if protein_path.exists():
                    protein_path.unlink()

        # Process nucleotide sequence
        if nucleotide_record and sequence_type in ['nucleotide', 'both']:
            try:
                result_data['nucleotide_accession'] = nucleotide_record.id
                result_data['nucleotide_length'] = len(nucleotide_record.seq)
                result_data['nucleotide_path'] = str(nucleotide_path.absolute())
                
                nucleotide_record.id = process_id
                nucleotide_record.description = ""
                SeqIO.write(nucleotide_record, nucleotide_path, "fasta")
                logger.info(f"Written nucleotide sequence to '{nucleotide_path}'")
                sequences_found = True
            except Exception as e:
                logger.error(f"Error writing nucleotide sequence: {e}")
                if nucleotide_path.exists():
                    nucleotide_path.unlink()

        if sequences_found:
            output_manager.write_sequence_reference(result_data)
        else:
            output_manager.log_failure(process_id, taxid, "No sequences found")
            logger.warning(f"No valid sequences found for TaxID {taxid}")

    except Exception as e:
        logger.error(f"Error processing sample {process_id}: {e}")
        output_manager.log_failure(process_id, taxid, f"Processing error: {str(e)}")


def main():
    parser = setup_argument_parser()
    args = parser.parse_args()

    gene_name = args.gene.lower()
    output_dir = Path(args.out)
    sequence_type = args.type.lower()

    # Ensure output directory exists before setting up logging
    ensure_directory(output_dir)
    logger = setup_logging(output_dir) 

    # Initialize components with required email/api_key
    try:
        config = Config(email=args.email, api_key=args.api_key)
        if args.single:
            # Set very low thresholds when in single mode (effectively no threshold)
            config.protein_length_threshold = 0
            config.nucleotide_length_threshold = 0
            logger.info("Single mode activated: sequence length thresholds disabled")
        else:
            config.update_thresholds(args.protein_size, args.nucleotide_size)
            
        search_type = config.set_gene_search_term(gene_name)

        if sequence_type not in config.valid_sequence_types:
            print(f"Invalid sequence type. Choose from: {', '.join(config.valid_sequence_types)}")
            sys.exit(1)
        
        logger.info(f"Using {search_type} search terms for {gene_name}")
        logger.info(f"Output directory: {output_dir}")
        logger.info(f"Sequence type: {sequence_type}")

        # Initialize remaining components
        entrez = EntrezHandler(config)
        processor = SequenceProcessor(config, entrez)

        # Check if we're in single taxid mode
        if args.single:
            logger.info(f"Single taxid mode activated for taxid: {args.single}")
            process_single_taxid(
                taxid=args.single,
                gene_name=gene_name,
                sequence_type=sequence_type,
                processor=processor,
                output_dir=output_dir
            )
            logger.info("Single taxid processing completed")
            sys.exit(0)

        # Regular CSV processing mode
        if not args.input_csv:
            logger.error("Input CSV file is required when not using --single mode")
            sys.exit(1)

        samples_csv = Path(args.input_csv)
        output_manager = OutputManager(output_dir)

        logger.info(f"Starting gene fetch for {gene_name}")
        logger.info(f"Samples file: {samples_csv}")

        try:
            with open(samples_csv, newline='', encoding='utf-8-sig') as f:
                reader = csv.DictReader(f)
                process_id_col = get_process_id_column(reader.fieldnames)
            
                if not process_id_col:
                    logger.error("Could not find process ID column in input CSV.")
                    sys.exit(1)

                # Count total samples
                total_samples = sum(1 for _ in reader)
                f.seek(0)
                next(reader)  # Skip header

                # Initialize progress tracking
                log_progress(0, total_samples)

                # Process each sample
                for i, row in enumerate(reader, 1):
                    try:
                        taxid = row['taxid'].strip()
                        process_id = row[process_id_col].strip()
                        
                        logger.info(f"====== Processing sample {i}/{total_samples}: {process_id} (TaxID: {taxid}) ======")
                        
                        process_sample(
                            process_id=process_id,
                            taxid=taxid,
                            sequence_type=sequence_type,
                            processor=processor,
                            output_manager=output_manager,
                            gene_name=gene_name
                        )
                        
                        # Log progress
                        log_progress(i, total_samples)
                        
                        # Add a small delay between samples
                        sleep(uniform(0.5, 1.0))
                        
                    except Exception as e:
                        logger.error(f"Error processing row {i}: {e}")
                        continue

                # Log final progress
                log_progress(total_samples, total_samples)

        except Exception as e:
            logger.error(f"Fatal error: {e}")
            sys.exit(1)

    except ValueError as e:
        logger.error(str(e))
        sys.exit(1)

    logger.info("Gene fetch completed successfully")

if __name__ == "__main__":
    main()
