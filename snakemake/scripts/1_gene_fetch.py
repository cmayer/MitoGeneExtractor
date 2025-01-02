#!/usr/bin/env python3
"""
Gene Fetch - NCBI Sequence Retrieval Tool

This script fetches gene sequences from NCBI databases based on taxonomy IDs (taxids).
It can retrieve both protein and nucleotide sequences, with support for various genes
including protein-coding genes (e.g., cox1, cox2, cytb) and rRNA genes (e.g., 16S, 18S).

Key Features:
- Taxonomic traversal: If sequences aren't found at the input taxonomic level (e.g. species), searches up higher ranks
- Robust error handling with automatic retries for NCBI API calls
- Length filtering for both protein and nucleotide sequences
- Support for both protein-coding and rRNA genes
- Automatic CDS extraction for protein-coding genes (if '--type both')
- Detailed logging of all operations with progress tracking
- Rate limiting to comply with NCBI API guidelines
- Efficient processing of repeated taxonomy queries
- Progress tracking with percentage completion

Input:
- CSV file containing taxonomy IDs (must have 'taxid' and 'ID' column)
- Gene name (e.g., 'cox1', '16s')
- Output directory path
- Sequence type ('protein', 'nucleotide', or 'both')
- Optional: Minimum sequence length thresholds

Output:
- FASTA files containing retrieved sequences
- Log file detailing all operations with progress updates
- CSV files tracking successful and failed retrievals
- Reference file mapping sequences to taxonomy

Dependencies:
- Biopython>=1.80: For sequence handling and NCBI interactions
- ratelimit>=2.2.1: For API rate limiting
- Standard library: csv, sys, os, time, logging, functools.lru_cache

Usage:
    python gene_fetch.py <gene_name> <output_directory> <samples.csv> --type <sequence_type>
                        [--protein_size <min_size>] [--nucleotide_size <min_size>]

Example:
    python gene_fetch.py cox1 ./output_dir samples.csv --type both --protein_size 500 --nucleotide_size 1500

Notes:
- For protein-coding genes, 'both' type first fetches protein then corresponding nucleotide
- Progress updates logged every 100 samples by default
- Failed retrievals and errors logged to separate CSV file

Author: D. Parsons
Version: 1.0.3
License: MIT
Last updated: 2024-12-29

To do:
- LRU caching for taxonomy lookups to reduce API calls and improve performance
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

def log_progress(current: int, total: int, interval: int = 100) -> None:
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
    
    parser.add_argument('gene_name',
                      help='Name of gene to search for in NCBI RefSeq database (e.g., cox1)')
    
    parser.add_argument('output_directory',
                      help='Path to directory to save output files')
    
    parser.add_argument('samples_csv',
                      help='Path to input CSV file containing TaxIDs')
    
    parser.add_argument('--type', required=True, choices=['protein', 'nucleotide', 'both'],
                      help='Specify sequence type to fetch')
    
    parser.add_argument('--protein_size', type=int, default=500,
                      help='Minimum protein sequence length (default: 500)')
    
    parser.add_argument('--nucleotide_size', type=int, default=1500,
                      help='Minimum nucleotide sequence length (default: 1500)')
    
    return parser



@dataclass
class Config:
    email: str = "b.price@nhm.ac.uk"
    api_key: str = "82df5a6f5cf735302d3cf1fcf48b206cfe09"
    max_calls_per_second: int = 10
    protein_length_threshold: int = 500
    nucleotide_length_threshold: int = 1500
    valid_sequence_types: set = frozenset({'protein', 'nucleotide', 'both'})
    gene_search_term: str = ""
    
    def __post_init__(self):
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
            # Add more genes as needed
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
                '"ribulose bisphosphate carboxylase large chain"[Protein Name]',
                '"RuBisCO large subunit"[Protein Name]',
                '"rbcL gene"[Gene]',
                '"RBCL gene"[Gene]'
            ],
            'matk': [
                'matK[Gene]',
                'MATK[Gene]',
                '"maturase K"[Protein Name]',
                '"maturase K"[Gene]',
                '"matK gene"[Gene]',
                '"MATK gene"[Gene]'
            ]
            # Add more genes as needed
        }

    def update_thresholds(self, protein_size: int, nucleotide_size: int):
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
        """Execute an Entrez esearch query with enhanced error handling."""
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

        return _do_search(**kwargs)


class SequenceProcessor:
    """Handles sequence processing and validation."""
    def __init__(self, config: Config, entrez: EntrezHandler):
        self.config = config
        self.entrez = entrez
        
    def extract_cds(self, record: SeqRecord, gene_name: str) -> Optional[SeqRecord]:
        """
        Extract CDS region from a sequence record.
        
        Args:
            record: The sequence record to process
            gene_name: Name of the gene to find
            
        Returns:
            Optional[SeqRecord]: Extracted CDS record if found and valid, None otherwise
        """
        logger.info(f"Attempting to extract CDS for gene {gene_name}")
        gene_variations = set()
        if gene_name.lower() in self.config._protein_coding_genes:
            # Strip the [Gene] and [Protein Name] tags for matching
            variations = self.config._protein_coding_genes[gene_name.lower()]
            gene_variations = {v.split('[')[0].strip('"') for v in variations}
        else:
            logger.warning(f"No defined variations for gene {gene_name}")
            gene_variations = {gene_name.lower()}  # Use just the gene name if no variations defined
        
        for feature in record.features:
            if feature.type != "CDS":
                continue
                
            qualifiers = []
            for field in ['gene', 'product', 'note']:
                qualifiers.extend(feature.qualifiers.get(field, []))
                
            logger.debug(f"Found qualifiers: {qualifiers}")
                
            if any(var in qualifier.lower() for qualifier in qualifiers 
                  for var in gene_variations):
                try:
                    cds_record = record[:]
                    cds_record.seq = feature.extract(record.seq)
                    
                    if len(cds_record.seq) >= 100:  # Sanity check for minimum length
                        logger.info(f"Successfully extracted CDS of length {len(cds_record.seq)}")
                        return cds_record
                    else:
                        logger.warning(f"Extracted CDS too short ({len(cds_record.seq)} bp)")
                except Exception as e:
                    logger.error(f"CDS extraction error: {e}")
                    logger.error("Full error details:", exc_info=True)
                    
        logger.debug(f"No valid CDS found for gene {gene_name}")
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
                        logger.error("No valid segments found in coded_by qualifier")
                        return None
                    
                    # Fetch full sequence for first accession
                    first_accession = segments[0][0]
                    logger.info(f"Fetching nucleotide sequence for accession: {first_accession}")
                    
                    handle = self.entrez.fetch(db="nucleotide", id=first_accession, rettype="gb", retmode="text")
                    if not handle:
                        logger.error(f"Failed to fetch nucleotide sequence for accession: {first_accession}")
                        return None
                    
                    nucleotide_record = next(SeqIO.parse(handle, "genbank"))
                    handle.close()
                    
                    # Extract and join all segments
                    complete_sequence = ""
                    for accession, coordinates in segments:
                        if accession != first_accession:
                            handle = self.entrez.fetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
                            if not handle:
                                logger.error(f"Failed to fetch additional sequence: {accession}")
                                continue
                            nucleotide_record = next(SeqIO.parse(handle, "genbank"))
                            handle.close()
                        
                        if coordinates:
                            start, end = coordinates
                            logger.debug(f"Sequence length before extraction: {len(nucleotide_record.seq)}")
                            logger.debug(f"Using array indices [{start-1}:{end}] for slicing")
                            # Fix: Ensure we use the full end coordinate 
                            segment_seq = str(nucleotide_record.seq[start-1:end])
                            logger.debug(f"Extracted sequence length: {len(segment_seq)}")
                            logger.debug(f"Expected length from coordinates: {end - (start-1)}")
                            if len(segment_seq) == 0:
                                logger.error(f"Zero-length sequence extracted using coordinates {start}..{end}")
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
                    logger.info(f"Successfully extracted nucleotide sequence of length {len(complete_sequence)} (Accession: {first_accession})")
                    
                    return new_record
                    
            logger.warning("No valid nucleotide reference found in protein record")
            return None
            
        except Exception as e:
            logger.error(f"Error in fetch_nucleotide_from_protein: {e}")
            logger.error("Full error details:", exc_info=True)
            return None

    @lru_cache(maxsize=1000) 
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

            # Initialize rank information dictionaries
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
                
        except Exception as e:
            if isinstance(e, HTTPError) and e.code == 400:
                logger.error(f"HTTP 400 error for TaxID {taxid}, skipping for now")
            else:
                logger.error(f"Error fetching taxonomy for TaxID {taxid}: {e}")
                logger.error("Full error details:", exc_info=True)
            return [], {}, "", {}


    def try_fetch_at_taxid(self, current_taxid: str, rank_name: str, taxon_name: str,
                          sequence_type: str, gene_name: str,
                          protein_record: Optional[SeqRecord],
                          nucleotide_record: Optional[SeqRecord],
                          best_taxonomy: List[str],
                          best_matched_rank: Optional[str]) -> Tuple[bool, bool, List[str], Optional[str], Optional[SeqRecord], Optional[SeqRecord]]:
        """
        Helper function to attempt sequence fetch at a specific taxonomic level.
        Now fetches the longest available sequence at each rank.
        """
        protein_found = False
        nucleotide_found = False

        try:
            # Try protein first if doing 'protein' or 'both'
            if sequence_type in ['protein', 'both'] and not protein_record:
                protein_search = (f"{self.config.gene_search_term} AND txid{current_taxid}[Organism:exp] "
                               f"AND {self.config.protein_length_threshold}:10000[SLEN]")
                logger.info(f"Searching protein database at rank {rank_name or 'species'} ({taxon_name}) with term: {protein_search}")
                
                protein_results = self.entrez.search(db="protein", term=protein_search)
                
                if protein_results and protein_results.get("IdList"):
                    id_list = protein_results.get("IdList")
                    logger.info(f"Found {len(id_list)} protein IDs: {id_list}")
                    
                    # Find longest protein sequence
                    longest_protein = None
                    max_length = 0
                    
                    for protein_id in id_list:
                        handle = self.entrez.fetch(db="protein", id=protein_id,
                                                 rettype="gb", retmode="text")
                        if handle:
                            temp_record = next(SeqIO.parse(handle, "genbank"))
                            handle.close()

                            # Skip problematic UniProt/Swiss-Prot protein accession numbers (unless they have coded_by)
                            if re.match(r'^[A-Z]\d+', temp_record.id) and not re.match(r'^[A-Z]{2,}', temp_record.id):
                                logger.info(f"Skipping UniProtKB/Swiss-Prot protein record {temp_record.id}")
                                continue
                            
                            seq_length = len(temp_record.seq)
                            if seq_length > max_length:
                                max_length = seq_length
                                longest_protein = temp_record
                                logger.info(f"Fetched longest protein sequence: {seq_length} aa (Accession: {temp_record.id})")
                    
                    if longest_protein:
                        protein_record = longest_protein
                        best_taxonomy = protein_record.annotations.get("taxonomy", [])
                        protein_found = True

            # For nucleotide, prioritise getting it from protein if we're doing 'both'
            if sequence_type in ['nucleotide', 'both'] and not nucleotide_record:
                if sequence_type == 'both':
                    if protein_record:
                        nucleotide_record = self.fetch_nucleotide_from_protein(protein_record, gene_name)
                        if nucleotide_record:
                            nucleotide_found = True
                    else:
                        # Skip nucleotide search until we try finding protein at all ranks
                        return protein_found, nucleotide_found, best_taxonomy, best_matched_rank, protein_record, nucleotide_record
                else:
                    # Direct nucleotide search only if specifically requested
                    nucleotide_search = (f"{self.config.gene_search_term} AND txid{current_taxid}[Organism:exp] "
                                       f"AND {self.config.nucleotide_length_threshold}:30000[SLEN]")
                    logger.info(f"Searching nucleotide database at rank {rank_name or 'species'} ({taxon_name}) with term: {nucleotide_search}")
                    
                    nucleotide_results = self.entrez.search(db="nucleotide", term=nucleotide_search)
                    
                    if nucleotide_results:
                        id_list = nucleotide_results.get("IdList", [])
                        if id_list:
                            logger.info(f"Found {len(id_list)} nucleotide sequence IDs: {id_list}")
                        else:
                            logger.info("Found 0 nucleotide sequence IDs")   
                      
                        if id_list:
                            longest_sequence = None
                            max_length = 0
                            best_temp_taxonomy = []

                            # Examine all sequences to find the longest
                            for seq_id in id_list:
                                try:
                                    handle = self.entrez.fetch(db="nucleotide", id=seq_id,
                                                             rettype="gb", retmode="text")
                                    if handle:
                                        temp_record = next(SeqIO.parse(handle, "genbank"))
                                        handle.close()
                                        
                                        # For rRNA genes, use full sequence
                                        if gene_name not in self.config._protein_coding_genes:
                                            seq_length = len(temp_record.seq)
                                            if seq_length > max_length:
                                                max_length = seq_length
                                                longest_sequence = temp_record
                                                best_temp_taxonomy = temp_record.annotations.get("taxonomy", [])
                                                logger.info(f"Found longer nucleotide sequence: {seq_length} bp (Accession: {temp_record.id})")
                                        else:
                                            # For protein-coding genes, check CDS length
                                            cds_record = self.extract_cds(temp_record, gene_name)
                                            if cds_record:
                                                seq_length = len(cds_record.seq)
                                                if seq_length > max_length:
                                                    max_length = seq_length
                                                    longest_sequence = cds_record
                                                    best_temp_taxonomy = temp_record.annotations.get("taxonomy", [])
                                                    logger.info(f"Found longer CDS sequence: {seq_length} bp (Accession: {temp_record.id})")

                                except Exception as e:
                                    logger.error(f"Error processing sequence {seq_id}: {e}")
                                    continue

                            # Use the longest sequence found
                            if longest_sequence:
                                nucleotide_record = longest_sequence
                                if not best_taxonomy:
                                    best_taxonomy = best_temp_taxonomy
                                nucleotide_found = True
                                logger.info(f"Selected longest nucleotide sequence of length {len(nucleotide_record.seq)} (Accession: {nucleotide_record.id})")

            if protein_found or nucleotide_found:
                current_match = f"{rank_name}:{taxon_name}" if rank_name else f"exact match:{taxon_name}"
                if not best_matched_rank or (rank_name and not best_matched_rank.startswith("exact")):
                    best_matched_rank = current_match

        except Exception as e:
            logger.error(f"Error in try_fetch_at_taxid for taxid {current_taxid}: {e}")
            logger.error("Full error details:", exc_info=True)

        return protein_found, nucleotide_found, best_taxonomy, best_matched_rank, protein_record, nucleotide_record

    def search_and_fetch_sequences(self, taxid: str, gene_name: str, sequence_type: str) -> Tuple[Optional[SeqRecord], Optional[SeqRecord], List[str], str]:
        """
        Search and fetch sequences for a given taxid, traversing up taxonomy if needed.
        Returns (protein_record, nucleotide_record, taxonomy, matched_rank)
        """      
        protein_record = None
        nucleotide_record = None
        best_taxonomy = []
        best_matched_rank = None

        # Fetch taxonomy first
        taxonomy, taxon_ranks, initial_rank, taxon_ids = self.fetch_taxonomy(taxid)
        if not taxonomy:
            logger.error(f"Could not fetch taxonomy for TaxID {taxid}")
            return None, None, [], "No taxonomy found"

        logger.info(f"Full taxonomy for TaxID {taxid}: {', '.join(taxonomy)}")
        logger.debug(f"Rank information: {taxon_ranks}")
        logger.debug(f"TaxID mapping: {taxon_ids}")

        # Get ordered list of ranks to traverse, INCLUDING the species level
        current_taxonomy = taxonomy[:]
        current_taxon = current_taxonomy.pop()  # Start with species
        current_rank = taxon_ranks.get(current_taxon, 'unknown')
        current_taxid = taxid  # Use original taxid for species level

        # Traverse taxonomy from species up
        while True:  # Changed from while current_taxonomy to allow final level check
            logger.info(f"Attempting search at {current_rank} level: {current_taxon} (taxid: {current_taxid})")
            
            # Try fetch at current taxonomic level first
            protein_found, nucleotide_found, best_taxonomy, best_matched_rank, protein_record, nucleotide_record = \
                self.try_fetch_at_taxid(
                    current_taxid, current_rank, current_taxon,
                    sequence_type, gene_name,
                    protein_record, nucleotide_record,
                    best_taxonomy, best_matched_rank
                )

            # Check if we have everything we need
            have_protein = protein_record is not None or sequence_type not in ['protein', 'both']
            have_nucleotide = nucleotide_record is not None or sequence_type not in ['nucleotide', 'both']
            
            # Stop if we have everything we need
            if have_protein and have_nucleotide:
                break

            # Stop if we've gone too high in taxonomy or no more levels to check
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
                logger.debug(f"No taxid found for {current_taxon}, continuing to next rank")
                continue

            # Add delay between attempts
            logger.debug("Adding delay between taxonomy level attempts")
            sleep(uniform(1, 2))

        # Set final matched rank
        matched_rank = best_matched_rank if best_matched_rank else "No match"
        
        if not protein_record and not nucleotide_record:
            logger.warning("No sequences found after traversing taxonomy")
            return None, None, [], "No match"

        logger.info(f"Search completed. Matched at rank: {matched_rank}")
        return protein_record, nucleotide_record, best_taxonomy, matched_rank

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

def get_process_id_column(header):
    """Identify the process ID column from possible variations."""
    valid_names = ['ID', 'process_id', 'Process ID', 'process id', 'Process id', 
                  'PROCESS ID', 'sample', 'SAMPLE', 'Sample']
    for name in valid_names:
        if name in header:
            return name
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

        # Fetch sequences
        protein_record, nucleotide_record, taxonomy, matched_rank = processor.search_and_fetch_sequences(
            taxid, gene_name, sequence_type)

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
    """Main function to coordinate the script execution."""
    parser = setup_argument_parser()
    args = parser.parse_args()

    gene_name = args.gene_name.lower()
    output_dir = Path(args.output_directory)
    samples_csv = Path(args.samples_csv)
    sequence_type = args.type.lower()

    # Initialize components
    config = Config()
    
    # Update thresholds if provided
    config.update_thresholds(args.protein_size, args.nucleotide_size)
    config.set_gene_search_term(gene_name) 

    if sequence_type not in config.valid_sequence_types:
        print(f"Invalid sequence type. Choose from: {', '.join(config.valid_sequence_types)}")
        sys.exit(1)

    # Ensure output directory exists before setting up logging
    ensure_directory(output_dir)
    
    logger = setup_logging(output_dir)
    
    search_type = config.set_gene_search_term(gene_name)
    logger.info(f"Using {search_type} search terms for {gene_name}")
    
    # Initialize remaining components
    entrez = EntrezHandler(config)
    processor = SequenceProcessor(config, entrez)
    output_manager = OutputManager(output_dir)

    logger.info(f"Starting gene fetch for {gene_name}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Samples file: {samples_csv}")
    logger.info(f"Sequence type: {sequence_type}")

    try:
        with open(samples_csv, newline='') as f:
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

    logger.info("Gene fetch completed successfully")

if __name__ == "__main__":
    main()
