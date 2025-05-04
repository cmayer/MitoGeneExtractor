#!/usr/bin/env python3
"""
FASTA Sequence Quality Analyser
===============================

This script analyses FASTA files containing DNA sequences, evaluates sequence quality based on 
various metrics, and selects the best sequences according to specific (BOLD BIN) quality criteria.

The script processes both full sequences and barcode regions of different genetic markers (cox1/rbcl/matk),
outputting filtered sequences to separate FASTA files and providing detailed analysis in a CSV file.

Features:
---------
- Multi-file processing: analyse multiple FASTA files in a single run
- Quality assessment: Evaluate sequences based on gaps, ambiguous bases, and continuous stretches
- Sequence selection: Select the 'best' sequences based on configurable ranking criteria
- Target-specific barcode region extraction:
  - cox1/COI: Positions 40-700
  - rbcl: Positions 1-700
  - matk: Positions 1-900
- Detailed reporting: Generate comprehensive CSV report of all sequences and their metrics
- Gap handling: Specialised handling of different types of gaps (-, ~)
- N-trimming: Automatic trimming of leading and trailing N characters

Ranking Criteria:
----------------
Barcode Rank (1-6, lower is better):
1: No ambiguous bases, longest stretch ≥ 650
2: No ambiguous bases, longest stretch ≥ 500
3: No ambiguous bases, 300 ≤ longest stretch ≤ 499
4: No ambiguous bases, 1 ≤ longest stretch ≤ 299
5: Has ambiguous bases
6: Other cases

Full Sequence Rank (1-3, lower is better):
1: No ambiguous bases
2: Has ambiguous bases
3: Other cases

Input Parameters:
----------------
--output-csv/-o: Path where the analysis results CSV file will be saved
--output-fasta/-of: Path where the best full sequences FASTA file will be saved
--output-barcode/-ob: Path where the best barcode sequences FASTA file will be saved
--input/-i: One or more input FASTA files to analyse (space-separated)
--target/-t: Target genetic marker (cox1/COI, rbcl, or matk)

Optional:
--log-file LOG_FILE: Specify a custom path for the log file (default: creates timestamped log in current directory)
--verbose, -v: Enable detailed debug logging

Examples:
--------
Basic usage with cox1 target:
    python fasta_compare.py --output-csv results.csv --output-fasta best.fasta --output-barcode barcode.fasta --input sample1.fasta sample2.fasta --target cox1

Dependencies:
------------
- BioPython: For parsing and manipulating FASTA files
- Standard library: csv, re, os, argparse, logging, etc.

Output Files:
------------
1. CSV report: Contains detailed metrics for all sequences across all input files
2. Best sequences FASTA: Contains the highest quality full sequences (one per process_id)
3. Best barcode FASTA: Contains the highest quality barcode regions (one per process_id)
4. Log file: Records the analysis process, warnings, and errors

Author: 
-------
Created by Ben Price & Dan Parsons @ NHMUK
"""

import sys
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import os
import argparse
from pathlib import Path
import logging
from datetime import datetime


# Define marker-specific barcode regions
BARCODE_REGIONS = {
    'cox1': (39, 700),  # 0-based indexing, so 40-700 becomes 39-700
    'rbcl': (0, 700),   # 1-700 becomes 0-700
    'matk': (0, 900)    # 1-900 becomes 0-900
}

# Define marker aliases
MARKER_ALIASES = {
    'cox1': ['cox1', 'coi', 'co1'],
    'rbcl': ['rbcl'],
    'matk': ['matk']
}


def normalise_marker(marker):
    """
    normalise marker name to handle case insensitivity and aliases.
    
    Parameters:
        marker (str): The marker name provided by the user
        
    Returns:
        str: normalised marker name or None if no match found
    """
    marker = marker.lower()
    
    for normalised, aliases in MARKER_ALIASES.items():
        if marker in aliases:
            return normalised
    
    return None


def setup_logging(log_file=None):
    """
    Set up logging configuration.
    
    Parameters:
        log_file (str): Optional path to log file
    """
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file if log_file else f"fasta_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
        ]
    )


def validate_files(input_files, output_csv, output_fasta):
    """
    Validate input and output file paths.
    
    Parameters:
        input_files (list): List of input FASTA file paths
        output_csv (str): Path to output CSV file
        output_fasta (str): Path to output FASTA file
        
    Returns:
        bool: True if validation passes, False otherwise
    """
    #Check input files exist
    for file in input_files:
        if not os.path.exists(file):
            logging.error(f"Input file does not exist: {file}")
            return False
        
    #Check output directories exist, or create them
    for output_file in [output_csv, output_fasta]:
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
                logging.info(f"Created output directory: {output_dir}")
            except OSError as e:
                logging.error(f"Cannot create output directory {output_dir}: {e}")
                return False
                
    #Check output files don't already exist, if so overwrite
    for output_file in [output_csv, output_fasta]:
        if os.path.exists(output_file):
            logging.warning(f"Output file already exists and will be overwritten: {output_file}")
            
    return True


def trim_n_characters(sequence):
    """
    Trim leading and trailing N characters from a sequence.
    
    Parameters:
        sequence (str): The sequence to trim
        
    Returns:
        str: Trimmed sequence
    """
    # First trim leading N's
    start = 0
    while start < len(sequence) and sequence[start] in ['N', 'n']:
        start += 1
    
    # If sequence is all N's, return empty string
    if start == len(sequence):
        return ""
    
    # Then trim trailing N's
    end = len(sequence) - 1
    while end >= 0 and sequence[end] in ['N', 'n']:
        end -= 1
    
    # Return trimmed sequence
    return sequence[start:end + 1]
    

def analyse_fasta(file_path, target_marker):
    """
    Analyse a FASTA file and extract various metrics from each sequence.

    Parameters:
        file_path (str): The path to the FASTA file to be analysed.
        target_marker (str): The target genetic marker (cox1, rbcl, matk)

    Returns:
        dict: A dictionary containing analysis results for each sequence in the FASTA file.
    """
    try:
        # Get barcode region bounds based on target marker
        barcode_start, barcode_end = BARCODE_REGIONS[target_marker]
        
        # Initialize dictionaries
        results = {}
        sequences = {}
        
        # Count total sequences for progress reporting
        total_sequences = sum(1 for _ in SeqIO.parse(file_path, "fasta"))
        processed = 0

        # Parse the FASTA file
        for record in SeqIO.parse(file_path, "fasta"):
            try:
                processed += 1
                if processed % 100 == 0:  # Log progress every 100 sequences
                    logging.info(f"Processing {file_path}: {processed}/{total_sequences} sequences")

                # Get the sequence ID (header) and the sequence itself
                seq_id = record.id
                seq = str(record.seq)
                
                # Store the complete record
                sequences[seq_id] = record

                # Clean and split seq_id into process_id and parameters
                seq_id_clean = seq_id.replace("Consensus_", "") if seq_id.startswith("Consensus_") else seq_id

                # Initialize process_id and parameters
                process_id = seq_id_clean
                parameters = ""

                # Extract process_id and parameters
                if '_r_' in seq_id_clean:
                    param_start_idx = seq_id_clean.index('_r_')
                    process_id = seq_id_clean[:param_start_idx]
                    parameters = seq_id_clean[param_start_idx + 1:]

                logging.debug(f"Detected process_id: {process_id}, parameters: {parameters}")

                # Calculate sequence metrics
                length = len(seq)
                leading_gaps = len(seq) - len(seq.lstrip('-'))
                trailing_gaps = len(seq) - len(seq.rstrip('-'))
                internal_gaps = seq.count('-') + seq.count('~') - leading_gaps - trailing_gaps
                ambiguous_bases = seq.count('N')

                # Calculate longest stretch without gaps
                longest_stretch = max(len(s) for s in re.split('-|~', seq))

                # Analyse barcode region based on target marker
                subseq = seq[barcode_start:barcode_end]
                barcode_length = len(subseq)
                barcode_ambiguous_bases = subseq.count('N')
                barcode_longest_stretch = max(len(s) for s in re.split('-|~', subseq))

                # Determine ranks
                barcode_rank = calculate_barcode_rank(barcode_ambiguous_bases, barcode_longest_stretch)
                full_rank = calculate_full_rank(ambiguous_bases)

                # Store results
                results[seq_id] = {
                    'file': file_path,
                    'length': length,
                    'leading_gaps': leading_gaps,
                    'trailing_gaps': trailing_gaps,
                    'internal_gaps': internal_gaps,
                    'ambiguous_bases': ambiguous_bases,
                    'longest_stretch': longest_stretch,
                    'barcode_length': barcode_length,
                    'barcode_ambiguous_bases': barcode_ambiguous_bases,
                    'barcode_longest_stretch': barcode_longest_stretch,
                    'barcode_rank': barcode_rank,
                    'full_rank': full_rank,
                    'process_id': process_id,
                    'parameters': parameters,
                    'sequence_record': sequences[seq_id]
                }

            except Exception as e:
                logging.error(f"Error processing sequence {seq_id} in {file_path}: {str(e)}")
                continue

        return results

    except Exception as e:
        logging.error(f"Error analysing file {file_path}: {str(e)}")
        return {}


def calculate_barcode_rank(barcode_ambiguous_bases, barcode_longest_stretch):
    """Calculate barcode rank based on criteria."""
    if barcode_ambiguous_bases == 0:
        if barcode_longest_stretch >= 650:
            return 1
        elif barcode_longest_stretch >= 500:
            return 2
        elif 300 <= barcode_longest_stretch <= 499:
            return 3
        elif 1 <= barcode_longest_stretch <= 299:
            return 4
    return 5 if barcode_ambiguous_bases >= 1 else 6


def calculate_full_rank(ambiguous_bases):
    """Calculate full rank based on criteria."""
    if ambiguous_bases == 0:
        return 1
    elif ambiguous_bases >= 1:
        return 2
    return 3


def format_sequence(seq_record, trim_gaps=True, convert_internal_gaps=True):
    """
    Format a sequence according to specified criteria.
    
    Parameters:
        seq_record (SeqRecord): The sequence record to format
        trim_gaps (bool): Whether to trim leading and trailing gaps
        convert_internal_gaps (bool): Whether to convert internal gaps to N
        
    Returns:
        SeqRecord: Formatted sequence record
    """
    sequence = str(seq_record.seq)
    
    if trim_gaps:
        # Trim leading and trailing gaps
        sequence = sequence.strip('-').strip('~')
    
    if convert_internal_gaps:
        # First remove ~ characters and stitch sequence
        sequence = sequence.replace('~', '')
        # Then replace remaining gaps (-) with N
        sequence = re.sub(r'-', 'N', sequence)
    
    # Trim leading and trailing N's
    sequence = trim_n_characters(sequence)
    
    # Create new sequence record with formatted sequence
    new_record = SeqRecord(
        Seq(sequence),
        id=seq_record.id,
        description=seq_record.description
    )
    
    return new_record


def format_barcode_sequence(seq_record, target_marker):
    """
    Extract and format the barcode region based on the target marker.
    
    Parameters:
        seq_record (SeqRecord): The sequence record to process
        target_marker (str): The target genetic marker (cox1, rbcl, matk)
        
    Returns:
        SeqRecord: Formatted barcode sequence record
    """
    # Get barcode region bounds based on target marker
    barcode_start, barcode_end = BARCODE_REGIONS[target_marker]
    
    # Extract barcode region
    sequence = str(seq_record.seq)[barcode_start:barcode_end]
    
    # Format gaps according to rules
    formatted_sequence = format_barcode_gaps(sequence)
    
    # Trim leading and trailing N's
    formatted_sequence = trim_n_characters(formatted_sequence)
    
    new_record = SeqRecord(
        Seq(formatted_sequence),
        id=seq_record.id,
        description=seq_record.description
    )
    
    return new_record


def format_barcode_gaps(sequence):
    """
    Format gaps in barcode sequence according to specified rules:
    - Remove internal ~ characters and stitch sequence back together
    - For gaps marked with -, if â‰¤ 6 bases fill with N
    - For gaps > 6 bases marked with -, keep the longest fragment
    
    Parameters:
        sequence (str): The sequence to process
        
    Returns:
        str: Formatted sequence containing the longest fragment with ~ removed
    """
    # First handle ~ characters by removing them and stitching sequence
    sequence = sequence.replace('~', '')
    
    # Now handle remaining gaps (-). Split sequence by gaps of length > 6
    fragments = []
    current_fragment = []
    gap_count = 0
    
    for char in sequence:
        if char == '-':
            gap_count += 1
            if gap_count > 6:
                if current_fragment:
                    fragments.append(''.join(current_fragment))
                current_fragment = []
                gap_count = 0
        else:
            if gap_count > 0 and gap_count <= 6:
                # Fill small gaps with N
                current_fragment.extend(['N'] * gap_count)
            gap_count = 0
            current_fragment.append(char)
    
    # Add last fragment if exists
    if current_fragment:
        fragments.append(''.join(current_fragment))
    
    # Find and get longest fragment
    if not fragments:
        return ""
    longest_fragment = max(fragments, key=len)
    return longest_fragment


def format_sequence_id(process_id, parameters):
    """
    Format sequence ID according to specified format.
    
    Parameters:
        process_id (str): Process ID
        parameters (str): Parameters
        
    Returns:
        str: Formatted sequence ID
    """
    return f"{process_id}_{parameters}" if parameters else process_id


def select_sequences_for_process(results):
    """
    Select the best sequences for a process based on ranking criteria.
    First tries to find sequences with full_rank 1, then falls back to full_rank 2 if necessary.
    
    Parameters:
        results (list): List of sequence results for a process
        
    Returns:
        tuple: (best_full_sequence, best_barcode_sequence) or (None, None) if no qualifying sequences
    """
    #First try to find sequences with full_rank 1
    rank1_sequences = [
        result for result in results
        if (result['full_rank'] == 1 and result['barcode_rank'] in [1, 2, 3])
    ]
    
    #If we found rank 1 sequences, use those
    if rank1_sequences:
        #Sort by criteria
        best_sequence = min(rank1_sequences, key=lambda x: (
            x['barcode_rank'],
            -x['barcode_longest_stretch'],
            -x['longest_stretch']
        ))
        return best_sequence, best_sequence
    
    #If no rank 1 sequences, try rank 2 sequences for barcode
    rank2_sequences = [
        result for result in results
        if (result['full_rank'] == 2 and result['barcode_rank'] in [1, 2, 3])
    ]
    
    if rank2_sequences:
        #Sort by criteria
        best_sequence = min(rank2_sequences, key=lambda x: (
            x['barcode_rank'],
            -x['barcode_longest_stretch'],
            -x['longest_stretch']
        ))
        return None, best_sequence
    
    return None, None


def write_best_sequences(best_sequences, output_fasta, output_barcode_fasta, target_marker):
    """
    Write the best sequences that meet specific criteria to new FASTA files.
    Sequences are written without line breaks.
    
    Parameters:
        best_sequences (dict): Dictionary containing the best sequences for each process_id
        output_fasta (str): Path to output FASTA file for full sequences
        output_barcode_fasta (str): Path to output FASTA file for barcode regions
        target_marker (str): The target genetic marker (cox1, rbcl, matk)
    """
    try:
        selected_full_records = []
        selected_barcode_records = []
        selection_records = {}
        
        for process_id, sequence_data in best_sequences.items():
            result = sequence_data[1]  # Get the single best result directly
            
            # Create unique identifier for the sequence
            unique_key = f"{result['file']}_{result['seq_id']}"
            
            #Select best sequences using criteria
            best_full, best_barcode = select_sequences_for_process([result])
            
            if best_full or best_barcode:
                #Process full sequence, if we have one
                if best_full:
                    seq_record = best_full['sequence_record']
                    new_id = format_sequence_id(best_full['process_id'], best_full['parameters'])
                    
                    full_seq_record = format_sequence(seq_record)
                    full_seq_record.id = new_id
                    full_seq_record.description = ""
                    selected_full_records.append(full_seq_record)
                    
                    selection_records[process_id] = {
                        'full_selected_seq': unique_key,
                        'barcode_selected_seq': unique_key if best_barcode else None
                    }
                
                #Process barcode sequence, if we have one
                if best_barcode:
                    seq_record = best_barcode['sequence_record']
                    new_id = format_sequence_id(best_barcode['process_id'], best_barcode['parameters'])
                    
                    barcode_seq_record = format_barcode_sequence(seq_record, target_marker)
                    barcode_seq_record.id = new_id
                    barcode_seq_record.description = ""
                    selected_barcode_records.append(barcode_seq_record)
                    
                    if process_id not in selection_records:
                        selection_records[process_id] = {
                            'full_selected_seq': None,
                            'barcode_selected_seq': unique_key
                        }
                
                logging.info(f"Process {process_id}:")
                if best_full:
                    logging.info(f"  Selected full sequence {best_full['seq_id']} from {best_full['file']} (full_rank: {best_full['full_rank']}, barcode_rank: {best_full['barcode_rank']})")
                if best_barcode and best_barcode != best_full:
                    logging.info(f"  Selected barcode sequence {best_barcode['seq_id']} from {best_barcode['file']} (full_rank: {best_barcode['full_rank']}, barcode_rank: {best_barcode['barcode_rank']})")
        
        #Write selected sequences to output files without line breaks
        def write_fasta_no_wrap(records, filename):
            with open(filename, 'w') as f:
                for record in records:
                    f.write(f">{record.id}\n{str(record.seq)}\n")
                    
        if selected_full_records:
            write_fasta_no_wrap(selected_full_records, output_fasta)
            logging.info(f"Wrote {len(selected_full_records)} sequences to {output_fasta}")
        else:
            logging.warning("No sequences met the criteria for full sequence output")
            
        if selected_barcode_records:
            write_fasta_no_wrap(selected_barcode_records, output_barcode_fasta)
            logging.info(f"Wrote {len(selected_barcode_records)} sequences to {output_barcode_fasta}")
        else:
            logging.warning("No sequences met the criteria for barcode sequence output")
            
        return selection_records
            
    except Exception as e:
        logging.error(f"Error writing sequences: {str(e)}")
        raise


def main():
    # Arg parser
    parser = argparse.ArgumentParser(description='Analyse FASTA files and select best sequences based on quality criteria.')
    
    # Required arguments with flags
    parser.add_argument('--output-csv', '-o', required=True, help='Path to output CSV file')
    parser.add_argument('--output-fasta', '-of', required=True, help='Path to output FASTA file for best sequences')
    parser.add_argument('--output-barcode', '-ob', required=True, help='Path to output FASTA file for barcode regions')
    parser.add_argument('--input', '-i', required=True, nargs='+', help='Input FASTA files to analyse')
    parser.add_argument('--target', '-t', required=True, help='Target genetic marker (cox1/COI/CO1, rbcl/RBCL, or matk/MATK)')
    
    # Optional arguments
    parser.add_argument('--log-file', help='Path to log file (optional)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Setup logging FIRST
    setup_logging(args.log_file)
    
    # Set verbosity level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Normalise target marker
    target_marker = normalise_marker(args.target)
    if not target_marker:
        logging.error(f"Invalid target marker: {args.target}. Must be one of: cox1/COI/CO1, rbcl/RBCL, or matk/MATK")
        sys.exit(1)
    
    # Now log info after logging is set up
    logging.info(f"Using target marker: {target_marker} (barcode region: {BARCODE_REGIONS[target_marker][0]+1}-{BARCODE_REGIONS[target_marker][1]})")
        
    #Setup logging
    setup_logging(args.log_file)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    #Validate input and output files
    if not validate_files(args.input, args.output_csv, args.output_fasta):
        sys.exit(1)
    
    try:
        #Initialise results list
        all_results = []
        
        #Analyse each FASTA file
        for file in args.input:
            logging.info(f"Processing file: {file}")
            results = analyse_fasta(file, target_marker)
            for seq_id, result in results.items():
                result['seq_id'] = seq_id
                all_results.append(result)

        #Select best sequences
        best_sequences = {}
        for result in all_results:
            process_id = result['process_id']
            # Create a unique key combining file path and seq_id
            unique_key = f"{result['file']}_{result['seq_id']}"
            
            criteria = (
                result['barcode_rank'],
                result['full_rank'],
                -result['barcode_longest_stretch'],
                -result['longest_stretch'],
                result['internal_gaps'],
                result['ambiguous_bases'],
                unique_key  # Use unique_key instead of just seq_id
            )
            
            if process_id not in best_sequences or criteria < best_sequences[process_id][0]:
                best_sequences[process_id] = [criteria, result]

        #Annotate best sequences 
        for result in all_results:
            process_id = result['process_id']
            result['best_sequence'] = 'yes' if (
                process_id in best_sequences and 
                result == best_sequences[process_id][1]  # Compare with single result instead of list
            ) else 'no'
        
        #Write best sequences to FASTA files and get selection records
        selection_records = write_best_sequences(best_sequences, args.output_fasta, args.output_barcode, target_marker)

        #Update results with new columns
        for result in all_results:
            process_id = result['process_id']
            unique_key = f"{result['file']}_{result['seq_id']}"
            
            if process_id in selection_records:
                result['selected_full_fasta'] = 'yes' if selection_records[process_id]['full_selected_seq'] == unique_key else 'no'
                result['selected_barcode_fasta'] = 'yes' if selection_records[process_id]['barcode_selected_seq'] == unique_key else 'no'
            else:
                result['selected_full_fasta'] = 'no'
                result['selected_barcode_fasta'] = 'no'

        #Now remove sequence records after FASTA files are written
        for result in all_results:
            if 'sequence_record' in result:
                result.pop('sequence_record', None)

        #Write CSV
        fieldnames = [
            'file', 'seq_id', 'process_id', 'parameters', 'length',
            'leading_gaps', 'trailing_gaps', 'internal_gaps',
            'ambiguous_bases', 'longest_stretch', 'barcode_length',
            'barcode_ambiguous_bases', 'barcode_longest_stretch',
            'barcode_rank', 'full_rank', 'best_sequence',
            'selected_full_fasta', 'selected_barcode_fasta'
        ]

        with open(args.output_csv, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_results)

        logging.info(f"Analysis results saved to {args.output_csv}")
     
    except Exception as e:
        logging.error(f"An error occurred during execution: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
