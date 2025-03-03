import os
import csv
import argparse
import re
import logging
from Bio import SeqIO
import numpy as np
import glob

"""
FASTA Sequence Analysis and Summary Tool
--------------------------------------

This script analyses a log file containing a list of alignment FASTA files, .out files
containing 'raw' summary stats from MGE, and cleaning log files from fasta_cleaner.
It generates comprehensive summary statistics in CSV format. 

Usage:
    python mge_stats.py -a/--alignment_log <log_file> -o/--output <output_csv> -od/--out_file_dir <out_file_dir> -c/--cleaning_logs <cleaning_logs>

Arguments:
    -a, --alignment_log : str
        Path to a text file containing a list of FASTA file paths (one per line)
    -o, --output : str
        Name of the output CSV file
    -od, --out_file_dir : str
        Directory containing .out files with additional sequence statistics
    -c, --cleaning_logs : list[str], optional
        One or more cleaning log files (can use wildcards)

Outputs:
    - <output_file>.csv: Main summary file containing all statistics
    - mge_stats.log: Log file with processing information

The script generates the following metrics for each sample:
    - ID: Identifier extracted from the filename
    - mge_params: Parameters used for MGE (e.g., r_1.3_s_50)
    - n_reads_in: Number of input sequences (from .out file)
    - n_reads_aligned: Number of aligned sequences in the FASTA file
    - n_reads_skipped: Number of sequences that were successfully aligned but not in the FASTA file
    - ref_length: Length of alignment (from .out file)
    - cov_min: Minimum coverage depth across alignment
    - cov_max: Maximum coverage depth across alignment
    - cov_avg: Average coverage depth across alignment
    - cov_med: Median coverage depth across alignment
    - cleaning_input_seqs: Number of input sequences for cleaning
    - cleaning_kept_seqs: Number of sequences kept after cleaning
    - cleaning_removed_human: Number of sequences removed due to human similarity
    - cleaning_removed_at: Number of sequences removed due to AT content
    - cleaning_removed_outlier: Number of sequences removed as statistical outliers
"""

# Set up logging
logger = logging.getLogger('mge_stats')
logger.setLevel(logging.INFO)
# Use mode='w' to overwrite the log file each time
file_handler = logging.FileHandler('mge_stats.log', mode='w')
file_handler.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
console_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.addHandler(console_handler)

def extract_sample_info(filename):
    """
    Extract the process ID and parameters from the filename.
    E.g., from "BGSNL096-23_r_1.3_s_50_align_BGSNL096-23.fas"
    Returns:
        - base_id: "BGSNL096-23"
        - full_id: "BGSNL096-23_r_1.3_s_50"
        - params: "r_1.3_s_50"
    """
    filename_no_ext = os.path.splitext(filename)[0]
    
    # Extract base ID (everything before first underscore)
    base_id_match = re.match(r'^([^_]+)', filename_no_ext)
    base_id = base_id_match.group(1) if base_id_match else None
    
    # Extract parameters (r_X_s_Y pattern)
    params_match = re.search(r'(r_[0-9.]+_s_[0-9.]+)', filename_no_ext)
    params = params_match.group(1) if params_match else ""
    
    # Combine to form the full ID used for matching files
    full_id = f"{base_id}_{params}" if params else base_id
    
    return base_id, full_id, params

def extract_process_id(filename):
    """Extract the process ID from the filename (for backward compatibility)."""
    base_id, _, _ = extract_sample_info(filename)
    return base_id

def parse_out_file(out_file_path):
    """Parse the .out file for additional statistics."""
    try:
        with open(out_file_path, 'r') as file:
            content = file.read()
    except Exception as e:
        logger.error(f"Error reading .out file {out_file_path}: {e}")
        return None

    process_id_data = {
        'n_reads_in': None,
        'ref_length': None,
        'successful_aligned': None  # New field for calculating n_reads_skipped
    }

    # Use regex for more reliable parsing
    reads_in_match = re.search(r'Number of input sequences considered:\s*(\d+)', content)
    if reads_in_match:
        process_id_data['n_reads_in'] = int(reads_in_match.group(1))
    
    # Parse successful aligned sequences for calculating n_reads_skipped
    aligned_match = re.search(r'Number of input sequences successful aligned with exonerate \(all\)\s*\n\s*to the amino acid sequence found in vulgar file:\s*(\d+)', content)
    if aligned_match:
        process_id_data['successful_aligned'] = int(aligned_match.group(1))
    
    ref_length_match = re.search(r'Length of alignment:\s*(\d+)', content)
    if ref_length_match:
        process_id_data['ref_length'] = int(ref_length_match.group(1))
 
    return process_id_data

def parse_cleaning_log(file_path):
    """Parse the fasta_cleaner log file for filtering statistics."""
    try:
        with open(file_path, 'r') as file:
            content = file.read()
            
            # Extract key statistics using regex
            input_match = re.search(r'Input sequences: (\d+)', content)
            kept_match = re.search(r'Kept sequences: (\d+)', content)
            human_match = re.search(r'Removed \(human_similar\): (\d+)', content)
            at_match = re.search(r'Removed \(at_difference\): (\d+)', content)
            # Add parsing for outlier statistics
            outlier_match = re.search(r'Removed \(statistical_outlier\): (\d+)', content)
            
            cleaning_stats = {
                'input_seqs': int(input_match.group(1)) if input_match else 0,
                'kept_seqs': int(kept_match.group(1)) if kept_match else 0,
                'removed_human': int(human_match.group(1)) if human_match else 0,
                'removed_at': int(at_match.group(1)) if at_match else 0,
                'removed_outlier': int(outlier_match.group(1)) if outlier_match else 0
            }
            
            return cleaning_stats
            
    except Exception as e:
        logger.error(f"Error reading cleaning log file {file_path}: {e}")
        return None

def process_fasta_file(file_path):
    """Process a FASTA file and extract sequence statistics."""
    if os.path.getsize(file_path) == 0:
        base_id, _, _ = extract_sample_info(os.path.basename(file_path))
        return {
            'ID': base_id,
            'n_reads_aligned': 0,
            'cov_min': 0,
            'cov_max': 0,
            'cov_avg': 0,
            'cov_med': 0,
        }

    try:
        with open(file_path, 'r', encoding='utf-8', errors='replace') as handle:
            sequences = list(SeqIO.parse(handle, 'fasta'))
    except Exception as e:
        logger.error(f"Error reading file {file_path}: {e}")
        return None

    if not sequences:
        base_id, _, _ = extract_sample_info(os.path.basename(file_path))
        return {
            'ID': base_id,
            'n_reads_aligned': 0,
            'cov_min': 0,
            'cov_max': 0,
            'cov_avg': 0,
            'cov_med': 0,
        }

    unique_sequences = {}
    for seq in sequences:
        if seq.id not in unique_sequences:
            unique_sequences[seq.id] = seq

    sequence_count = len(unique_sequences)
    coverage = np.zeros(len(next(iter(unique_sequences.values())).seq))
    for seq in unique_sequences.values():
        coverage += np.array([1 if base != '-' else 0 for base in seq.seq])

    min_coverage = np.min(coverage)
    max_coverage = np.max(coverage)
    mean_coverage = np.mean(coverage)
    median_coverage = np.median(coverage)

    base_id, _, _ = extract_sample_info(os.path.basename(file_path))

    return {
        'ID': base_id,
        'n_reads_aligned': sequence_count,
        'cov_min': min_coverage,
        'cov_max': max_coverage,
        'cov_avg': mean_coverage,
        'cov_med': median_coverage,
    }

def summarise_fasta(log_file, output_file, out_file_dir, cleaning_log_patterns=None):
    """Summarise FASTA files based on log file with cleaning statistics."""
    try:
        with open(log_file, 'r') as f:
            file_paths = [line.strip() for line in f if line.strip().endswith(('.fasta', '.fas'))]
    except Exception as e:
        logger.error(f"Error reading log file {log_file}: {e}")
        return

    # Log file information
    logger.info(f"The input log_file contains paths to {len(file_paths)} files for processing")
    
    if not file_paths:
        logger.warning(f"No valid FASTA files found in log file: {log_file}")
        return

    # Extract and log sample names
    sample_info = [extract_sample_info(os.path.basename(path)) for path in file_paths]
    base_ids = [info[0] for info in sample_info if info[0]]
    
    logger.info(f"Sample base IDs found: {', '.join(set(base_ids))}")
    logger.info(f"Total number of samples being processed: {len(file_paths)}")

    # Get alignment stats from .out files
    out_files = [f for f in os.listdir(out_file_dir) if f.endswith('.out')]
    logger.info(f"The out_file_dir contains {len(out_files)} .out files for processing")
    
    out_file_data = {}
    for file in out_files:
        _, full_id, _ = extract_sample_info(file)
        if full_id:
            out_data = parse_out_file(os.path.join(out_file_dir, file))
            if out_data:
                out_file_data[full_id] = out_data
    
    # Get cleaning stats from log files
    cleaning_data = {}
    
    # Function to process cleaning log files
    def process_cleaning_logs(log_files):
        logger.info(f"{len(log_files)} cleaning logs being processed")
        for file in log_files:
            _, full_id, _ = extract_sample_info(os.path.basename(file))
            if full_id:
                cleaning_stats = parse_cleaning_log(file)
                if cleaning_stats:
                    cleaning_data[full_id] = cleaning_stats
    
    # Process cleaning logs from patterns if provided
    if cleaning_log_patterns:
        all_cleaning_logs = []
        for pattern in cleaning_log_patterns:
            matching_files = glob.glob(pattern)
            if matching_files:
                all_cleaning_logs.extend(matching_files)
            else:
                logger.warning(f"No files found matching pattern {pattern}")
        
        if all_cleaning_logs:
            process_cleaning_logs(all_cleaning_logs)
    else:
        # Look in default location
        cleaning_logs_dir = os.path.join(os.path.dirname(out_file_dir), "fasta_cleaner/logs")
        if os.path.exists(cleaning_logs_dir):
            log_files = [f for f in os.listdir(cleaning_logs_dir) if f.endswith('_log.txt')]
            process_cleaning_logs([os.path.join(cleaning_logs_dir, f) for f in log_files])
        else:
            logger.warning(f"Cleaning logs directory not found at {cleaning_logs_dir}")

    # Define fieldnames with updated column headings
    fieldnames = [
        'Filename', 'ID', 'mge_params', 'n_reads_in', 'n_reads_aligned', 'n_reads_skipped', 'ref_length', 
        'cov_min', 'cov_max', 'cov_avg', 'cov_med',
        'cleaning_input_seqs', 'cleaning_kept_seqs', 'cleaning_removed_human', 'cleaning_removed_at', 'cleaning_removed_outlier'
    ]

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for file_path in file_paths:
            base_id, full_id, params = extract_sample_info(os.path.basename(file_path))
            if not base_id:
                continue

            result = process_fasta_file(file_path)
            if result:
                result['Filename'] = os.path.basename(file_path).replace('.fasta', '').replace('.fas', '').replace('_align_', '_')
                result['ID'] = base_id
                result['mge_params'] = params

                # Add alignment stats with updated keys - use full_id for matching
                out_data = out_file_data.get(full_id, {})
                result['n_reads_in'] = out_data.get('n_reads_in', '')
                result['ref_length'] = out_data.get('ref_length', '')
                
                # Calculate n_reads_skipped (successful_aligned - n_reads_aligned)
                successful_aligned = out_data.get('successful_aligned', None)
                n_reads_aligned = result.get('n_reads_aligned', 0)
                
                if successful_aligned is not None and n_reads_aligned is not None:
                    result['n_reads_skipped'] = max(0, successful_aligned - n_reads_aligned)
                else:
                    result['n_reads_skipped'] = ''

                # Add cleaning stats - use full_id for matching
                cleaning_stats = cleaning_data.get(full_id, {})
                result['cleaning_input_seqs'] = cleaning_stats.get('input_seqs', '')
                result['cleaning_kept_seqs'] = cleaning_stats.get('kept_seqs', '')
                result['cleaning_removed_human'] = cleaning_stats.get('removed_human', '')
                result['cleaning_removed_at'] = cleaning_stats.get('removed_at', '')
                result['cleaning_removed_outlier'] = cleaning_stats.get('removed_outlier', '')

                # Write all results to the CSV file
                writer.writerow(result)

    output_file_abs = os.path.abspath(output_file)
    logger.info(f"CSV summary created: {output_file_abs}")
    print(f"CSV summary created: '{output_file}'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Summarise FASTA file statistics including sequence counts, coverage, and cleaning results.')
    parser.add_argument('-a', '--alignment_log', required=True, type=str, help='The log file containing the paths to FASTA files.')
    parser.add_argument('-o', '--output', required=True, type=str, help='The output CSV file name.')
    parser.add_argument('-od', '--out_file_dir', required=True, type=str, help='The directory containing .out files with additional statistics.')
    parser.add_argument('-c', '--cleaning_logs', type=str, nargs='+', help='Optional: path patterns to cleaning log files')

    args = parser.parse_args()

    if not os.path.isfile(args.alignment_log):
        parser.error(f"The log file '{args.alignment_log}' does not exist.")
    if not os.path.isdir(args.out_file_dir):
        parser.error(f"The directory '{args.out_file_dir}' does not exist or is not a directory.")

    summarise_fasta(args.alignment_log, args.output, args.out_file_dir, args.cleaning_logs)
