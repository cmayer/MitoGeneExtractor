import os
import csv
import argparse
import re
from Bio import SeqIO
import numpy as np



"""
FASTA Sequence Analysis and Summary Tool
--------------------------------------

This script analyses a log file containing a list of alignment FASTA files and .out files
containing 'raw' summary stats from MGE, and generates summary statistics in CSV format. 
It processes both regular alignments and potential contaminant sequences,
outputting results to separate CSV files.

Usage:
    python mge_stats.py <log_file> <output_file> <out_file_dir>

Arguments:
    log_file : str
        Path to a text file containing a list of FASTA file paths (one per line)
    output_file : str
        Name of the output CSV file (will also generate a -contaminants.csv variant)
    out_file_dir : str
        Directory containing .out files with additional sequence statistics

Outputs:
    - <output_file>.csv: Main summary file containing all statistics
    - <output_file>-contaminants.csv: Summary of potential contaminant sequences (if used)

The script generates the following metrics for each sample (FASTA alignment file):
    - Process ID: Identifier extracted from the filename
    - n_reads: Number of input sequences (from .out file)
    - n_aligned: Number of aligned sequences in the FASTA file
    - skipped_reads_low_rel: Number of skipped reads (from .out file)
    - length: Length of alignment (from .out file)
    - cov_min: Minimum coverage depth across alignment
    - cov_max: Maximum coverage depth across alignment
    - cov_avg: Average coverage depth across alignment
    - cov_med: Median coverage depth across alignment

Requirements:
    - Python 3.6+
    - BioPython
    - NumPy

Example:
    python mge_stats.py path/to/alignment_log.log path/to/summary_stats.csv path/to/out_files/

Notes:
    - Empty FASTA files or files with no sequences will be processed with zero values
    - The script automatically handles both .fasta and .fas file extensions
    - Existing output files will be overwritten with a warning message
"""







def extract_process_id(filename):
    """Extract the process ID from the filename."""
    filename_no_ext = os.path.splitext(filename)[0]
    match = re.match(r'^([^_]+)', filename_no_ext)
    return match.group(1) if match else None
def parse_out_file(out_file_path):
    """Parse the .out file for additional statistics."""
    try:
        with open(out_file_path, 'r') as file:
            lines = file.readlines()
    except Exception as e:
        print(f"Error reading .out file {out_file_path}: {e}")
        return None

    process_id_data = {
        'n_reads': None,
        'skipped_reads_low_rel': None,
        'length': None
    }

    for line in lines:
        if process_id_data['n_reads'] is None and 'Number of input sequences considered:' in line:
            process_id_data['n_reads'] = int(line.split(':')[-1].strip())
        elif process_id_data['skipped_reads_low_rel'] is None and '# skipped reads due to low rel. score:' in line:
            process_id_data['skipped_reads_low_rel'] = int(line.split(':')[-1].strip())
        elif process_id_data['length'] is None and 'Length of alignment:' in line:
            process_id_data['length'] = int(line.split(':')[-1].strip())
        
        #Break out of the loop if all required values are found
        if all(val is not None for val in process_id_data.values()):
            break
 
    return process_id_data


def process_fasta_file(file_path):
    """Process a FASTA file and extract sequence statistics."""
    if os.path.getsize(file_path) == 0:
        # Return placeholder for empty files
        process_id = extract_process_id(os.path.basename(file_path))
        return {
            'Process ID': process_id,
            'n_aligned': 0,
            'cov_min': 0,
            'cov_max': 0,
            'cov_avg': 0,
            'cov_med': 0,
        }

    try:
        with open(file_path, 'r', encoding='utf-8', errors='replace') as handle:
            sequences = list(SeqIO.parse(handle, 'fasta'))
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None

    if not sequences:
        # No sequences found; treat similarly to empty files
        process_id = extract_process_id(os.path.basename(file_path))
        return {
            'Process ID': process_id,
            'n_aligned': 0,
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

    process_id = extract_process_id(os.path.basename(file_path))

    return {
        'Process ID': process_id,
        'n_aligned': sequence_count,
        'cov_min': min_coverage,
        'cov_max': max_coverage,
        'cov_avg': mean_coverage,
        'cov_med': median_coverage,
    }


def summarise_fasta(log_file, output_file, out_file_dir):
    """Summarize FASTA files based on log file."""
    try:
        with open(log_file, 'r') as f:
            file_paths = [line.strip() for line in f if line.strip().endswith(('.fasta', '.fas'))]
    except Exception as e:
        print(f"Error reading log file {log_file}: {e}")
        return

    if not file_paths:
        print(f"No valid FASTA files found in log file: {log_file}")
        return

    out_file_data = {}
    for file in os.listdir(out_file_dir):
        if file.endswith('.out'):
            process_id = extract_process_id(file)
            if process_id:
                out_data = parse_out_file(os.path.join(out_file_dir, file))
                if out_data:
                    out_file_data[process_id] = out_data

    # Automatically overwrite the file if it exists
    if os.path.exists(output_file):
        print(f"Overwriting the existing file: '{output_file}'.")
        os.remove(output_file)

    # Define two sets of fieldnames: one for the main output and one for contaminants.csv
    main_fieldnames = [
        'Filename', 'Process ID', 'n_reads', 'n_aligned', 'skipped_reads_low_rel', 'length', 
        'cov_min', 'cov_max', 'cov_avg', 'cov_med'
    ]
    
    contaminants_fieldnames = [
        'Filename', 'Process ID', 'n_reads', 'n_aligned', 'cov_min', 'cov_max', 'cov_avg', 'cov_med'
    ]

    # Open both the main output file and contaminants file
    with open(output_file, 'w', newline='') as csvfile, \
         open(f"{os.path.splitext(output_file)[0]}-contaminants.csv", 'w', newline='') as contfile:

        # Write headers for the main file (with all fields)
        writer = csv.DictWriter(csvfile, fieldnames=main_fieldnames)
        writer.writeheader()

        # Write headers for contaminants file (without 'skipped_reads_low_rel' and 'length')
        cont_writer = csv.DictWriter(contfile, fieldnames=contaminants_fieldnames)
        cont_writer.writeheader()

        for file_path in file_paths:
            process_id = extract_process_id(os.path.basename(file_path))
            if not process_id:
                continue

            result = process_fasta_file(file_path)
            if result:
                result['Filename'] = os.path.basename(file_path).replace('.fasta', '').replace('.fas', '').replace('_align_', '_')

                # Fetch the first occurrence of data from the out file data
                out_data = out_file_data.get(process_id, {})
                result['n_reads'] = out_data.get('n_reads', '')
                result['skipped_reads_low_rel'] = out_data.get('skipped_reads_low_rel', '')
                result['length'] = out_data.get('length', '')

                current_count = file_path.count(process_id)
                output_type = "contaminants.csv" if current_count == 1 else "output.csv"

                # Write to contaminants.csv without 'skipped_reads_low_rel' and 'length'
                if current_count == 1:
                    result_for_contaminants = {k: result[k] for k in contaminants_fieldnames}
                    cont_writer.writerow(result_for_contaminants)  # Write only relevant fields
                else:
                    writer.writerow(result)  # Write all fields to the main output file

    print(f"CSV summaries created: '{output_file}' and '{os.path.splitext(output_file)[0]}-contaminants.csv'.")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Summarise the number of sequences and coverage in each FASTA file specified in the log file and extract additional statistics from .out files.')
    parser.add_argument('log_file', type=str, help='The log file containing the paths to FASTA files.')
    parser.add_argument('output_file', type=str, help='The output CSV file name.')
    parser.add_argument('out_file_dir', type=str, help='The directory containing .out files with additional statistics.')

    args = parser.parse_args()

    if not os.path.isfile(args.log_file):
        parser.error(f"The log file '{args.log_file}' does not exist.")
    if os.path.exists(args.output_file):
        parser.error(f"The file '{args.output_file}' already exists. Choose a different name or remove the existing file.")
    if not os.path.isdir(args.out_file_dir):
        parser.error(f"The directory '{args.out_file_dir}' does not exist or is not a directory.")

    summarise_fasta(args.log_file, args.output_file, args.out_file_dir)
