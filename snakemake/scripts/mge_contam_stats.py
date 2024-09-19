import os
import csv
import argparse
from Bio import SeqIO
import numpy as np
from concurrent.futures import ProcessPoolExecutor




"""
Usage:
    mge_contam_stats-list_file.py LOG_FILE OUTPUT_FILE

Description:
    Summarise the number of sequences and coverage in each FASTA file specified in the log file.
    The log file should contain paths to FASTA (.fasta or .fas) files.

Arguments:
    LOG_FILE
        The file containing paths to the FASTA files to be processed. Each line in the file should be the path to a FASTA file.

    OUTPUT_FILE
        The name of the output CSV file where the summary statistics will be written.

Example:
    python mge_contam_stats-list_file.py alignment_files.log summary.csv
"""





def process_file(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8', errors='replace') as handle:
            sequences = list(SeqIO.parse(handle, 'fasta'))
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None

    if not sequences:
        print(f"No sequences found in file: {file_path}")
        return None

    # Use a set to store unique sequence headers
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

    return {
        'File Name': os.path.basename(file_path),
        'Number of Sequences': sequence_count,
        'Min Coverage': min_coverage,
        'Max Coverage': max_coverage,
        'Mean Coverage': mean_coverage
    }




def summarise_fasta(log_file, output_file):
    # Read paths to FASTA files from the log file
    try:
        with open(log_file, 'r') as f:
            file_paths = [line.strip() for line in f.readlines() if line.strip().endswith(('.fasta', '.fas'))]
    except Exception as e:
        print(f"Error reading log file {log_file}: {e}")
        return

    if not file_paths:
        print(f"No valid FASTA files found in log file: {log_file}")
        return

    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['File Name', 'Number of Sequences', 'Min Coverage', 'Max Coverage', 'Mean Coverage']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        with ProcessPoolExecutor() as executor:
            results = executor.map(process_file, file_paths)

        for result in results:
            if result:
                writer.writerow(result)

    print(f"CSV summary of the number of sequences and coverage in each FASTA file has been created as '{output_file}'.")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Summarise the number of sequences and coverage in each FASTA file specified in the log file.')
    parser.add_argument('log_file', type=str, help='The log file containing the paths to FASTA files.')
    parser.add_argument('output_file', type=str, help='The output CSV file name.')

    args = parser.parse_args()

    # Basic validation to ensure the log file exists and the output file does not already exist
    if not os.path.isfile(args.log_file):
        parser.error(f"The log file '{args.log_file}' does not exist.")
    if os.path.exists(args.output_file):
        parser.error(f"The file '{args.output_file}' already exists. Choose a different name or remove the existing file.")

    summarise_fasta(args.log_file, args.output_file)
