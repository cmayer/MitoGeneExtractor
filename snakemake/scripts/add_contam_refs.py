#!/usr/bin/env python3
import os
import argparse
from Bio import SeqIO



"""
Usage:
    add_sequences.py --target_dir TARGET_DIR --additions_file ADDITIONS_FILE [--output_dir OUTPUT_DIR]

Description:
    Add sequences from one FASTA file to multiple other FASTA files. The script reads all FASTA files
    in the specified target directory, appends sequences from a specified additions file to each of them,
    and saves the modified FASTA files to the specified output directory.

Arguments:
    --target_dir TARGET_DIR
        Directory containing the target FASTA files to which sequences will be added.

    --additions_file ADDITIONS_FILE
        Path to the FASTA file containing the sequences to be added to each target FASTA file.

    --output_dir OUTPUT_DIR
        Directory where the modified FASTA files will be saved. Defaults to 'output' if not provided.

Example:
    python add_contaminants_fasta.py --target_dir ./data/ --additions_file ./additions.fasta --output_dir ./results/
"""



def add_sequences(target_dir, additions_file, output_dir):
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read the sequences from the additions file
    additions = list(SeqIO.parse(additions_file, "fasta"))

    # Iterate over the files in the target directory
    for filename in os.listdir(target_dir):
        if filename.endswith(".fasta"):
            # Read the sequences from the target file
            target = list(SeqIO.parse(os.path.join(target_dir, filename), "fasta"))

            # Add the sequences from the additions file
            target.extend(additions)

            # Write the modified sequences to a new file in the output directory
            SeqIO.write(target, os.path.join(output_dir, filename), "fasta")


# Define the command-line arguments
parser = argparse.ArgumentParser(description="Add sequences from one fasta file to multiple other fasta files.")
parser.add_argument("--target_dir", required=True, help="Directory containing the target fasta files.")
parser.add_argument("--additions_file", required=True, help="Path to the fasta file containing the sequences to add.")
parser.add_argument("--output_dir", default="output", help="Directory to save the modified fasta files. Defaults to 'output'.")


# Parse the command-line arguments
args = parser.parse_args()


# Call the function with the command-line arguments
add_sequences(args.target_dir, args.additions_file, args.output_dir)
