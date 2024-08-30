import os
import sys
import pandas as pd
import numpy as np



def extract_data_from_out_files(folder_path, output_csv):
    try:
        data_list = []

        for filename in os.listdir(folder_path):
            if filename.endswith(".out"):
                file_path = os.path.join(folder_path, filename)
                with open(file_path, 'r') as file:
                    data = {}
                    for line in file:
                        terms = [
                            "consensus sequence output",
                            "Number of input sequences considered",
                            "Length of alignment",
                            "sequence found in vulgar file",
                            "number of aligned reads",
                            "Coverage minimum",
                            "Coverage maximum",
                            "Coverage mean",
                            "Coverage median",
                            "# skipped reads due to low rel. score"
                        ]

                        for term in terms:
                            if term in line:
                                value = line.split(":")[1].strip()

                                if value.startswith('"') and value.endswith('"'):
                                    value = value[1:-1]

                                if term == "consensus sequence output":
                                    start_index = value.find("BSNHM")
                                    if start_index != -1:
                                        value = value[start_index:start_index + 11]
                                        data[term] = value
                                        break
                                elif term == "# skipped reads due to low rel. score":
                                    try:
                                        data[term] = int(value)
                                    except ValueError:
                                        data[term] = value
                                else:
                                    try:
                                        numeric_value = float(value)
                                        if numeric_value.is_integer():
                                            data[term] = int(numeric_value)
                                        else:
                                            data[term] = numeric_value
                                    except ValueError:
                                        data[term] = value

                    if data:
                        data_list.append(data)

        if data_list:
            summary_df = pd.DataFrame(data_list)
            print("Summary DataFrame before dropping rows:")
            print(summary_df)

            if "Number of input sequences considered" in summary_df.columns:
                summary_df.dropna(subset=["Number of input sequences considered"], inplace=True)
            else:
                print("Warning: 'Number of input sequences considered' not found in columns.")

            for column in summary_df.select_dtypes(include=['float']).columns:
                if all(summary_df[column].dropna().apply(float.is_integer)):
                    summary_df[column] = summary_df[column].astype('Int64')

            summary_df.to_csv(output_csv, index=False)
            print(f"Cleaned summary stats saved to {output_csv}")
        else:
            print("No data extracted. The output CSV will not be created.")

    except FileNotFoundError:
        print(f"Folder '{folder_path}' not found.")




def process_fasta_file(fasta_file):
    sequences = {}
    try:
        with open(fasta_file, 'r') as file:
            header = None
            sequence = []
            for line in file:
                if line.startswith(">"):
                    if header and sequence:
                        sequences[header] = "".join(sequence)
                    elif header:
                        sequences[header] = ""
                    header = line.strip().split(" ")[0][1:].split("_")[0]
                    sequence = []
                else:
                    sequence.append(line.strip())
            if header and sequence:
                sequences[header] = "".join(sequence)
            elif header:
                sequences[header] = ""
    except FileNotFoundError:
        print(f"FASTA file '{fasta_file}' not found.")
    return sequences





def add_fasta_info_to_df(df, fasta_sequences, consensus_column):
    lengths = []
    non_gtac_counts = []
    n_counts = []
    dash_tilde_counts = []

    for _, row in df.iterrows():
        consensus_id = row.get(consensus_column, "")
        sequence = fasta_sequences.get(consensus_id, "")

        if not sequence or sequence.strip() == "":
            print(f"No valid sequence found for consensus ID: {consensus_id}")
            sequence_length = 0
            non_gtac_count = 0
            n_count = 0
            dash_tilde_count = 0
        else:
            sequence_length = sum(1 for char in sequence if char in 'GTAC')
            non_gtac_count = sum(1 for char in sequence if char not in 'GTAC')
            n_count = sequence.count('N')
            dash_tilde_count = sequence.count('-') + sequence.count('~')

        lengths.append(sequence_length)
        non_gtac_counts.append(non_gtac_count)
        n_counts.append(n_count)
        dash_tilde_counts.append(dash_tilde_count)

    df["Sequence Length"] = lengths
    df["N Count"] = n_counts
    df["Dash and Tilde Count"] = dash_tilde_counts
    df["Non-GTAC Count"] = non_gtac_counts

    return df




if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("""
        Usage: python 5_mge_stats.py /path/to/.out /path/to/[summary_stats].csv /path/to/multi-fasta

        Extract data from MGE .out files and MGE concatenated cox1 consensus sequences in multi-FASTA file to the output .csv
        
        MGE output:
        Length of alignment = length of subject/reference sequence 
        Sequence found in vulgar file = number of sequences used by exonerate
        # skipped reads due to low rel. score = number of reads omitted by MGE due to a low relative exonerate alignment score (-r)
        """)
        sys.exit(1)




    folder_path = sys.argv[1]
    output_csv = sys.argv[2]
    fasta_file = sys.argv[3]

    extract_data_from_out_files(folder_path, output_csv)

    fasta_sequences = process_fasta_file(fasta_file)

    if os.path.exists(output_csv) and os.path.getsize(output_csv) > 0:
        summary_df = pd.read_csv(output_csv)

        summary_df = add_fasta_info_to_df(summary_df, fasta_sequences, "consensus sequence output")

        summary_df.to_csv(output_csv, index=False)
        print(f"Summary stats saved to {output_csv}")
    else:
        print(f"Error: No data extracted. The output CSV '{output_csv}' is empty or does not exist.")
