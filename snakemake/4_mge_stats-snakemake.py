import os
import sys
import pandas as pd
import time
from Bio.Blast import NCBIWWW, NCBIXML
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
    blast_results = []

    for _, row in df.iterrows():
        consensus_id = row.get(consensus_column, "")
        sequence = fasta_sequences.get(consensus_id, "")

        if not sequence or sequence.strip() == "":
            print(f"No valid sequence found for consensus ID: {consensus_id}")
            sequence_length = 0
            non_gtac_count = 0
            n_count = 0
            dash_tilde_count = 0
            blast_result = {
                'Title': np.nan,
                'Score': np.nan,
                'Alignment Length': np.nan,
                'Identity': np.nan,
                'Percent Identity': np.nan,
                'Gaps': np.nan,
                'Query Coverage': np.nan,
                'Subject Coverage': np.nan,
                'Accession': np.nan
            }
        else:
            sequence_length = sum(1 for char in sequence if char in 'GTAC')
            non_gtac_count = sum(1 for char in sequence if char not in 'GTAC')
            n_count = sequence.count('N')
            dash_tilde_count = sequence.count('-') + sequence.count('~')
            blast_result = run_blast(sequence)

        lengths.append(sequence_length)
        non_gtac_counts.append(non_gtac_count)
        n_counts.append(n_count)
        dash_tilde_counts.append(dash_tilde_count)
        blast_results.append(blast_result)

    df["Sequence Length"] = lengths
    df["N Count"] = n_counts
    df["Dash and Tilde Count"] = dash_tilde_counts
    df["Non-GTAC Count"] = non_gtac_counts
    df["BLAST Title"] = [result.get('Title', np.nan) for result in blast_results]
    df["BLAST Score"] = [result.get('Score', np.nan) for result in blast_results]
    df["BLAST Alignment Length"] = [result.get('Alignment Length', np.nan) for result in blast_results]
    df["BLAST Identity"] = [result.get('Identity', np.nan) for result in blast_results]
    df["BLAST Percent Identity"] = [result.get('Percent Identity', np.nan) for result in blast_results]
    df["BLAST Gaps"] = [result.get('Gaps', np.nan) for result in blast_results]
    df["BLAST Query Coverage"] = [result.get('Query Coverage', np.nan) for result in blast_results]
    df["BLAST Subject Coverage"] = [result.get('Subject Coverage', np.nan) for result in blast_results]
    df["BLAST Accession"] = [result.get('Accession', np.nan) for result in blast_results]

    return df


def run_blast(sequence):
    if not sequence or sequence.strip() == "":
        print("Error: Attempting to run BLAST with an empty sequence.")
        return {
            'Title': np.nan,
            'Score': np.nan,
            'Alignment Length': np.nan,
            'Identity': np.nan,
            'Percent Identity': np.nan,
            'Gaps': np.nan,
            'Query Coverage': np.nan,
            'Subject Coverage': np.nan,
            'Accession': np.nan
        }

    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
        blast_records = NCBIXML.parse(result_handle)

        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    percent_identity = (hsp.identities / hsp.align_length) * 100
                    return {
                        'Title': alignment.title,
                        'Score': hsp.score,
                        'Alignment Length': hsp.align_length,
                        'Identity': hsp.identities,
                        'Percent Identity': percent_identity,
                        'Gaps': hsp.gaps,
                        'Query Coverage': f"{hsp.query_start}-{hsp.query_end}",
                        'Subject Coverage': f"{hsp.sbjct_start}-{hsp.sbjct_end}",
                        'Accession': alignment.accession
                    }

        return {
            'Title': np.nan,
            'Score': np.nan,
            'Alignment Length': np.nan,
            'Identity': np.nan,
            'Percent Identity': np.nan,
            'Gaps': np.nan,
            'Query Coverage': np.nan,
            'Subject Coverage': np.nan,
            'Accession': np.nan
        }

    except Exception as e:
        print(f"Error running BLAST: {e}")
        return {
            'Title': np.nan,
            'Score': np.nan,
            'Alignment Length': np.nan,
            'Identity': np.nan,
            'Percent Identity': np.nan,
            'Gaps': np.nan,
            'Query Coverage': np.nan,
            'Subject Coverage': np.nan,
            'Accession': np.nan
        }

    finally:
#Sleep to handle rate limiting
        time.sleep(2)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("""
        Usage: python 5_mge_stats.py /path/to/.out /path/to/[summary_stats].csv /path/to/multi.fasta

        Extract data from MGE .out files, MGE concatenated cox1 consensus sequences in multi-FASTA file, and BLAST results to the output .csv
        
        MGE output:
        Length of alignment = length of subject/reference sequence 
        Sequence found in vulgar file = number of sequences used by exonerate
        # skipped reads due to low rel. score = number of reads omitted by MGE due to a low relative exonerate alignment score (-r)

        BLAST output:
        Title = title of BLAST alignment
        Score = HSP score assigned to the alignment (higher = better)
        Expect = E-value (expectation value), which indicates the likelihood of the match occurring by chance (lower = better)
        Identity = number of identical matches between the query and subject sequences
        Query Coverage = range of positions in query (input) sequence covered by alignment
        Subject Coverage = range of positions in the subject sequence covered by alignment
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
        print(f"Updated summary stats saved to {output_csv}")
    else:
        print(f"Error: No valid data extracted. The output CSV '{output_csv}' is empty or does not exist.")
