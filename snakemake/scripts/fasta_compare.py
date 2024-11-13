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


###This script analyzes FASTA files containing DNA sequences, evaluates sequence quality based on various metrics, 
###and selects the best sequences based on specific quality criteria. It processes both full sequences and barcode regions (positions 40-700), 
###outputting filtered sequences to separate FASTA files and providing detailed analysis in a CSV file.

##OUTPUT_CSV: Path where the analysis results CSV file will be saved
##OUTPUT_FASTA: Path where the best full sequences FASTA file will be saved
##OUTPUT_BARCODE_FASTA: Path where the best barcode sequences FASTA file will be saved
##INPUT_FILES: One or more input FASTA files to analyze (space-separated)
##Optional:
#--log-file LOG_FILE: Specify a custom path for the log file (default: creates timestamped log in current directory)
#--verbose, -v: Enable detailed debug logging

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




def analyse_fasta(file_path):
    """
    Analyse a FASTA file and extract various metrics from each sequence.

    Parameters:
        file_path (str): The path to the FASTA file to be analysed.

    Returns:
        dict: A dictionary containing analysis results for each sequence in the FASTA file.
    """
    try:
        #Initialise dictionaries
        results = {}
        sequences = {}
        
        #Count total sequences for progress reporting
        total_sequences = sum(1 for _ in SeqIO.parse(file_path, "fasta"))
        processed = 0

        #Parse the FASTA file
        for record in SeqIO.parse(file_path, "fasta"):
            try:
                processed += 1
                if processed % 100 == 0:  #Log progress every 100 sequences
                    logging.info(f"Processing {file_path}: {processed}/{total_sequences} sequences")

                #Get the sequence ID (header) and the sequence itself
                seq_id = record.id
                seq = str(record.seq)
                
                #Store the complete record
                sequences[seq_id] = record

                #Clean and split seq_id into process_id and parameters
                seq_id_clean = seq_id.replace("Consensus_", "") if seq_id.startswith("Consensus_") else seq_id

                #Initialise process_id and parameters
                process_id = seq_id_clean
                parameters = ""

                #Extract process_id and parameters
                if '_r_' in seq_id_clean:
                    param_start_idx = seq_id_clean.index('_r_')
                    process_id = seq_id_clean[:param_start_idx]
                    parameters = seq_id_clean[param_start_idx + 1:]

                logging.debug(f"Detected process_id: {process_id}, parameters: {parameters}")

                #Calculate sequence metrics
                length = len(seq)
                leading_gaps = len(seq) - len(seq.lstrip('-'))
                trailing_gaps = len(seq) - len(seq.rstrip('-'))
                internal_gaps = seq.count('-') + seq.count('~') - leading_gaps - trailing_gaps
                ambiguous_bases = seq.count('N')

                #Calculate longest stretch without gaps
                longest_stretch = max(len(s) for s in re.split('-|~', seq))

                #Analyse barcode region (positions 40-700)
                subseq = seq[39:700]
                barcode_length = len(subseq)
                barcode_ambiguous_bases = subseq.count('N')
                barcode_longest_stretch = max(len(s) for s in re.split('-|~', subseq))

                #Determine ranks
                barcode_rank = calculate_barcode_rank(barcode_ambiguous_bases, barcode_longest_stretch)
                full_rank = calculate_full_rank(ambiguous_bases)

                #Store results
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
        logging.error(f"Error analyzing file {file_path}: {str(e)}")
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
        #Trim leading and trailing gaps
        sequence = sequence.strip('-').strip('~')
    
    if convert_internal_gaps:
        #First remove ~ characters and stitch sequence
        sequence = sequence.replace('~', '')
        #Then replace remaining gaps (-) with N
        sequence = re.sub(r'-', 'N', sequence)
    
    #Create new sequence record with formatted sequence
    new_record = SeqRecord(
        Seq(sequence),
        id=seq_record.id,
        description=seq_record.description
    )
    
    return new_record





def format_barcode_sequence(seq_record):
    """
    Extract and format the barcode region (positions 40-700).
    
    Parameters:
        seq_record (SeqRecord): The sequence record to process
        
    Returns:
        SeqRecord: Formatted barcode sequence record
    """
    #Extract barcode region (positions 40-700, 0-based indexing)
    sequence = str(seq_record.seq)[39:700]
    
    #Format gaps according to rules
    formatted_sequence = format_barcode_gaps(sequence)
    
    #Since we're keeping the longest fragment, no need to trim N's
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
    - For gaps marked with -, if ≤ 6 bases fill with N
    - For gaps > 6 bases marked with -, keep the longest fragment
    
    Parameters:
        sequence (str): The sequence to process
        
    Returns:
        str: Formatted sequence containing the longest fragment with ~ removed
    """
    #First handle ~ characters by removing them and stitching sequence
    sequence = sequence.replace('~', '')
    
    #Now handle remaining gaps (-). Split sequence by gaps of length > 6
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
                #Fill small gaps with N
                current_fragment.extend(['N'] * gap_count)
            gap_count = 0
            current_fragment.append(char)
    
    #Add last fragment if exists
    if current_fragment:
        fragments.append(''.join(current_fragment))
    
    #Find and get longest fragment
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
    return f"{process_id}|{parameters}" if parameters else process_id





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




def write_best_sequences(best_sequences, output_fasta, output_barcode_fasta):
    """
    Write the best sequences that meet specific criteria to new FASTA files.
    
    Parameters:
        best_sequences (dict): Dictionary containing the best sequences for each process_id
        output_fasta (str): Path to output FASTA file for full sequences
        output_barcode_fasta (str): Path to output FASTA file for barcode regions
    """
    try:
        selected_full_records = []
        selected_barcode_records = []
        selection_records = {}
        
        for process_id, sequence_data in best_sequences.items():
            #Get all results for this process_id
            results = sequence_data[1]  #Access the list of results directly
            
            #Select best sequences using criteria
            best_full, best_barcode = select_sequences_for_process(results)
            
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
                        'full_selected_seq': best_full['seq_id'],
                        'barcode_selected_seq': best_barcode['seq_id'] if best_barcode else None
                    }
                
                #Process barcode sequence, if we have one
                if best_barcode:
                    seq_record = best_barcode['sequence_record']
                    new_id = format_sequence_id(best_barcode['process_id'], best_barcode['parameters'])
                    
                    barcode_seq_record = format_barcode_sequence(seq_record)
                    barcode_seq_record.id = new_id
                    barcode_seq_record.description = ""
                    selected_barcode_records.append(barcode_seq_record)
                    
                    if process_id not in selection_records:
                        selection_records[process_id] = {
                            'full_selected_seq': None,
                            'barcode_selected_seq': best_barcode['seq_id']
                        }
                
                logging.info(f"Process {process_id}:")
                if best_full:
                    logging.info(f"  Selected full sequence {best_full['seq_id']} (full_rank: {best_full['full_rank']}, barcode_rank: {best_full['barcode_rank']})")
                if best_barcode and best_barcode != best_full:
                    logging.info(f"  Selected barcode sequence {best_barcode['seq_id']} (full_rank: {best_barcode['full_rank']}, barcode_rank: {best_barcode['barcode_rank']})")
        
        #Write selected sequences to output files
        if selected_full_records:
            SeqIO.write(selected_full_records, output_fasta, "fasta")
            logging.info(f"Wrote {len(selected_full_records)} sequences to {output_fasta}")
        else:
            logging.warning("No sequences met the criteria for full sequence output")
            
        if selected_barcode_records:
            SeqIO.write(selected_barcode_records, output_barcode_fasta, "fasta")
            logging.info(f"Wrote {len(selected_barcode_records)} sequences to {output_barcode_fasta}")
        else:
            logging.warning("No sequences met the criteria for barcode sequence output")
            
        return selection_records
            
    except Exception as e:
        logging.error(f"Error writing sequences: {str(e)}")
        raise



def main():
    #Arg parser
    parser = argparse.ArgumentParser(description='analyse FASTA files and select best sequences based on quality criteria.')
    parser.add_argument('output_csv', help='Path to output CSV file')
    parser.add_argument('output_fasta', help='Path to output FASTA file for best sequences')
    parser.add_argument('output_barcode_fasta', help='Path to output FASTA file for barcode regions')
    parser.add_argument('input_files', nargs='+', help='Input FASTA files to analyse')
    parser.add_argument('--log-file', help='Path to log file (optional)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    #Setup logging
    setup_logging(args.log_file)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    #Validate input and output files
    if not validate_files(args.input_files, args.output_csv, args.output_fasta):
        sys.exit(1)
    
    try:
        #Initialise results list
        all_results = []
        
        #Analyse each FASTA file
        for file in args.input_files:
            logging.info(f"Processing file: {file}")
            results = analyse_fasta(file)
            for seq_id, result in results.items():
                result['seq_id'] = seq_id
                all_results.append(result)
        
        #Select best sequences
        best_sequences = {}
        for result in all_results:
            process_id = result['process_id']
            criteria = (
                result['barcode_rank'],
                result['full_rank'],
                -result['barcode_longest_stretch'],
                -result['longest_stretch']
            )
            
            if process_id not in best_sequences or criteria < best_sequences[process_id][0]:
                best_sequences[process_id] = [criteria, [result]]
            elif criteria == best_sequences[process_id][0]:
                best_sequences[process_id][1].append(result)
        
        #Annotate best sequences
        for result in all_results:
            process_id = result['process_id']
            result['best_sequence'] = 'yes' if (
                process_id in best_sequences and 
                result in best_sequences[process_id][1]
            ) else 'no'
        
        #Write best sequences to FASTA files and get selection records
        selection_records = write_best_sequences(best_sequences, args.output_fasta, args.output_barcode_fasta)

        #Update results with new columns
        for result in all_results:
            process_id = result['process_id']
            seq_id = result['seq_id']
            
            if process_id in selection_records:
                result['selected_full_fasta'] = 'yes' if selection_records[process_id]['full_selected_seq'] == seq_id else 'no'
                result['selected_barcode_fasta'] = 'yes' if selection_records[process_id]['barcode_selected_seq'] == seq_id else 'no'
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
