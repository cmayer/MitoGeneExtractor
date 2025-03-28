# =============================================================================
# Imports and Setup
# =============================================================================
import os
import sys
import re
import csv
import argparse
from datetime import datetime
from typing import List, Dict, Tuple, Optional, Any, NamedTuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Align import AlignInfo
import numpy as np
from collections import Counter, defaultdict
import traceback
import multiprocessing as mp
from functools import partial
import concurrent.futures
from itertools import islice



"""
FASTA Sequence Cleaner and Analysis Tool

A comprehensive tool for cleaning and analysing FASTA sequence alignments using multiple filtering approaches:
- Human COX1 similarity detection
- AT content analysis
- Statistical outlier detection
- Reference sequence comparison

The tool processes each alignment file through a sequential filtering pipeline:
1. Remove sequences with high human COX1 similarity
2. Filter sequences with divergent AT content
3. Remove statistical outliers
4. Compare against reference sequences (if provided)

Features:
---------
- Parallel processing with progress tracking and ETA
- Robust error handling and logging
- Multiple AT content filtering modes (absolute, higher, lower)
- Statistical outlier detection with adjustable percentile thresholds
- Optional reference sequence comparison
- Comprehensive metrics calculation and reporting
- Ordered and annotated sequence output

Input:
------
Required:
    - Directory containing FASTA alignment files (.fasta, .fas, or .fa)
    - Output directory for results

Optional:
    - Directory containing reference sequences (named as original_filename_reference.fasta)
    - Number of threads for parallel processing
    - Various filtering thresholds and parameters

Parameters:
-----------
    --input_dir (-i):          Input directory containing FASTA files
    --output_dir (-o):         Output directory for results
    --threads (-t):            Number of threads for parallel processing [default: system CPU count]
    --consensus_threshold (-c): Threshold for consensus generation [default: 0.5]
    --human_threshold (-u):    Human COX1 similarity threshold [default: 0.95]
    --at_difference (-d):      Maximum allowed AT content difference [default: 0.1]
    --at_mode (-m):           AT content filtering mode (absolute/higher/lower) [default: absolute]
    --percentile_threshold (-p): Percentile for outlier detection [default: 90.0]
    --reference_dir (-r):      Optional directory containing reference sequences
    --disable_human:          Disable human COX1 similarity filtering
    --disable_at:            Disable AT content filtering
    --disable_outliers:      Disable statistical outlier detection

Output:
-------
For each input file, generates files organised in subdirectories:
    filter_pass_seqs/
        - {basename}_cleaned.fasta: Filtered sequences that passed all criteria
    consensus_seqs/
        - {basename}_consensus.fasta: Consensus sequence generated from cleaned alignment
    metrics/
        - {basename}_metrics.csv: Detailed metrics for all sequences
    filter_annotated_seqs/
        - {basename}_ordered_annotated.fasta: All sequences with filtering annotations
    logs/
        - {basename}_log.txt: Detailed processing log for each file
        - processing_summary.txt: Overall processing summary

Progress Tracking:
----------------
- Real-time progress bar with percentage completion
- Estimated time remaining (ETA)
- Error reporting for failed files (graceful handling of empty alignments)
- Final summary with total processing time and error count

Authors: B. Price & D. Parsons @ NHMUK
Version: 1.0.2
Licence: MIT
"""


# =============================================================================
# Core Constants
# =============================================================================
# Reference human COX1 sequence
HUMAN_COX1 = """ATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCG
CATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCT
TCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTA
ATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATGGCGT
TTCCCCGCATAAACAACATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTAT
AGTGGAGGCCGGAGCAGGAACAGGTTGAACAGTCTACCCTCCCTTAGCAGGGAACTACTCCCACCCTGGA
GCCTCCGTAGACCTAACCATCTTCTCCTTACACCTAGCAGGTGTCTCCTCTATCTTAGGGGCCATCAATT
TCATCACAACAATTATCAATATAAAACCCCCTGCCATAACCCAATACCAAACGCCCCTCTTCGTCTGATC
CGTCCTAATCACAGCAGTCCTACTTCTCCTATCTCTCCCAGTCCTAGCTGCTGGCATCACTATACTACTA
ACAGACCGCAACCTCAACACCACCTTCTTCGACCCCGCCGGAGGAGGAGACCCCATTCTATACCAACACC
TATTCTGATTTTTCGGTCACCCTGAAGTTTATATTCTTATCCTACCAGGCTTCGGAATAATCTCCCATAT
TGTAACTTACTACTCCGGAAAAAAAGAACCATTTGGATACATAGGTATGGTCTGAGCTATGATATCAATT
GGCTTCCTAGGGTTTATCGTGTGAGCACACCATATATTTACAGTAGGAATAGACGTAGACACACGAGCAT
ATTTCACCTCCGCTACCATAATCATCGCTATCCCCACCGGCGTCAAAGTATTTAGCTGACTCGCCACACT
CCACGGAAGCAATATGAAATGATCTGCTGCAGTGCTCTGAGCCCTAGGATTCATCTTTCTTTTCACCGTA
GGTGGCCTGACTGGCATTGTATTAGCAAACTCATCACTAGACATCGTACTACACGACACGTACTACGTTG
TAGCCCACTTCCACTATGTCCTATCAATAGGAGCTGTATTTGCCATCATAGGAGGCTTCATTCACTGATT
TCCCCTATTCTCAGGCTACACCCTAGACCAAACCTACGCCAAAATCCATTTCACTATCATATTCATCGGC
GTAAATCTAACTTTCTTCCCACAACACTTTCTCGGCCTATCCGGAATGCCCCGACGTTACTCGGACTACC
CCGATGCATACACCACATGAAACATCCTATCATCTGTAGGCTCATTCATTTCTCTAACAGCAGTAATATT
AATAATTTTCATGATTTGAGAAGCCTTCGCTTCGAAGCGAAAAGTCCTAATAGTAGAAGAACCCTCCATA
AACCTGGAGTGACTATATGGATGCCCCCCACCCTACCACACATTCGAAGAACCCGTATACATAAAATCTA
GA""".replace("\n", "")

# Custom type definitions
SequenceRecord = SeqRecord
SequenceAlignment = MultipleSeqAlignment





# =============================================================================
# Utility Functions
# =============================================================================
def get_base_filename(filepath: str) -> str:
    """
    Extract base filename without extension and without '_align' substring.
    """
    basename = os.path.basename(filepath)
    # First remove the file extension
    for ext in ['.fasta', '.fas', '.fa']:
        if basename.lower().endswith(ext):
            basename = basename[:-len(ext)]
            break
    
    # Then remove '_align' if present
    if '_align' in basename:
        basename = basename.replace('_align', '')
    
    return basename

def get_default_threads() -> int:
    """Get default number of threads from SLURM environment or system CPU count."""
    # First check SLURM environment
    slurm_cpus = os.environ.get('SLURM_CPUS_PER_TASK')
    if slurm_cpus:
        try:
            return int(slurm_cpus)
        except ValueError:
            pass
    
    # Fallback to system CPU count
    return mp.cpu_count()

def chunk_list(lst: list, chunk_size: int):
    """
    Yield successive chunks from list.
    """
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]

def create_output_subdirectories(base_output_dir: str) -> Dict[str, str]:
    """
    Create subdirectories for different output file types.
    Omits directories for removed sequences per user request.
    """
    subdirs = {
        'filter_pass_seqs': os.path.join(base_output_dir, 'filter_pass_seqs'),
        'consensus_seqs': os.path.join(base_output_dir, 'consensus_seqs'),
        'metrics': os.path.join(base_output_dir, 'metrics'),
        'filter_annotated_seqs': os.path.join(base_output_dir, 'filter_annotated_seqs'),
        'logs': os.path.join(base_output_dir, 'logs')
    }
    
    # Create all subdirectories
    for subdir in subdirs.values():
        os.makedirs(subdir, exist_ok=True)
    
    return subdirs

def log_message(message: str, log_file=None, stdout=False, error=False) -> None:
    """
    Log message to file and optionally to stdout.
    """
    # Write to log file if provided
    if log_file:
        log_file.write(message + '\n')
        log_file.flush()
    
    # Print to stdout if requested or if it's an error
    if stdout or error:
        if error:
            print(f"\nERROR: {message}", flush=True)
        else:
            print(message, flush=True)

def check_file_already_processed(base_name: str, output_subdirs: Dict[str, str]) -> bool:
    """
    Check if a FASTA file has already been fully processed by checking if expected output files exist.
    """
    # Define expected output files based on base_name
    expected_files = {
        'cleaned': os.path.join(output_subdirs['filter_pass_seqs'], f"{base_name}_cleaned.fasta"),
        'consensus': os.path.join(output_subdirs['consensus_seqs'], f"{base_name}_.fasta"),
        'metrics': os.path.join(output_subdirs['metrics'], f"{base_name}_metrics.csv"),
        'ordered_annotated': os.path.join(output_subdirs['filter_annotated_seqs'], f"{base_name}_ordered_annotated.fasta"),
        'log': os.path.join(output_subdirs['logs'], f"{base_name}_log.txt")
    }
    
    # Check if all files exist
    for file_path in expected_files.values():
        if not os.path.exists(file_path):
            return False
        
        # Additional check: make sure files are not zero-sized (possibly corrupted)
        if os.path.getsize(file_path) == 0:
            return False
    
    return True

# =============================================================================
# Sequence Analysis Functions
# =============================================================================
def sort_sequences_by_position(records: List[SeqRecord]) -> List[SeqRecord]:
    """
    Sort sequence records by their starting position in the alignment.
    """
    # Create tuples of (start_position, original_index, record) to maintain stable sort
    indexed_records = [
        (get_sequence_start_position(str(record.seq)), i, record)
        for i, record in enumerate(records)
    ]
    
    # Sort by start position, using original index as tiebreaker
    sorted_records = [
        record for _, _, record in sorted(indexed_records)
    ]
    
    return sorted_records

def find_best_alignment(query: str, reference: str, min_overlap: int = 20) -> tuple[int, int, int]:
    """
    Find the best local alignment between a query sequence and a reference sequence.
     """
    query_len = len(query)
    ref_len = len(reference)
    
    best_start = 0
    best_end = 0
    max_matches = 0
    
    # Try all possible alignments of query against reference
    for start in range(ref_len - min_overlap + 1):
        # Don't look past end of reference
        end = min(start + query_len, ref_len)
        
        # Count matches in this alignment window
        ref_segment = reference[start:end]
        query_segment = query[:end-start]
        
        matches = sum(1 for q, r in zip(query_segment, ref_segment) 
                     if q == r and q != '-' and r != '-')
        
        # Update best alignment if we found more matches
        if matches > max_matches:
            max_matches = matches
            best_start = start
            best_end = end

    return best_start, best_end, max_matches

def get_sequence_start_position(sequence: str) -> int:
    """Find the position of the first non-gap character in a sequence."""
    for i, char in enumerate(sequence):
        if char != '-':
            return i
    return len(sequence)

def calculate_sequence_similarity(seq1: str, seq2: str, min_overlap: int = 20) -> float:
    """Calculate sequence similarity between two sequences."""
    # Remove gaps and convert to uppercase
    query = seq1.replace('-', '').upper()
    reference = seq2.replace('-', '').upper()
    
    # Check sequence lengths
    if len(query) < min_overlap or len(reference) < min_overlap:
        return 0.0
        
    # Find best local alignment
    start, end, matches = find_best_alignment(query, reference, min_overlap)
    
    # Calculate similarity over aligned region
    aligned_length = end - start
    if aligned_length < min_overlap:
        return 0.0
        
    return matches / aligned_length

def calculate_at_content(sequence: str) -> float:
    """Calculate AT content for a sequence."""
    sequence_no_gaps = sequence.replace('-', '').upper()
    
    if not sequence_no_gaps:
        return 0.0
    
    at_count = sequence_no_gaps.count('A') + sequence_no_gaps.count('T')
    return at_count / len(sequence_no_gaps)

def compare_at_content(consensus_seq: str, query_seq: str) -> Tuple[float, float, float]:
    """Compare AT content between consensus and query sequence."""
    consensus_seq = consensus_seq.upper()
    query_seq = query_seq.upper()
    
    # Find positions where both sequences have bases
    overlapping_positions = [
        (c, q) for c, q in zip(consensus_seq, query_seq)
        if c != '-' and q != '-'
    ]
    
    if not overlapping_positions:
        return 0.0, 0.0, 0.0
    
    # Split into consensus and query sequences
    consensus_bases = ''.join(c for c, _ in overlapping_positions)
    query_bases = ''.join(q for _, q in overlapping_positions)
    
    # Calculate AT content for overlapping regions
    consensus_at = calculate_at_content(consensus_bases)
    query_at = calculate_at_content(query_bases)
    
    return query_at, consensus_at, abs(query_at - consensus_at)

def calculate_position_frequencies(sequences: List[str]) -> List[Dict[str, float]]:
    """Calculate residue frequencies at each position."""
    if not sequences:
        raise ValueError("No sequences provided")
    
    align_length = len(sequences[0])
    position_freqs = []
    
    for i in range(align_length):
        # Get residues at this position, excluding gaps
        residues = [seq[i] for seq in sequences if seq[i] != '-']
        
        if residues:
            counts = Counter(residues)
            total = sum(counts.values())
            frequencies = {res: count/total for res, count in counts.items()}
        else:
            frequencies = {}
            
        position_freqs.append(frequencies)
    
    return position_freqs

def calculate_weighted_deviation(sequence: str, reference: str, 
                              frequencies: List[Dict[str, float]]) -> float:
    """Calculate weighted deviation score."""
    if len(sequence) != len(reference) != len(frequencies):
        raise ValueError("Sequence, reference, and frequencies must all have same length")
    
    total_score = 0.0
    total_weight = 0.0
    
    for seq_res, ref_res, pos_freqs in zip(sequence, reference, frequencies):
        if seq_res != '-' and ref_res != '-':
            # Weight is based on reference residue conservation
            conservation_weight = pos_freqs.get(ref_res, 0)
            total_weight += conservation_weight
            
            if seq_res != ref_res:
                total_score += conservation_weight
    
    return total_score / total_weight if total_weight > 0 else 0.0

def calculate_unweighted_deviation(sequence: str, reference: str) -> float:
    """Calculate unweighted deviation score."""
    if len(sequence) != len(reference):
        raise ValueError("Sequence and reference must be same length")
    
    differences = 0
    valid_positions = 0
    
    for seq_res, ref_res in zip(sequence, reference):
        if seq_res != '-' and ref_res != '-':
            valid_positions += 1
            if seq_res != ref_res:
                differences += 1
    
    return differences / valid_positions if valid_positions > 0 else 0.0

def analyse_sequence_metrics(record: SeqRecord, 
                           consensus_seq: Optional[str] = None,
                           frequencies: Optional[List[Dict[str, float]]] = None,
                           alignment_length: Optional[int] = None) -> Dict[str, float]:
    """
    Calculate comprehensive sequence metrics.
    """
    sequence = str(record.seq).upper()
    sequence_no_gaps = sequence.replace('-', '')
    
    metrics = {
        'length': len(sequence_no_gaps),
        'at_content': calculate_at_content(sequence_no_gaps),
        'human_similarity': calculate_sequence_similarity(sequence, HUMAN_COX1)
    }
    
    # Add coverage metrics if alignment length is provided
    if alignment_length is not None:
        coverage_metrics = calculate_sequence_coverage(record, alignment_length)
        metrics.update(coverage_metrics)
    
    if consensus_seq is not None:
        # Calculate AT content comparison
        query_at, cons_at, at_diff = compare_at_content(consensus_seq, sequence)
        metrics.update({
            'consensus_at': cons_at,
            'at_difference': at_diff,
            'unweighted_deviation': calculate_unweighted_deviation(sequence, consensus_seq)
        })
        
        # Add weighted deviation if frequencies are provided
        if frequencies is not None:
            metrics['weighted_deviation'] = calculate_weighted_deviation(
                sequence, consensus_seq, frequencies
            )
    
    return metrics

def generate_consensus_sequence(alignment: MultipleSeqAlignment, 
                              threshold: float = 0.5) -> Tuple[str, List[Dict[str, float]]]:
    """
    Generate consensus sequence and calculate position-specific frequencies.
    """
    if not alignment:
        raise ValueError("Empty alignment provided")
    
    # Extract sequences as strings
    sequences = [str(record.seq) for record in alignment]
    
    # Calculate position-specific frequencies
    frequencies = calculate_position_frequencies(sequences)
    
    # Generate consensus sequence
    consensus = []
    for pos_freqs in frequencies:
        if pos_freqs:
            # Get most common residue and its frequency
            most_common = max(pos_freqs.items(), key=lambda x: x[1])
            if most_common[1] >= threshold:
                consensus.append(most_common[0])
            else:
                consensus.append('-')
        else:
            consensus.append('-')
    
    return ''.join(consensus), frequencies

def calculate_coverage_depth(alignment: MultipleSeqAlignment) -> Tuple[np.ndarray, Dict[str, float]]:
    """
    Calculate coverage depth at each position of the alignment.
    """
    if not alignment:
        return np.array([]), {
            'coverage_percent': 0.0,
            'min_coverage': 0.0,
            'max_coverage': 0.0,
            'avg_coverage': 0.0,
            'median_coverage': 0.0
        }
    
    # Get alignment length
    align_length = alignment.get_alignment_length()
    
    # Initialize depth array
    depth_array = np.zeros(align_length, dtype=int)
    
    # Count non-gap characters at each position
    for record in alignment:
        seq_str = str(record.seq).upper()
        for i, char in enumerate(seq_str):
            if char != '-':
                depth_array[i] += 1
    
    # Calculate coverage statistics
    covered_positions = np.count_nonzero(depth_array)
    coverage_percent = (covered_positions / align_length) * 100 if align_length > 0 else 0
    
    # Handle empty alignments
    if covered_positions == 0:
        return depth_array, {
            'coverage_percent': 0.0,
            'min_coverage': 0.0,
            'max_coverage': 0.0,
            'avg_coverage': 0.0,
            'median_coverage': 0.0
        }
    
    # Only consider positions that have coverage
    non_zero_depths = depth_array[depth_array > 0]
    
    coverage_stats = {
        'coverage_percent': coverage_percent,
        'min_coverage': np.min(non_zero_depths),
        'max_coverage': np.max(depth_array),
        'avg_coverage': np.mean(non_zero_depths),
        'median_coverage': np.median(non_zero_depths)
    }
    
    return depth_array, coverage_stats

def calculate_sequence_coverage(record: SeqRecord, alignment_length: int) -> Dict[str, float]:
    """
    Calculate coverage metrics for an individual sequence.
    """
    sequence = str(record.seq).upper()
    
    # Count non-gap positions
    non_gap_count = sum(1 for char in sequence if char != '-')
    
    # Calculate coverage percentage
    coverage_percent = (non_gap_count / alignment_length) * 100 if alignment_length > 0 else 0
    
    return {
        'sequence_coverage_percent': coverage_percent,
        'sequence_length_no_gaps': non_gap_count
    }


# =============================================================================
# Reference Analysis Functions      
# =============================================================================
def get_reference_sequence(reference_file: str) -> str:
    """
    Extract and validate reference sequence from a FASTA file.
    """
    try:
        references = list(SeqIO.parse(reference_file, "fasta"))
        if not references:
            raise ValueError("No sequences found in reference file")
        if len(references) > 1:
            print(f"Warning: Multiple sequences found in {reference_file}, using first one")
        return str(references[0].seq)
    except Exception as e:
        raise ValueError(f"Error reading reference file: {str(e)}")

def calculate_reference_metrics(sequence: str, 
                             reference_seq: str,
                             frequencies: List[Dict[str, float]]) -> Dict[str, float]:
    """
    Calculate comprehensive reference-based metrics for a sequence.
    """
    return {
        'reference_unweighted_deviation': calculate_unweighted_deviation(sequence, reference_seq),
        'reference_weighted_deviation': calculate_weighted_deviation(sequence, reference_seq, frequencies),
        'reference_at_difference': abs(calculate_at_content(sequence) - calculate_at_content(reference_seq))
    }



# =============================================================================
# Filtering Functions
# =============================================================================
def filter_human_sequences(records: List[SeqRecord], 
                         threshold: float = 0.95) -> Tuple[List[SeqRecord], List[SeqRecord]]:
    """Filter out sequences with high human COX1 similarity."""
    kept = []
    removed = []
    
    for record in records:
        similarity = calculate_sequence_similarity(str(record.seq), HUMAN_COX1)
        if similarity >= threshold:
            removed.append(record)
        else:
            kept.append(record)
    
    return kept, removed

def filter_at_content(records: List[SeqRecord],
                     consensus_seq: str,
                     threshold: float = 0.1,
                     mode: str = 'absolute') -> Tuple[List[SeqRecord], List[SeqRecord]]:
    """Filter sequences based on AT content."""
    kept = []
    removed = []
    
    for record in records:
        query_at, cons_at, at_diff = compare_at_content(consensus_seq, str(record.seq))
        # Calculate signed difference (positive means query has higher AT content)
        signed_diff = query_at - cons_at
        
        should_remove = False
        if mode == 'absolute':
            should_remove = abs(signed_diff) > threshold
        elif mode == 'higher':
            should_remove = signed_diff > threshold
        elif mode == 'lower':
            should_remove = signed_diff < -threshold
        
        if should_remove:
            removed.append(record)
        else:
            kept.append(record)
    
    return kept, removed

def filter_statistical_outliers(records: List[SeqRecord],
                              consensus_seq: str,
                              frequencies: List[Dict[str, float]],
                              percentile: float = 90.0) -> Tuple[List[SeqRecord], List[SeqRecord]]:
    """Filter sequences based on statistical outlier detection."""
    # Calculate deviation scores
    unweighted_scores = []
    weighted_scores = []
    
    for record in records:
        sequence = str(record.seq)
        unweighted = calculate_unweighted_deviation(sequence, consensus_seq)
        weighted = calculate_weighted_deviation(sequence, consensus_seq, frequencies)
        
        unweighted_scores.append(unweighted)
        weighted_scores.append(weighted)
    
    # Calculate thresholds (using non-zero scores)
    non_zero_unweighted = [s for s in unweighted_scores if s > 0]
    non_zero_weighted = [s for s in weighted_scores if s > 0]
    
    # Check if we have any variation to analyse
    if not non_zero_unweighted and not non_zero_weighted:
        log_message("Warning: No sequence variation detected for outlier analysis", 
                   log_file=None, stdout=False)  # Log to file only
        return records, []  # Return all sequences as kept, none removed
    
    # Only calculate thresholds if we have non-zero scores
    if non_zero_unweighted:
        unweighted_threshold = np.percentile(non_zero_unweighted, percentile)
    else:
        unweighted_threshold = float('inf')
        
    if non_zero_weighted:
        weighted_threshold = np.percentile(non_zero_weighted, percentile)
    else:
        weighted_threshold = float('inf')
    
    # Filter sequences
    kept = []
    removed = []
    
    for record, unw_score, w_score in zip(records, unweighted_scores, weighted_scores):
        if unw_score > unweighted_threshold or w_score > weighted_threshold:
            removed.append(record)
        else:
            kept.append(record)
    
    return kept, removed

def filter_reference_outliers(records: List[SeqRecord],
                            reference_seq: str,
                            frequencies: List[Dict[str, float]],
                            percentile: float = 90.0) -> Tuple[List[SeqRecord], List[SeqRecord], Dict[str, float]]:
    """Filter sequences based on reference comparison."""
    # Calculate metrics for all sequences
    sequence_metrics = {
        record.id: calculate_reference_metrics(
            str(record.seq), reference_seq, frequencies
        ) for record in records
    }
    
    # Calculate thresholds using non-zero scores
    unweighted_scores = [
        metrics['reference_unweighted_deviation'] 
        for metrics in sequence_metrics.values()
        if metrics['reference_unweighted_deviation'] > 0
    ]
    weighted_scores = [
        metrics['reference_weighted_deviation']
        for metrics in sequence_metrics.values()
        if metrics['reference_weighted_deviation'] > 0
    ]
    
    # Calculate percentile thresholds
    unweighted_threshold = np.percentile(unweighted_scores, percentile)
    weighted_threshold = np.percentile(weighted_scores, percentile)
    
    # Filter sequences
    kept = []
    removed = []
    
    for record in records:
        metrics = sequence_metrics[record.id]
        if (metrics['reference_unweighted_deviation'] > unweighted_threshold or
            metrics['reference_weighted_deviation'] > weighted_threshold):
            removed.append(record)
        else:
            kept.append(record)
    
    # Calculate aggregate metrics
    aggregate_metrics = {
        'mean_unweighted_deviation': np.mean(unweighted_scores),
        'mean_weighted_deviation': np.mean(weighted_scores),
        'unweighted_threshold': unweighted_threshold,
        'weighted_threshold': weighted_threshold
    }
    
    return kept, removed, aggregate_metrics

def count_ambiguous_bases(consensus_sequence: str) -> int:
    """
    Count the number of ambiguous bases (not G, T, A, or C) in a consensus sequence.
    
    Args:
        consensus_sequence: Consensus sequence string
    
    Returns:
        Number of ambiguous bases
    """
    # Convert to uppercase and remove gaps
    sequence = consensus_sequence.upper().replace('-', '')
    
    # Count non-GTAC bases
    non_standard_bases = sum(1 for base in sequence if base not in 'GTAC')
    
    return non_standard_bases


# =============================================================================
# Parallel Processing Functions
# =============================================================================
def parallel_sequence_analysis(records: List[SeqRecord], 
                             consensus_seq: Optional[str],
                             frequencies: Optional[List[Dict[str, float]]],
                             alignment_length: Optional[int] = None,
                             chunk_size: int = 100) -> Dict[str, Dict[str, float]]:
    """
    Parallel version of sequence metric calculation.
    """
    def process_chunk(chunk):
        return {record.id: analyse_sequence_metrics(record, consensus_seq, frequencies) 
                for record in chunk}
    
    # Split records into chunks
    record_chunks = list(chunk_list(records, chunk_size))
    
    # Process chunks in parallel
    with concurrent.futures.ProcessPoolExecutor() as executor:
        chunk_results = executor.map(process_chunk, record_chunks)
    
    # Combine results
    all_metrics = {}
    for chunk_result in chunk_results:
        all_metrics.update(chunk_result)
    
    return all_metrics

def parallel_filter_human_sequences(records: List[SeqRecord], 
                                  threshold: float = 0.95,
                                  chunk_size: int = 100) -> Tuple[List[SeqRecord], List[SeqRecord]]:
    """Parallel version of human sequence filtering."""
    def process_chunk(chunk):
        kept, removed = [], []
        for record in chunk:
            similarity = calculate_sequence_similarity(str(record.seq), HUMAN_COX1)
            (removed if similarity >= threshold else kept).append(record)
        return kept, removed
    
    # Split records into chunks
    record_chunks = list(chunk_list(records, chunk_size))
    
    # Process chunks in parallel
    with concurrent.futures.ProcessPoolExecutor() as executor:
        chunk_results = list(executor.map(process_chunk, record_chunks))
    
    # Combine results
    kept, removed = [], []
    for k, r in chunk_results:
        kept.extend(k)
        removed.extend(r)
    
    return kept, removed

def optimise_sequence_similarity(seq1: str, seq2: str, min_overlap: int = 20) -> float:
    """optimised sequence similarity calculation using numpy."""
    # Remove gaps and convert to uppercase
    query = np.array(list(seq1.replace('-', '').upper()))
    reference = np.array(list(seq2.replace('-', '').upper()))
    
    if len(query) < min_overlap or len(reference) < min_overlap:
        return 0.0
    
    # Create a boolean mask for matches
    matches = query[:, None] == reference[None, :]
    
    # Use sliding window to find best alignment
    max_matches = 0
    best_start = 0
    
    for start in range(len(reference) - min_overlap + 1):
        end = min(start + len(query), len(reference))
        window_matches = np.sum(matches[:end-start, start:end])
        if window_matches > max_matches:
            max_matches = window_matches
            best_start = start
    
    aligned_length = min(len(query), len(reference) - best_start)
    return max_matches / aligned_length if aligned_length >= min_overlap else 0.0




# =============================================================================
# File Writing Functions
# =============================================================================
class FilterResult(NamedTuple):
    """Container for sequence filtering results"""
    kept_records: List[SeqRecord]
    removed_records: Dict[str, List[SeqRecord]]  # Maps reason -> list of removed records
    metrics: Dict[str, Dict[str, float]]  # Maps sequence ID -> metric values
    consensus_seq: str  # Final consensus sequence
    frequencies: List[Dict[str, float]]  # Position-specific residue frequencies
    reference_metrics: Optional[Dict[str, float]] = None  # Optional reference-based metrics

def write_sequences_to_fasta(records: List[SeqRecord], filename: str) -> None:
    """Write sequences to FASTA file."""
    with open(filename, 'w') as handle:
        for record in records:
            handle.write(f">{record.id}\n{str(record.seq)}\n")

def write_metrics_report(filter_result: FilterResult, output_path: str) -> None:
    """Write comprehensive sequence metrics report."""
    with open(output_path, 'w', newline='') as csvfile:
        # Gather all possible metric fields
        metric_fields = set()
        for seq_metrics in filter_result.metrics.values():
            metric_fields.update(seq_metrics.keys())
        
        # Add standard fields
        base_fields = [
            'sequence_id',
            'length',
            'at_content',
            'human_similarity',
            'sequence_coverage_percent',
            'sequence_length_no_gaps',
            'consensus_at',
            'at_difference',
            'unweighted_deviation',
            'weighted_deviation'
        ]
        
        # Add reference fields if available
        reference_fields = []
        if filter_result.reference_metrics:
            reference_fields = [
                'reference_unweighted_deviation',
                'reference_weighted_deviation',
                'reference_at_difference'
            ]
        
        # Position-specific metrics if available
        position_fields = []
        if filter_result.frequencies:
            position_fields = [
                'conservation_score',
                'gap_frequency'
            ]
        
        # Final fields list
        fields = base_fields + reference_fields + position_fields + ['removal_reason']
        
        writer = csv.writer(csvfile)
        writer.writerow(fields)
        
        # Create removal reason mapping
        removal_mapping = {}
        for reason, records in filter_result.removed_records.items():
            for record in records:
                removal_mapping[record.id] = reason
        
        # Write data for each sequence
        for seq_id, seq_metrics in filter_result.metrics.items():
            row = [seq_id]
            
            # Add base metrics
            for field in base_fields[1:]:  # Skip sequence_id
                value = seq_metrics.get(field, '')
                row.append(f"{value:.4f}" if isinstance(value, float) else str(value))
            
            # Add reference metrics
            if reference_fields:
                for field in reference_fields:
                    value = seq_metrics.get(field, '')
                    row.append(f"{value:.4f}" if isinstance(value, float) else str(value))
            
            # Add position-specific metrics
            if position_fields:
                for field in position_fields:
                    value = seq_metrics.get(field, '')
                    row.append(f"{value:.4f}" if isinstance(value, float) else str(value))
            
            # Add removal reason
            row.append(removal_mapping.get(seq_id, 'kept'))
            
            writer.writerow(row)


def write_summary(total_stats: Dict[str, Any], log_dir: str) -> None:
    """Write summary statistics to a file."""
    summary_log_path = os.path.join(log_dir, 'processing_summary.txt')
    with open(summary_log_path, 'w') as summary_file:
        summary_file.write("Processing Summary:\n")
        summary_file.write("===================\n\n")
        summary_file.write(f"Total files processed: {total_stats['processed_files']}\n")
        summary_file.write(f"Total input sequences: {total_stats['total_sequences']}\n")
        summary_file.write(f"Total sequences kept: {total_stats['kept_sequences']}\n")
        summary_file.write("\nRemoved sequences by category:\n")
        for key, value in total_stats.items():
            if key.startswith('removed_'):
                reason = key.replace('removed_', '')
                summary_file.write(f"- {reason}: {value}\n")
    
    print(f"\nDetailed summary written to: {summary_log_path}")

def create_annotated_sequence_record(record: SeqRecord, fate: str) -> SeqRecord:
    """Create a new sequence record with fate annotation."""
    # Create a deep copy to avoid modifying the original
    new_record = SeqRecord(
        seq=record.seq,
        id=f"{fate}_{record.id}",
        name=record.name,
        description=record.description
    )
    return new_record

def write_ordered_annotated_alignment(kept_records: List[SeqRecord],
                                    removed_records: Dict[str, List[SeqRecord]],
                                    output_path: str) -> None:
    """Write all sequences with annotations."""
    # Create annotated records for all sequences
    all_records = []
    
    # Add kept sequences with 'kept' prefix
    all_records.extend(
        create_annotated_sequence_record(record, 'kept')
        for record in kept_records
    )
    
    # Add removed sequences with appropriate prefixes
    for reason, records in removed_records.items():
        prefix = f"removed_{reason}"
        all_records.extend(
            create_annotated_sequence_record(record, prefix)
            for record in records
        )
    
    # Sort all records by start position
    sorted_records = sort_sequences_by_position(all_records)
    
    # Write to file
    write_sequences_to_fasta(sorted_records, output_path)

def concatenate_consensus_sequences(consensus_dir: str, output_file: str, preprocessing_mode: str = 'concat') -> None:
    """
    Concatenate all consensus sequences from the consensus directory into a single FASTA file.
    
    Args:
        consensus_dir: Directory containing individual consensus FASTA files
        output_file: Path to output concatenated FASTA file
        preprocessing_mode: If 'merge', '_merge' will be appended to sequence headers
    """
    consensus_files = [
        f for f in os.listdir(consensus_dir)
        if f.lower().endswith(('_consensus.fasta', '_consensus.fas', '_consensus.fa'))
    ]
    
    if not consensus_files:
        print("No consensus sequences found to concatenate")
        return
    
    # Sort files to ensure consistent order
    consensus_files.sort()
    
    with open(output_file, 'w') as outfile:
        for filename in consensus_files:
            filepath = os.path.join(consensus_dir, filename)
            with open(filepath) as infile:
                for line in infile:
                    # Modify headers
                    if line.startswith('>'):
                        # Replace '_consensus' with '_fcleaner'
                        modified_line = line.replace('_consensus', '_fcleaner')
                        
                        # Append '_merge' if preprocessing mode is 'merge'
                        if preprocessing_mode == 'merge' and not modified_line.rstrip().endswith('_merge'):
                            modified_line = modified_line.rstrip() + '_merge\n'
                        
                        outfile.write(modified_line)
                    else:
                        outfile.write(line)

def extract_log_statistics(log_file_path: str) -> Dict[str, Any]:
    """
    Extract statistics from a sample log file.
    
    Args:
        log_file_path: Path to the log file
    
    Returns:
        Dictionary containing extracted statistics
    """
    stats = {
        'sample_name': '',
        'input_reads': 0,
        'removed_human': 0,
        'removed_at_distance': 0,
        'removed_outliers': 0,
        'cleaned_reads': 0,
        'cov_percent': 0.0,
        'cov_mean': 0.0,
        'cov_max': 0,
        'cov_min': 0
    }
    
    try:
        with open(log_file_path, 'r') as f:
            log_content = f.read()
            
            # Extract sample name from "Processing Sample:" line
            sample_match = re.search(r"Processing Sample: (.+)$", log_content, re.MULTILINE)
            if sample_match:
                stats['sample_name'] = sample_match.group(1)
            
            # Extract input sequences
            input_match = re.search(r"Input sequences: (\d+)", log_content)
            if input_match:
                stats['input_reads'] = int(input_match.group(1))
            
            # Extract removed sequences by category
            human_match = re.search(r"Removed \(human_similar\): (\d+)", log_content)
            if human_match:
                stats['removed_human'] = int(human_match.group(1))
                
            at_match = re.search(r"Removed \(at_difference\): (\d+)", log_content)
            if at_match:
                stats['removed_at_distance'] = int(at_match.group(1))
                
            outlier_match = re.search(r"Removed \(statistical_outlier\): (\d+)", log_content)
            if outlier_match:
                stats['removed_outliers'] = int(outlier_match.group(1))
            
            # Extract kept sequences
            kept_match = re.search(r"Kept sequences: (\d+)", log_content)
            if kept_match:
                stats['cleaned_reads'] = int(kept_match.group(1))
            
            # Extract coverage statistics
            cov_percent_match = re.search(r"Coverage percentage: ([\d.]+)%", log_content)
            if cov_percent_match:
                stats['cov_percent'] = float(cov_percent_match.group(1))
                
            cov_min_match = re.search(r"Minimum coverage depth: (\d+)", log_content)
            if cov_min_match:
                stats['cov_min'] = int(cov_min_match.group(1))
                
            cov_max_match = re.search(r"Maximum coverage depth: (\d+)", log_content)
            if cov_max_match:
                stats['cov_max'] = int(cov_max_match.group(1))
                
            cov_avg_match = re.search(r"Average coverage depth: ([\d.]+)", log_content)
            if cov_avg_match:
                stats['cov_mean'] = float(cov_avg_match.group(1))
    
    except Exception as e:
        print(f"Error extracting statistics from {log_file_path}: {str(e)}")
    
    return stats

def generate_combined_statistics(output_dir: str) -> None:
    """
    Generate a combined statistics CSV file for all processed samples,
    including those that were skipped during the current run.
    """
    
    # Define paths
    logs_dir = os.path.join(output_dir, 'logs')
    consensus_dir = os.path.join(output_dir, 'consensus_seqs')
    output_file = os.path.join(output_dir, 'combined_statistics.csv')
    
    # Check if directories exist
    if not os.path.exists(logs_dir) or not os.path.exists(consensus_dir):
        print(f"Required directories not found in {output_dir}")
        return
    
    # Get list of log files (excludes the main fasta_cleaner.log)
    log_files = [f for f in os.listdir(logs_dir) 
                if f.endswith('_log.txt') and f != 'fasta_cleaner.log']
    
    # Initialize results list
    all_stats = []
    
    # Process each log file
    for log_file in log_files:
        log_path = os.path.join(logs_dir, log_file)
        
        # Extract statistics from log file
        stats = extract_log_statistics(log_path)
        
        # Get sample name from log statistics
        sample_name = stats['sample_name']
        
        # Find corresponding consensus file
        consensus_file = f"{sample_name}_consensus.fasta"
        consensus_path = os.path.join(consensus_dir, consensus_file)
        
        # Count ambiguous bases in consensus sequence if file exists
        ambig_bases = 0
        if os.path.exists(consensus_path):
            try:
                # Read consensus sequence
                consensus_records = list(SeqIO.parse(consensus_path, "fasta"))
                if consensus_records:
                    consensus_seq = str(consensus_records[0].seq)
                    ambig_bases = count_ambiguous_bases(consensus_seq)
            except Exception as e:
                print(f"Error counting ambiguous bases in {consensus_path}: {str(e)}")
        
        # Add ambiguous bases count to statistics
        stats['final_ambig_bases'] = ambig_bases
        
        # Add to results list
        all_stats.append(stats)
    
    # Write results to CSV file
    if all_stats:
        with open(output_file, 'w', newline='') as csvfile:
            fieldnames = [
                'sample_name', 'input_reads', 'removed_human', 'removed_at_distance',
                'removed_outliers', 'cleaned_reads', 'final_ambig_bases',
                'cov_percent', 'cov_mean', 'cov_max', 'cov_min'
            ]
            
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for stats in all_stats:
                writer.writerow(stats)
        
        print(f"Combined statistics written to: {output_file}")
    else:
        print("No statistics to write")



# =============================================================================
# Main Process Control Functions
# =============================================================================
def apply_filters(alignment: MultipleSeqAlignment,
                 consensus_threshold: float = 0.5,
                 human_threshold: float = 0.95,
                 at_threshold: float = 0.1,
                 at_mode: str = 'absolute',
                 outlier_percentile: float = 90.0,
                 reference_seq: Optional[str] = None,
                 enable_human: bool = True,
                 enable_at: bool = True,
                 enable_outliers: bool = True,
                 enable_reference: bool = True) -> FilterResult:
    """
    Apply filtering methods in strict order with consensus recalculation between steps.
    
    The filtering pipeline follows this exact order:
    1. Human COX1 similarity (if enabled)
    2. AT content difference (if enabled)
    3. Statistical outliers (if enabled)
    4. Reference sequence comparison (if enabled and provided)
    
    After each filtering step, if sequences were removed:
    - A new consensus sequence is generated from remaining sequences
    - New position-specific frequencies are calculated
    - New metrics are computed for all remaining sequences
    
    Args:
        alignment: Input sequence alignment
        consensus_threshold: Threshold for consensus generation (0-1)
        human_threshold: Human similarity threshold for removal (0-1)
        at_threshold: Maximum allowed AT content difference (0-1)
        at_mode: AT content filtering mode ('absolute', 'higher', or 'lower')
        outlier_percentile: Percentile threshold for outlier detection (0-100)
        reference_seq: Optional reference sequence string
        enable_human: Enable human similarity filtering
        enable_at: Enable AT content filtering
        enable_outliers: Enable statistical outlier detection
        enable_reference: Enable reference sequence filtering
    
    Returns:
        FilterResult containing:
        - kept_records: Sequences that passed all filters
        - removed_records: Dict mapping removal reasons to removed sequences
        - metrics: Dict mapping sequence IDs to their metric values
        - consensus_seq: Final consensus sequence after all filtering
    """
    # Initialise our tracking structures
    current_records = list(alignment)
    removed_records = {
        'human_similar': [],
        'at_difference': [],
        'statistical_outlier': [],
        'reference_outlier': []
    }
    
    # Get alignment length for coverage calculations
    alignment_length = alignment.get_alignment_length()
    
    # We'll track metrics for ALL sequences, including removed ones
    metrics = {}
    
    # Function to update consensus and frequencies
    def update_consensus(records: List[SeqRecord]) -> Tuple[str, List[Dict[str, float]]]:
        if not records:
            return '', []
        temp_alignment = MultipleSeqAlignment(records)
        consensus, freqs = generate_consensus_sequence(temp_alignment, consensus_threshold)
        return str(consensus), freqs
    
    # Initialise consensus and frequencies - used only for initial metrics
    initial_consensus_seq, initial_frequencies = update_consensus(current_records)

    # Calculate initial metrics for all sequences
    for record in current_records:
        metrics[record.id] = analyse_sequence_metrics(
            record, initial_consensus_seq, initial_frequencies, alignment_length
        )
    
    # 1. Human COX1 similarity filtering
    if enable_human and current_records:
        kept, removed = filter_human_sequences(current_records, human_threshold)
        removed_records['human_similar'].extend(removed)
        current_records = kept
        
        # Always recalculate consensus after human filtering, before AT content comparison
        consensus_seq, frequencies = update_consensus(current_records)
        # Update metrics for remaining sequences
        for record in current_records:
            metrics[record.id].update(
                analyse_sequence_metrics(record, consensus_seq, frequencies, alignment_length)
            )
    
    # 2. AT content filtering with mode support
    if enable_at and current_records:
        kept, removed = filter_at_content(
            current_records,
            consensus_seq,
            threshold=at_threshold,
            mode=at_mode
        )
        removed_records['at_difference'].extend(removed)
        current_records = kept
        
        if removed:
            consensus_seq, frequencies = update_consensus(current_records)
            for record in current_records:
                metrics[record.id].update(
                    analyse_sequence_metrics(record, consensus_seq, frequencies, alignment_length)
                )
    
    # 3. Statistical outlier filtering
    if enable_outliers and current_records:
        kept, removed = filter_statistical_outliers(
            current_records, consensus_seq, frequencies, outlier_percentile
        )
        removed_records['statistical_outlier'].extend(removed)
        current_records = kept
        
        if removed:
            consensus_seq, frequencies = update_consensus(current_records)
            for record in current_records:
                metrics[record.id].update(
                    analyse_sequence_metrics(record, consensus_seq, frequencies, alignment_length)
                )
    
    # 4. Reference sequence filtering (if provided)
    if enable_reference and reference_seq and current_records:
        kept, removed, ref_metrics = filter_reference_outliers(
            current_records, reference_seq, frequencies, outlier_percentile
        )
        removed_records['reference_outlier'].extend(removed)
        current_records = kept
        
        # Note: We don't recalculate consensus after reference filtering
        # as it's our last step and not needed for further filtering
    
    # Calculate reference metrics if reference sequence provided
    reference_metrics = None
    if reference_seq and enable_reference:
        reference_metrics = {
            'mean_deviation': np.mean([
                calculate_unweighted_deviation(str(record.seq), reference_seq)
                for record in current_records
            ]),
            'weighted_deviation': np.mean([
                calculate_weighted_deviation(str(record.seq), reference_seq, frequencies)
                for record in current_records
            ])
        }
    
    return FilterResult(
        kept_records=current_records,
        removed_records=removed_records,
        metrics=metrics,
        consensus_seq=consensus_seq,
        frequencies=frequencies,
        reference_metrics=reference_metrics
    )

def process_fasta_file(input_file: str, 
                      output_dir: str,
                      reference_file: Optional[str] = None,
                      consensus_threshold: float = 0.5,
                      human_threshold: float = 0.95,
                      at_threshold: float = 0.1,
                      at_mode: str = 'absolute',
                      outlier_percentile: float = 90.0,
                      enable_human: bool = True,
                      enable_at: bool = True,
                      enable_outliers: bool = True,
                      enable_reference: bool = True) -> Dict[str, Any]:
    """Process a single FASTA file with all enabled filtering methods."""
    # Get base name for clear identification
    base_name = get_base_filename(input_file)
    original_filename = os.path.basename(input_file)
    
    # Create output subdirectories
    output_subdirs = create_output_subdirectories(output_dir)
    
    # Set up all output paths
    output_paths = {
        'cleaned': os.path.join(output_subdirs['filter_pass_seqs'], f"{base_name}_cleaned.fasta"),
        'consensus': os.path.join(output_subdirs['consensus_seqs'], f"{base_name}_fcleaner.fasta"),
        'metrics': os.path.join(output_subdirs['metrics'], f"{base_name}_metrics.csv"),
        'ordered_annotated': os.path.join(output_subdirs['filter_annotated_seqs'], f"{base_name}_ordered_annotated.fasta"),
        'log': os.path.join(output_subdirs['logs'], f"{base_name}_log.txt")
    }
    
    # Initialise statistics
    stats = {
        'input_sequences': 0,
        'kept_sequences': 0,
        'removed_sequences': defaultdict(int),
        'processing_time': None,
        'filter_result': None
    }
    
    start_time = datetime.now()
    
    with open(output_paths['log'], 'w') as log_file:
        try:
            # Write separator to log file
            separator = "=" * 80
            log_message(separator, log_file)
            log_message(f"Processing Sample: {base_name}", log_file)
            log_message(separator, log_file)
            log_message(f"Started at: {start_time}", log_file)
            log_message(f"Input file: {input_file}\n", log_file)
            
            # Configuration section
            log_message(separator, log_file)
            log_message("Configuration:", log_file)
            log_message(separator, log_file)
            log_message("Parameters:", log_file)
            log_message(f"- Consensus threshold: {consensus_threshold}", log_file)
            log_message(f"- Human threshold: {human_threshold}", log_file)
            log_message(f"- AT threshold: {at_threshold}", log_file)
            log_message(f"- AT mode: {at_mode}", log_file)
            log_message(f"- Outlier percentile: {outlier_percentile}", log_file)
            
            log_message("\nEnabled filters:", log_file)
            log_message(f"- Human similarity: {'Yes' if enable_human else 'No'}", log_file)
            log_message(f"- AT content: {'Yes' if enable_at else 'No'}", log_file)
            log_message(f"- Statistical outliers: {'Yes' if enable_outliers else 'No'}", log_file)
            log_message(f"- Reference comparison: {'Yes' if enable_reference else 'No'}\n", log_file)
            
            # Read input alignment
            alignment = AlignIO.read(input_file, "fasta")
            stats['input_sequences'] = len(alignment)
            log_message(f"Input sequences: {stats['input_sequences']}\n", log_file)

            # Calculate coverage metrics for the entire alignment
            coverage_depth, coverage_stats = calculate_coverage_depth(alignment)
            log_message("\nCoverage Statistics:", log_file)
            log_message(f"- Coverage percentage: {coverage_stats['coverage_percent']:.2f}%", log_file)
            log_message(f"- Minimum coverage depth: {coverage_stats['min_coverage']}", log_file)
            log_message(f"- Maximum coverage depth: {coverage_stats['max_coverage']}", log_file)
            log_message(f"- Average coverage depth: {coverage_stats['avg_coverage']:.2f}", log_file)
            log_message(f"- Median coverage depth: {coverage_stats['median_coverage']:.2f}\n", log_file)

            # Apply filters
            log_message(separator, log_file)
            log_message("Filtering Progress:", log_file)
            log_message(separator, log_file)
            
            filter_result = apply_filters(
                alignment,
                consensus_threshold=consensus_threshold,
                human_threshold=human_threshold,
                at_threshold=at_threshold,
                at_mode=at_mode,
                outlier_percentile=outlier_percentile,
                reference_seq=None if not reference_file else get_reference_sequence(reference_file),
                enable_human=enable_human,
                enable_at=enable_at,
                enable_outliers=enable_outliers,
                enable_reference=enable_reference
            )
            
            stats['filter_result'] = filter_result
            
            # Write outputs
            log_message(separator, log_file)
            log_message("Writing Output Files:", log_file)
            log_message(separator, log_file)
            
            if filter_result.kept_records:
                sorted_kept_records = sort_sequences_by_position(filter_result.kept_records)
                write_sequences_to_fasta(sorted_kept_records, output_paths['cleaned'])
                log_message(f"Wrote {len(sorted_kept_records)} kept sequences to: {output_paths['cleaned']}", log_file)
                
                consensus_record = SeqRecord(
                    Seq(filter_result.consensus_seq),
                    id=f"{base_name}_fcleaner",
                    description=f"c{consensus_threshold}_h{human_threshold}_a{at_threshold}_p{outlier_percentile}"
                )
                write_sequences_to_fasta([consensus_record], output_paths['consensus'])
                log_message(f"Wrote consensus sequence to: {output_paths['consensus']}", log_file)
            
            # Track removed sequence counts for stats without writing files
            all_removed_count = 0
            for reason, records in filter_result.removed_records.items():
                if records:
                    all_removed_count += len(records)
                    stats['removed_sequences'][reason] = len(records)
                    log_message(f"Filtered out {len(records)} {reason} sequences", log_file)
            
            if all_removed_count:
                log_message(f"Total of {all_removed_count} sequences were filtered out", log_file)
            
            write_ordered_annotated_alignment(
                filter_result.kept_records,
                filter_result.removed_records,
                output_paths['ordered_annotated']
            )
            log_message(f"Wrote ordered and annotated alignment to: {output_paths['ordered_annotated']}", log_file)
            
            write_metrics_report(filter_result, output_paths['metrics'])
            log_message(f"Wrote metrics report to: {output_paths['metrics']}", log_file)
            
            # Update statistics and log final results
            stats['kept_sequences'] = len(filter_result.kept_records)
            stats['processing_time'] = (datetime.now() - start_time).total_seconds()
            
            log_message("\nFinal Results:", log_file)
            log_message(f"Input sequences: {stats['input_sequences']}", log_file)
            log_message(f"Kept sequences: {stats['kept_sequences']}", log_file)
            for reason, count in stats['removed_sequences'].items():
                log_message(f"Removed ({reason}): {count}", log_file)
            
            log_message(f"\nProcessing completed in {stats['processing_time']:.2f} seconds", log_file)
            log_message(separator + "\n", log_file)
            
        except Exception as e:
            # Simplify the error - DO NOT print to stdout here, just log it
            error_msg = f"Error processing {original_filename}: {str(e)}"
            log_message(error_msg, log_file, stdout=False)  # Change stdout=False
            log_message(f"Traceback:\n{traceback.format_exc()}", log_file)
            raise RuntimeError(error_msg)  # Pass the error up without printing
    
    return stats

def parallel_process_directory(args: argparse.Namespace) -> None:
    """
    Process multiple FASTA files in parallel with progress tracking.
    """
    
    # Create output subdirectories
    output_subdirs = create_output_subdirectories(args.output_dir)
    
    # Get list of FASTA files
    fasta_files = [
        f for f in os.listdir(args.input_dir)
        if f.lower().endswith(('.fasta', '.fas', '.fa'))
    ]
    
    total_files = len(fasta_files)
    print(f"\nProcessing {total_files} FASTA files using {args.threads} threads")
    print("=" * 80)
    
    # Track progress
    completed_files = 0
    error_files = []
    skipped_files = []
    processed_files = []
    start_time = datetime.now()
    
    # Open a main log file for recording skipped files
    main_log_path = os.path.join(output_subdirs['logs'], 'fasta_cleaner.log')
    with open(main_log_path, 'w') as main_log:
        main_log.write(f"FASTA Cleaner Processing Log - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        main_log.write("=" * 80 + "\n\n")
        
        # Check for already processed files
        for filename in fasta_files:
            base_name = get_base_filename(filename)
            
            if check_file_already_processed(base_name, output_subdirs):
                skipped_files.append(base_name)
                log_message = f"Skipping already processed file: {filename}"
                main_log.write(f"{log_message}\n")
                main_log.flush()
                print(f"- {log_message}")
        
        if skipped_files:
            main_log.write(f"\nTotal files skipped: {len(skipped_files)}\n")
            main_log.write("=" * 80 + "\n\n")
            main_log.flush()
            print(f"\nSkipped {len(skipped_files)} already processed files.")
    
    def update_progress(filename: str, success: bool = True):
        nonlocal completed_files
        completed_files += 1
        if not success:
            error_files.append(filename)
        else:
            processed_files.append(filename)
        
        # Calculate progress percentage and estimated time
        progress = (completed_files / total_files) * 100
        elapsed_time = (datetime.now() - start_time).total_seconds()
        files_per_second = completed_files / elapsed_time if elapsed_time > 0 else 0
        estimated_total_time = total_files / files_per_second if files_per_second > 0 else 0
        time_remaining = max(0, estimated_total_time - elapsed_time)
        
        # Create progress bar
        bar_width = 50
        filled = int(bar_width * completed_files // total_files)
        bar = '=' * filled + '-' * (bar_width - filled)
        
        # Format time remaining
        hours = int(time_remaining // 3600)
        minutes = int((time_remaining % 3600) // 60)
        seconds = int(time_remaining % 60)
        time_str = f"{hours:02d}:{minutes:02d}:{seconds:02d}"
        
        # Clear line and update progress
        sys.stdout.write('\r' + ' ' * 100 + '\r')  # Clear the line
        progress_msg = f"Progress: [{bar}] {progress:.1f}% ({completed_files}/{total_files}) ETA: {time_str}"
        sys.stdout.write(progress_msg)
        sys.stdout.flush()
        
        # Print newline if there are errors
        if not success:
            sys.stdout.write('\n')
            sys.stdout.flush()
    
    # Create partial function with common arguments
    process_func = partial(
        process_fasta_file,
        output_dir=args.output_dir,
        consensus_threshold=args.consensus_threshold,
        human_threshold=args.human_threshold,
        at_threshold=args.at_difference,
        at_mode=args.at_mode,
        outlier_percentile=args.percentile_threshold,
        enable_human=not args.disable_human,
        enable_at=not args.disable_at,
        enable_outliers=not args.disable_outliers,
        enable_reference=True if args.reference_dir else False
    )
    
    # Process files in parallel, skipping already processed ones
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {}
        for filename in fasta_files:
            input_path = os.path.join(args.input_dir, filename)
            base_name = get_base_filename(filename)
            
            # Skip if already processed
            if base_name in skipped_files:
                continue
            
            # Get reference file if applicable
            reference_file = None
            if args.reference_dir:
                reference_path = os.path.join(args.reference_dir, f"{base_name}_reference.fasta")
                if os.path.exists(reference_path):
                    reference_file = reference_path
            
            future = executor.submit(process_func, input_path, reference_file=reference_file)
            futures[future] = filename
        
        # Process results as they complete
        try:
            for future in concurrent.futures.as_completed(futures):
                filename = futures[future]
                try:
                    stats = future.result()
                    update_progress(filename, success=True)
                except Exception as e:
                    update_progress(filename, success=False)
                    # Print only the simple error message, avoiding nesting
                    error_msg = str(e)
                    # If we have a nested error message, simplify it
                    if "Error processing" in error_msg and "Error processing" in error_msg[20:]:
                        # Extract just the innermost error
                        innermost_error = error_msg.split(": ")[-1]
                        print(f"ERROR: Error processing {filename}: {innermost_error}")
                    else:
                        print(f"ERROR: {error_msg}")
        except KeyboardInterrupt:
            print("\nProcessing interrupted by user. Waiting for running tasks to complete...")
            executor.shutdown(wait=True)
            print("Shutdown complete.")
            sys.exit(1)
    
    # Final newline and summary
    print("\n" + "=" * 80)
    elapsed_time = (datetime.now() - start_time).total_seconds()
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = int(elapsed_time % 60)
    
    # Update the main log with final statistics
    with open(main_log_path, 'a') as main_log:
        main_log.write(f"\nProcessing Summary - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        main_log.write("=" * 80 + "\n\n")
        main_log.write(f"Processing complete in {hours:02d}:{minutes:02d}:{seconds:02d}\n")
        main_log.write(f"Total files processed: {len(processed_files)}/{total_files}\n")
        main_log.write(f"Total files skipped: {len(skipped_files)}/{total_files}\n")
        if error_files:
            main_log.write(f"Errors occurred in {len(error_files)} files:\n")
            for filename in error_files:
                main_log.write(f"- {filename}\n")
    
    print(f"Processing complete in {hours:02d}:{minutes:02d}:{seconds:02d}")
    print(f"Total files processed: {len(processed_files)}/{total_files}")
    print(f"Total files skipped: {len(skipped_files)}/{total_files}")
    if error_files:
        print(f"Errors occurred in {len(error_files)} files:")
        for filename in error_files:
            print(f"- {filename}")
    print("=" * 80)
    
    print("\nConcatenating consensus sequences...")
    concatenated_consensus_file = os.path.join(args.output_dir, 'all_consensus_sequences.fasta')
    concatenate_consensus_sequences(
        output_subdirs['consensus_seqs'],
        concatenated_consensus_file,
        args.preprocessing
    )
    print(f"Concatenated consensus sequences written to: {concatenated_consensus_file}")
    
    # Write the summary file
    total_stats = defaultdict(int)
    total_stats['processed_files'] = len(processed_files) + len(skipped_files)  # Include skipped files
    total_stats['error_files'] = len(error_files)
    write_summary(total_stats, output_subdirs['logs'])
    
    # Generate combined statistics CSV file
    print("\nGenerating combined statistics CSV file...")
    generate_combined_statistics(args.output_dir)



# =============================================================================
# Command Line Interface
# =============================================================================
def parse_arguments() -> argparse.Namespace:
    """Set up command line argument parsing."""
    parser = argparse.ArgumentParser(
        description='Advanced FASTA sequence cleaner with multiple filtering methods: '
                   'human COX1 similarity, AT content, and statistical outlier detection.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Get default number of threads
    default_threads = get_default_threads()
    
    # Required arguments
    parser.add_argument(
        '-i', '--input_dir',
        required=True,
        help='Directory containing input FASTA alignment files'
    )

    # Add threads argument
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=default_threads,
        help='Number of threads to use for parallel processing. '
             'Defaults to SLURM_CPUS_PER_TASK if available, '
             'otherwise uses all available CPU cores.'
    )
    
    parser.add_argument(
        '-o', '--output_dir',
        required=True,
        help='Directory where output files will be saved'
    )
    
    # Filter enablement flags
    parser.add_argument(
        '--disable_human',
        action='store_true',
        help='Disable human COX1 similarity filtering'
    )
    
    parser.add_argument(
        '--disable_at',
        action='store_true',
        help='Disable AT content difference filtering'
    )
    
    parser.add_argument(
        '--disable_outliers',
        action='store_true',
        help='Disable statistical outlier detection'
    )
    
    # Human COX1 filtering parameters
    parser.add_argument(
        '-u', '--human_threshold',
        type=float,
        default=0.95,
        help='Human COX1 similarity threshold (0.0-1.0). Sequences with similarity >= threshold are removed'
    )
    
    # AT content filtering parameters
    parser.add_argument(
        '-d', '--at_difference',
        type=float,
        default=0.1,
        help='Maximum allowed AT content difference from consensus (0.0-1.0)'
    )
    
    parser.add_argument(
        '-m','--at_mode',
        type=str,
        choices=['absolute', 'higher', 'lower'],
        default='absolute',
        help='AT content filtering mode: "absolute" removes sequences if AT content differs '
             'from consensus by more than threshold in either direction, "higher" removes only '
             'sequences with AT content above consensus + threshold, "lower" removes only '
             'sequences with AT content below consensus - threshold'
    )
    
    # Statistical outlier detection parameters
    parser.add_argument(
        '-p', '--percentile_threshold',
        type=float,
        default=90.0,
        help='Percentile threshold for statistical outlier detection (0.0-100.0)'
    )
    
    # Consensus generation parameters
    parser.add_argument(
        '-c', '--consensus_threshold',
        type=float,
        default=0.5,
        help='Threshold for consensus sequence generation (0.0-1.0)'
    )
    
    # Optional reference sequence directory
    parser.add_argument(
        '-r', '--reference_dir',
        help='Optional directory containing reference FASTA files (named same as input files but with "_reference.fasta" suffix)'
    )
    
    parser.add_argument(
        '--preprocessing',
        type=str,
        choices=['concat', 'merge'],
        default='concat',
        help='Preprocessing mode that was used in the workflow. If "merge", "_merge" will be appended to fasta headers'
    )
    
    args = parser.parse_args()
    
    # Validate directory paths
    if not os.path.isdir(args.input_dir):
        parser.error(f"Input directory does not exist: {args.input_dir}")
    
    if args.reference_dir and not os.path.isdir(args.reference_dir):
        parser.error(f"Reference directory does not exist: {args.reference_dir}")
    
    # Validate thresholds
    if not 0 <= args.human_threshold <= 1:
        parser.error("Human threshold must be between 0.0 and 1.0")
    
    if not 0 <= args.at_difference <= 1:
        parser.error("AT difference threshold must be between 0.0 and 1.0")
    
    if not 0 <= args.percentile_threshold <= 100:
        parser.error("Percentile threshold must be between 0.0 and 100.0")
    
    if not 0 < args.consensus_threshold <= 1:
        parser.error("Consensus threshold must be between 0.0 and 1.0")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    return args



# =============================================================================
# Main Entry Point
# =============================================================================
if __name__ == "__main__":
    try:
        # Set the multiprocessing start method
        mp.set_start_method('spawn')  # Required for Bio* objects
        
        # Parse arguments
        args = parse_arguments()
                
        # Run processing
        parallel_process_directory(args)
        
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)
