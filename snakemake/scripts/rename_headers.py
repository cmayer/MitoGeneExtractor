#!/usr/bin/env python3
"""
rename_fasta_headers.py - A script to rename FASTA headers in consensus files

This script renames the headers in FASTA consensus files based on standardized naming conventions.
It's designed to be called by a Snakemake rule.
"""

import os
import sys
import signal
import argparse
import threading
import time
import atexit
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor
from threading import Lock, Event

# Lock for thread-safe file writing
log_lock = Lock()

# Global variables for process management
all_tasks_complete = Event()
early_termination = False
executor = None
timeout_timer = None

def log_message(msg, log_file):
    """Write message to both stdout and log file in a thread-safe manner"""
    print(msg, flush=True)
    with log_lock:
        with open(log_file, 'a') as f:
            f.write(f"{msg}\n")

def cleanup_resources():
    """Cleanup function to ensure all resources are properly released"""
    global executor, timeout_timer
    
    if timeout_timer is not None and timeout_timer.is_alive():
        log_message("Cancelling timeout timer", "cleanup.log")
        timeout_timer.cancel()
    
    if executor is not None and not executor._shutdown:
        log_message("Shutting down thread pool executor", "cleanup.log")
        executor.shutdown(wait=False)

def signal_handler(signum, frame):
    """Handle termination signals gracefully"""
    log_message(f"Received signal {signum}, performing cleanup", "cleanup.log")
    cleanup_resources()
    sys.exit(0)

def process_file(file_info):
    """Process a single file, to be called by the thread pool"""
    input_file, output_file, preprocessing_mode, log_file = file_info
    
    try:
        base_name = os.path.basename(input_file)
        parts = base_name.split('_')
        
        # Extract relevant parts from filename
        sample = parts[0]
        r = parts[2]
        s = parts[4]
        con_suffix = base_name.split('con_')[-1].replace('.fas', '')
        
        # Add suffix based on preprocessing mode
        suffix = "_merge" if preprocessing_mode == "merge" else ""
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    new_header = f">{sample}_r_{r}_s_{s}_{con_suffix}{suffix}\n"
                    outfile.write(new_header)
                else:
                    outfile.write(line)
        
        return output_file, None  # Return output file path and no error
    except Exception as e:
        return output_file, str(e)  # Return output file path and error message

def check_completion(output_files, complete_file, log_file):
    """
    Check if all expected outputs have been created and terminate if they have
    
    Args:
        output_files (list): List of expected output files
        complete_file (str): Path to completion status file
        log_file (str): Path to log file
    """
    global early_termination
    
    log_message("Timeout reached, checking if all tasks are complete...", log_file)
    
    # Check if all output files exist
    missing_files = []
    for output_file in output_files:
        if not os.path.exists(output_file):
            missing_files.append(output_file)
    
    # Check if completion file exists or is about to be created
    completion_file_exists = os.path.exists(complete_file) or all_tasks_complete.is_set()
    
    if not missing_files and (completion_file_exists or all_tasks_complete.is_set()):
        log_message("All expected outputs found. Script appears to have completed successfully.", log_file)
        log_message("Terminating process gracefully due to timeout check.", log_file)
        
        # Ensure completion file exists
        if not os.path.exists(complete_file):
            create_completion_file(complete_file, "Timeout-triggered completion", 0)
        
        early_termination = True
        cleanup_resources()
        os._exit(0)  # Force immediate exit
    else:
        if missing_files:
            log_message(f"Timeout check: Still missing {len(missing_files)} output files. Continuing to wait...", log_file)
            if len(missing_files) <= 10:
                for missing in missing_files:
                    log_message(f"  Missing: {missing}", log_file)
        if not completion_file_exists:
            log_message("Completion file not yet created", log_file)

def create_completion_file(complete_file, status_message, error_count, start_time=None, end_time=None):
    """Create a standardized completion file"""
    if start_time is None:
        start_time = datetime.now()
    if end_time is None:
        end_time = datetime.now()
        
    duration = end_time - start_time
    
    os.makedirs(os.path.dirname(complete_file), exist_ok=True)
    with open(complete_file, 'w') as f:
        f.write(f"Rename operation started at: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Rename operation completed at: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total duration: {duration}\n")
        f.write(f"Status: {status_message}\n")
        f.write(f"Errors: {error_count}")

def rename_fasta_headers(input_files, output_files, complete_file, log_file, preprocessing_mode="concat", num_threads=1, timeout_minutes=10):
    """
    Rename headers in FASTA consensus files using parallel processing.
    
    Args:
        input_files (list): List of input consensus FASTA files
        output_files (list): List of output renamed FASTA files
        complete_file (str): Path to write completion status
        log_file (str): Path to log file
        preprocessing_mode (str): Mode of preprocessing ('merge' or 'concat')
        num_threads (int): Number of threads to use for parallel processing
        timeout_minutes (int): Number of minutes after which to check if tasks are complete
    """
    global executor, timeout_timer, all_tasks_complete
    
    start_time = datetime.now()
    log_message(f"Starting rename_fasta_headers at: {start_time.strftime('%Y-%m-%d %H:%M:%S')}", log_file)
    log_message(f"Processing {len(input_files)} files using {num_threads} threads", log_file)
    log_message(f"Timeout check set for {timeout_minutes} minutes", log_file)
    
    # Set up timeout timer to check completion status
    timeout_seconds = timeout_minutes * 60
    timeout_timer = threading.Timer(
        timeout_seconds, 
        check_completion, 
        args=[output_files, complete_file, log_file]
    )
    timeout_timer.daemon = True  # Allow the timer to be terminated when the main thread exits
    timeout_timer.start()
    
    # Prepare arguments for parallel processing
    file_args = [(input_file, output_file, preprocessing_mode, log_file) 
                 for input_file, output_file in zip(input_files, output_files)]
    
    processed_count = 0
    errors = []
    
    # Process files in parallel with managed executor
    executor = ThreadPoolExecutor(max_workers=num_threads)
    try:
        for output_file, error in executor.map(process_file, file_args):
            processed_count += 1
            
            if error:
                log_message(f"Error processing {output_file}: {error}", log_file)
                errors.append((output_file, error))
            
            # Log progress periodically
            if processed_count % 50 == 0 or processed_count == len(input_files):
                log_message(f"Processed {processed_count}/{len(input_files)} files", log_file)
    finally:
        # Properly shutdown the executor
        executor.shutdown(wait=True)
        executor = None
    
    # Cancel the timer if we've completed normally
    if timeout_timer and timeout_timer.is_alive():
        timeout_timer.cancel()
        timeout_timer = None
    
    # Report any missing files
    missing_files = []
    for output_file in output_files:
        if not os.path.exists(output_file):
            missing_files.append(output_file)
            log_message(f"Warning: Output file missing: {output_file}", log_file)
    
    if missing_files:
        log_message(f"Warning: {len(missing_files)} output files are missing!", log_file)
    else:
        log_message("All output files successfully created", log_file)

    if errors:
        log_message(f"Completed with {len(errors)} errors", log_file)
    else:
        log_message("Rename operation completed successfully with no errors", log_file)

    end_time = datetime.now()
    duration = end_time - start_time
    log_message(f"Operation duration: {duration}", log_file)
    
    # Create checkpoint file
    create_completion_file(complete_file, "Completed normally", len(errors), start_time, end_time)
    log_message("Checkpoint file created", log_file)
    
    # Set completion flag
    all_tasks_complete.set()

def main():
    """Main function to parse arguments and run the script"""
    # Register signal handlers for graceful termination
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)
    
    # Register cleanup function to run on exit
    atexit.register(cleanup_resources)
    
    parser = argparse.ArgumentParser(description='Rename headers in FASTA consensus files')
    
    # Arguments for file list
    parser.add_argument('--input-file-list', help='File containing list of input FASTA files')
    parser.add_argument('--output-file-list', help='File containing list of output FASTA files')
    
    # Arguments for direct file specification
    parser.add_argument('--input-files', nargs='*', help='Space-separated list of input FASTA files')
    parser.add_argument('--output-files', nargs='*', help='Space-separated list of output FASTA files')
    
    # Required arguments
    parser.add_argument('--complete-file', required=True, help='Path to write completion status')
    parser.add_argument('--log-file', required=True, help='Path to log file')
    parser.add_argument('--preprocessing-mode', default='concat', choices=['concat', 'merge'], 
                        help='Preprocessing mode (concat or merge)')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    
    # Timeout argument
    parser.add_argument('--timeout-minutes', type=int, default=10, 
                        help='Number of minutes after which to check for completion')
    
    args = parser.parse_args()
    
    input_files = []
    output_files = []
    
    # Check if we have input file list
    if args.input_file_list and args.output_file_list:
        # Read input files list
        with open(args.input_file_list, 'r') as f:
            input_files = [line.strip() for line in f if line.strip()]
        
        # Read output files list
        with open(args.output_file_list, 'r') as f:
            output_files = [line.strip() for line in f if line.strip()]
    
    # Check if we have direct input files
    elif args.input_files and args.output_files:
        input_files = args.input_files
        output_files = args.output_files
    
    # If neither option is provided
    else:
        error_msg = "Error: You must provide either --input-file-list and --output-file-list OR --input-files and --output-files"
        print(error_msg)
        with open(args.log_file, 'a') as f:
            f.write(f"{error_msg}\n")
        sys.exit(1)
    
    # Validate that input and output file lists have the same length
    if len(input_files) != len(output_files):
        error_msg = f"Error: Number of input files ({len(input_files)}) must match number of output files ({len(output_files)})"
        print(error_msg)
        with open(args.log_file, 'a') as f:
            f.write(f"{error_msg}\n")
        sys.exit(1)
    
    try:
        rename_fasta_headers(
            input_files,
            output_files,
            args.complete_file,
            args.log_file,
            args.preprocessing_mode,
            args.threads,
            args.timeout_minutes
        )
        
        # Final cleanup
        cleanup_resources()
        
        # Log successful completion
        log_message("Script completed successfully and is terminating", args.log_file)
        
    except Exception as e:
        if early_termination:
            log_message("Script terminated early after verifying successful completion", args.log_file)
            sys.exit(0)
        else:
            log_message(f"Error: {str(e)}", args.log_file)
            # Log the traceback
            import traceback
            log_message(traceback.format_exc(), args.log_file)
            
            # Create failure completion file
            create_completion_file(args.complete_file, f"Failed with error: {str(e)}", 1)
            
            sys.exit(1)
    finally:
        # Make absolutely sure we exit properly by using os._exit if necessary
        time.sleep(1)  # Brief pause to allow logging to complete
        os._exit(0)  # Force termination of all threads

if __name__ == '__main__':
    main()
