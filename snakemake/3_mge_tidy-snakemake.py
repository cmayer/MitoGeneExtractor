import os
import sys
import shutil



def process_files_in_directory(directory, alignment_dir, consensus_dir, out_dir, err_dir, logs_dir, consensus_files):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)

        if filename.lower().endswith('.fas'):
            if 'align' in filename:
                shutil.move(file_path, os.path.join(alignment_dir, filename))
            elif 'con' in filename:
                rename_sequence_headers(file_path)
                shutil.move(file_path, os.path.join(consensus_dir, filename))
                consensus_files.append(os.path.join(consensus_dir, filename))
        elif filename.lower().endswith('.out'):
            shutil.move(file_path, os.path.join(out_dir, filename))
        elif filename.lower().endswith('.err'):
            shutil.move(file_path, os.path.join(err_dir, filename))
        elif filename.lower().endswith(('.log', '.txt')):
            shutil.move(file_path, os.path.join(logs_dir, filename))



def rename_sequence_headers(fasta_file):
    temp_file = fasta_file + ".tmp"
    filename = os.path.basename(fasta_file)
    prefix = filename[:11]

    with open(fasta_file, 'r') as infile, open(temp_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                header = line.strip().replace('Consensus__', '')
                new_header = f">{prefix}_{header[1:]}"
                outfile.write(new_header + '\n')
            else:
                outfile.write(line)

    os.replace(temp_file, fasta_file)



def concatenate_fasta_files(files_to_concatenate, output_file):
    with open(output_file, 'w') as outfile:
        for file in files_to_concatenate:
            with open(file, 'r') as infile:
                outfile.write(infile.read())



def main():
    if len(sys.argv) != 4:
        print(
        """
        Usage: python 4_mge_tidy-snakemake.py /path/to/cox1 cox1_concat_cons.fasta user_named_directory
        
        user_named_directory = directory to be created within cox1/ 
        """
        )

        sys.exit(1)

    directory_path = sys.argv[1]
    concatenated_filename = sys.argv[2]
    user_named_directory = sys.argv[3]

    results_dir_path = os.path.join(directory_path, user_named_directory)

    if not os.path.exists(results_dir_path):
        os.makedirs(results_dir_path)

    sub_dirs = ["alignment", "consensus", "out", "logs", "err"]
    sub_dir_paths = {sub_dir: os.path.join(results_dir_path, sub_dir) for sub_dir in sub_dirs}

    for sub_dir in sub_dir_paths.values():
        if not os.path.exists(sub_dir):
            os.makedirs(sub_dir)
            print(f"Created directory: {sub_dir}")

    consensus_files = []

    process_files_in_directory(directory_path, sub_dir_paths['alignment'], sub_dir_paths['consensus'], sub_dir_paths['out'], sub_dir_paths['err'], sub_dir_paths['logs'], consensus_files)

    if consensus_files:
        concatenate_fasta_files(consensus_files, os.path.join(sub_dir_paths['consensus'], concatenated_filename))
        print(f"Concatenated consensus files saved to {os.path.join(sub_dir_paths['consensus'], concatenated_filename)}")

    print("Organised files into directories.")


if __name__ == "__main__":
    main()
