import os
import shutil
from concurrent.futures import ThreadPoolExecutor

def delete_zero_size_file(file_path):
    try:
        if os.path.isfile(file_path) and os.path.getsize(file_path) == 0:
            os.remove(file_path)
            return f"Deleted: {file_path}"
    except Exception as e:
        print(f"Error deleting file {file_path}: {e}")
    return None

def get_atom_count(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            if lines:
                try:
                    atom_count = int(lines[0].strip())
                    return atom_count
                except ValueError:
                    return 0
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
    return 0

def process_file(file_path, destination_dir, atom_threshold=75):
    try:
        atom_count = get_atom_count(file_path)
        if atom_count > atom_threshold:
            target_path = os.path.join(destination_dir, os.path.basename(file_path))
            if file_path != target_path and (not os.path.exists(target_path) or os.path.getsize(target_path) == 0):
                shutil.move(file_path, target_path)
                return f"Moved: {file_path} to {destination_directory}"
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
    return None

def main(source_directory, destination_directory, atom_threshold=75):
    if not os.path.exists(destination_directory):
        os.makedirs(destination_directory)

    with ThreadPoolExecutor() as executor:
        # Delete zero size files and process files in parallel
        all_files = [os.path.join(source_directory, f) for f in os.listdir(source_directory)]
        delete_results = executor.map(delete_zero_size_file, all_files)
        process_results = executor.map(lambda f: process_file(f, destination_directory, atom_threshold), all_files)

    # Output the results
    for result in delete_results:
        if result:
            print(result)
    for result in process_results:
        if result:
            print(result)

source_directory = r'C:\Users\szou5\Desktop\work\XYZ'
destination_directory = r'C:\Users\szou5\Desktop\work\XYZ\Big'

main(source_directory, destination_directory)