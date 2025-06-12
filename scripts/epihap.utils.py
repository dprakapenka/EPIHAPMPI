# epihap.utils.py
# Utility functions for epihap scripts

import os
import sys
import re
import argparse
import pandas as pd
from collections import defaultdict

# ==============================================================================
# String/File Manipulation Utilities
# ==============================================================================

def get_digits(s):
    """
    Extracts integer digits from a string. Returns 0 if no digits found.

    Args:
        s (str): The input string.

    Returns:
        int: The integer formed by the digits in the string, or 0.
    """
    digits = ''.join(filter(str.isdigit, str(s)))
    return int(digits) if digits else 0

def sort_files(file_list, enable_sorting=True):
    """
    Sorts a list of file paths based on digits found within the filenames.

    Args:
        file_list (list): A list of file path strings.
        enable_sorting (bool, optional): If True (default), sorts the list.
                                         If False, returns the list as is.

    Returns:
        list: The sorted (or original) list of file paths.
    """
    if enable_sorting and file_list:
        # Sort primarily by digits, secondarily by the full string for tie-breaking
        return sorted(file_list, key=lambda f: (get_digits(os.path.basename(f)), f))
    else:
        return file_list

def parse_block_position_line(line):
    """
    Parses a line string in the format 'chr:begin_pos:end_pos'.

    Args:
        line (str): The input string line.

    Returns:
        tuple: A tuple containing (chromosome_str, begin_pos_int, end_pos_int),
               or None if parsing fails.
    """
    try:
        # Use regex for flexibility with separators (:, -, etc.) and whitespace
        match = re.match(r'^\s*(\w+)\D+(\d+)\D+(\d+)\s*$', line)
        if match:
            chrom = match.group(1)
            begin_pos = int(match.group(2))
            end_pos = int(match.group(3))
            return chrom, begin_pos, end_pos
        else:
            print(f"Warning: Could not parse block position line: {line.strip()}", file=sys.stderr)
            return None
    except Exception as e:
        print(f"Error parsing block position line '{line.strip()}': {e}", file=sys.stderr)
        return None

# ==============================================================================
# File I/O Utilities
# ==============================================================================

def prepare_output_directory(dir_path, verbose=False):
    """
    Creates an output directory if it doesn't exist.

    Args:
        dir_path (str): The path to the directory to create.
        verbose (bool, optional): If True, prints a message upon creation. Defaults to False.

    Returns:
        bool: True if the directory exists or was created successfully, False otherwise.
    """
    try:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            if verbose:
                print(f"Created output directory: {dir_path}")
        return True
    except OSError as e:
        print(f"Error creating directory {dir_path}: {e}", file=sys.stderr)
        return False

def check_file_validity(file_path):
    """
    Check if a file exists and is not empty.

    Args:
        file_path (str): The path to the file.

    Returns:
        bool: True if the file exists and is not empty, False otherwise.
    """
    if not os.path.exists(file_path):
        print(f"Warning: File not found: {file_path}", file=sys.stderr)
        return False
    if os.path.getsize(file_path) == 0:
        print(f"Warning: File is empty: {file_path}", file=sys.stderr)
        return False
    return True

def read_map_file(map_filepath, noheader=False):
    """
    Reads a map file (SNPID Chr Pos ...) and organizes it by chromosome.

    Args:
        map_filepath (str): Path to the map file.
        noheader (bool, optional): True if the map file has no header line. Defaults to False.

    Returns:
        dict: A dictionary where keys are chromosome identifiers (str) and
              values are lists of [snp_index_within_chr, snp_id, position_int].
              Returns an empty dict if the file cannot be opened or is empty.
    """
    snp_dict = defaultdict(list)
    if not check_file_validity(map_filepath):
        return snp_dict

    try:
        with open(map_filepath, "r") as f:
            if not noheader:
                f.readline()  # Skip header

            snp_indices = defaultdict(int) # Track index within each chromosome
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    snp_id = parts[0]
                    chrom = parts[1]
                    try:
                        pos = int(parts[2])
                        snp_idx = snp_indices[chrom]
                        snp_dict[chrom].append([snp_idx, snp_id, pos])
                        snp_indices[chrom] += 1
                    except ValueError:
                        print(f"Warning: Skipping line due to non-integer position in {map_filepath}: {line.strip()}", file=sys.stderr)
                else:
                     print(f"Warning: Skipping malformed line in {map_filepath}: {line.strip()}", file=sys.stderr)

    except IOError:
        print(f"Error: Can't open map file {map_filepath}", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred while reading {map_filepath}: {e}", file=sys.stderr)

    if not snp_dict:
         print(f"Warning: No valid SNP data read from map file: {map_filepath}", file=sys.stderr)

    return snp_dict

def read_hap_info_file(hap_info_filepath):
    """
    Reads a haplotype info file (e.g., blk_id start_snp_idx end_snp_idx).

    Args:
        hap_info_filepath (str): Path to the hap_info file.

    Returns:
        list: A list of lists, where each inner list contains
              [start_snp_index_int, end_snp_index_int].
              Returns an empty list if the file cannot be opened or is empty.
    """
    blocks = []
    if not check_file_validity(hap_info_filepath):
        return blocks

    try:
        with open(hap_info_filepath, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3: # Expecting at least blk_id, start, end
                    try:
                        start_idx = int(parts[-2]) # Second to last element
                        end_idx = int(parts[-1])   # Last element
                        if start_idx <= end_idx:
                             blocks.append([start_idx, end_idx])
                        else:
                             print(f"Warning: Skipping block with start > end in {hap_info_filepath}: {line.strip()}", file=sys.stderr)
                    except ValueError:
                        print(f"Warning: Skipping line due to non-integer indices in {hap_info_filepath}: {line.strip()}", file=sys.stderr)
                else:
                     print(f"Warning: Skipping malformed line in {hap_info_filepath}: {line.strip()}", file=sys.stderr)
    except IOError:
        print(f"Error: Can't open hap_info file {hap_info_filepath}", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred while reading {hap_info_filepath}: {e}", file=sys.stderr)

    if not blocks:
         print(f"Warning: No valid block data read from hap_info file: {hap_info_filepath}", file=sys.stderr)

    return blocks

def read_hap_file(hap_filepath, header=False):
    """
    Reads a haplotype file (ID allele_pairs...).

    Args:
        hap_filepath (str): Path to the haplotype file.
        header (bool, optional): True if the file has a header line. Defaults to False.

    Returns:
        list: A list of tuples, where each tuple is (individual_id_str, haplotype_str).
              Returns an empty list if the file cannot be opened or is empty.
    """
    haps_data = []
    if not check_file_validity(hap_filepath):
        return haps_data

    try:
        with open(hap_filepath, 'r') as f:
            if header:
                f.readline() # Skip header
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    ind_id = parts[0]
                    # Join remaining parts in case ID or hap string contains spaces
                    hap_string = "".join(parts[1:])
                    if len(hap_string) % 2 != 0:
                         print(f"Warning: Haplotype string for ID {ind_id} in {hap_filepath} has odd length. Skipping.", file=sys.stderr)
                         continue
                    haps_data.append((ind_id, hap_string))
                else:
                    print(f"Warning: Skipping malformed line in {hap_filepath}: {line.strip()}", file=sys.stderr)

    except IOError:
        print(f"Error: Can't open haplotype file {hap_filepath}", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred while reading {hap_filepath}: {e}", file=sys.stderr)

    if not haps_data:
         print(f"Warning: No valid haplotype data read from file: {hap_filepath}", file=sys.stderr)

    return haps_data

def read_file_with_auto_delimiter(file_path):
    """
    Read a file and automatically detect its delimiter (comma, tab, or space).
    Uses pandas for robust reading.

    Args:
        file_path (str): The path to the file.

    Returns:
        pandas.DataFrame: A DataFrame containing the file data, or None if reading fails.
    """
    if not check_file_validity(file_path):
        return None

    try:
        # Read the first line to guess delimiter
        with open(file_path, 'r') as f:
            first_line = f.readline()
            if ',' in first_line:
                delimiter = ','
            elif '\t' in first_line:
                 delimiter = '\t'
            else:
                 delimiter = r'\s+' # Regex for one or more whitespace characters

        # Read using pandas with the detected delimiter
        df = pd.read_csv(file_path, sep=delimiter)
        return df
    except pd.errors.EmptyDataError:
        print(f"Warning: File is empty or contains no data: {file_path}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error reading file {file_path} with pandas: {e}", file=sys.stderr)
        return None

# ==============================================================================
# Logging Utility
# ==============================================================================

def log_message(message, verbose=False, logfile=None):
    """
    Prints a message to stdout if verbose is True and appends to a logfile if provided.

    Args:
        message (str): The message to log.
        verbose (bool, optional): If True, print to stdout. Defaults to False.
        logfile (str, optional): Path to the log file. If provided, appends the message.
                                 Defaults to None.
    """
    if verbose:
        print(message)
    if logfile:
        try:
            with open(logfile, 'a') as logf:
                logf.write(str(message) + "\n")
        except IOError as e:
            print(f"Error: Could not write to log file {logfile}: {e}", file=sys.stderr)

# ==============================================================================
# Example usage
# ==============================================================================
# if __name__ == '__main__':
#     # Example usage
#     print(f"Digits in 'chr123file': {get_digits('chr123file')}")
#     print(f"Digits in 'no_digits_here': {get_digits('no_digits_here')}")

#     files = ['file1.txt', 'file10.txt', 'file2.txt']
#     print(f"Original files: {files}")
#     print(f"Sorted files: {sort_files(files)}")
#     print(f"Unsorted files: {sort_files(files, enable_sorting=False)}")

#     # Create dummy files for testing read functions 
#     # prepare_output_directory("temp_test_dir", verbose=True)
#     # with open("temp_test_dir/test_map.txt", "w") as f:
#     #     f.write("SNPID\tChr\tPos\n")
#     #     f.write("snp1\t1\t100\n")
#     #     f.write("snp2\t1\t200\n")
#     #     f.write("snp3\t2\t50\n")
#     # map_data = read_map_file("temp_test_dir/test_map.txt")
#     # print(f"Map data: {map_data}")

# ==============================================================================
# Argument Parsing Utilities
# ==============================================================================

def create_common_block_parser():
    """
    Creates an ArgumentParser with common arguments for block-based scripts.

    Returns:
        argparse.ArgumentParser: An ArgumentParser instance with predefined
                                 --map, --output, and --verbose arguments.
    """
    parser = argparse.ArgumentParser(add_help=False)  # add_help=False to allow overriding in child parsers
    
    # Required arguments
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument("-m", "--map",
                                required=True,
                                help="path to map file [required]\n"
                                     "format: 'SNPID Chr Pos' with header")

    # Optional arguments
    optional_group = parser.add_argument_group('optional arguments')
    optional_group.add_argument("-o", "--output",
                                default="hap_info",
                                help="output file name (default: hap_info)")
    optional_group.add_argument("-V", "--verbose",
                                action="store_true",
                                default=False,
                                help="verbose output")
    return parser
