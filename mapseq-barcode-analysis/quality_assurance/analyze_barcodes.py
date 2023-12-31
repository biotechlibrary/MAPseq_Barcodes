"""
This script performs the following tasks:
1. Reads the extracted barcodes from the output files generated by the fastqDataParsing.py script.
2. Calculates basic statistics, such as the number of unique barcodes and the frequency distribution of the barcodes.
3. Checks the length and format of the barcodes to ensure they are 30 nucleotides long and only contain valid nucleotide characters (A, T, G, and C).

You can run this script after the preprocessing and quality assurance steps to analyze the results and verify the quality of the extracted barcodes.

Usage:
1. Update the `input_directory` variable in the `main()` function to the path of the directory containing the output files generated by the fastqDataParsing.py script.
2. Run this script using a Python interpreter: python3 analyze_barcodes.py
"""

import os
from collections import Counter

def read_barcodes(file_path):
    with open(file_path, 'r') as f:
        barcodes = [line.strip() for line in f.readlines()]
    return barcodes

def basic_statistics(barcodes, top_n=10):
    unique_barcodes = set(barcodes)
    num_unique_barcodes = len(unique_barcodes)
    print(f"Number of unique barcodes: {num_unique_barcodes}")

    barcode_counts = Counter(barcodes)
    print(f"Top {top_n} most frequent barcodes:")
    for barcode, count in barcode_counts.most_common(top_n):
        print(f"{barcode}: {count}")

def check_length_and_format(barcodes):
    incorrect_barcodes = [barcode for barcode in barcodes if len(barcode) != 30 or not set(barcode).issubset(set("ATGC"))]
    if not incorrect_barcodes:
        print("All barcodes have the correct length and format.")
    else:
        print(f"{len(incorrect_barcodes)} barcodes have incorrect length or format.")

def main():
    input_directory = "path/to/output_directory"  # This should be the same output directory used in fastq_data_parsing.py

    for file_name in os.listdir(input_directory):
        if file_name.endswith("_barcodes.txt"):
            file_path = os.path.join(input_directory, file_name)
            print(f"Analyzing {os.path.splitext(file_name)[0]}:")
            
            barcodes = read_barcodes(file_path)
            basic_statistics(barcodes)
            check_length_and_format(barcodes)
            print("\n")

if __name__ == "__main__":
    main()
