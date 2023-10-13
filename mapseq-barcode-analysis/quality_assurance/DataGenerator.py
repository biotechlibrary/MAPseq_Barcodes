'''This Script generates synthetic MAPseq data for testing.''' 

import random
import os

def generate_random_sequence(length):
    nucleotides = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(nucleotides) for _ in range(length))

def generate_fastq_file(filename, num_sequences, barcode_length, anchor_sequence, read_length):
    with open(filename, 'w') as file:
        for i in range(num_sequences):
          
            # Generate a random barcode sequence
            barcode = generate_random_sequence(barcode_length)

            # Generate a random read sequence with the barcode and anchor sequence
            read_sequence = barcode + anchor_sequence
            remaining_bases = read_length - len(read_sequence)
            read_sequence += generate_random_sequence(remaining_bases)

            # Write the sequence in fastq format
            file.write(f"@seq{i}\n")
            file.write(f"{read_sequence}\n")
            file.write(f"+\n")
            file.write(f"{'I' * read_length}\n")  # Use constant quality score 'I'

if __name__ == "__main__":
    output_file = "synthetic_data.fastq"
    num_sequences = 100000
    barcode_length = 30
    anchor_sequence = "GTACTGCGGCCGCTACCTA"
    read_length = 50

    generate_fastq_file(output_file, num_sequences, barcode_length, anchor_sequence, read_length)
