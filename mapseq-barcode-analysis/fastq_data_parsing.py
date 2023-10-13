"""
THIS SCRIPT PARSES FASTQ FILES AND EXTRACTS THE BARCODE SEQUENCES USING THE ANCHOR SEQUENCES.
"""

import os
import gzip
from Bio import SeqIO
from multiprocessing import Pool

class BarcodeExtractor:
    """
    This class provides methods for extracting barcodes from FASTQ files based on a given anchor sequence.
    """

    def __init__(self, input_directory, output_directory, anchor_sequence='GTACTGCGGCCGCTACCTA', quality_threshold=30):
        """
        Initialize the BarcodeExtractor with the input and output directories, anchor sequence, and quality threshold.
        
        :param input_directory: str, path to the input directory containing FASTQ files
        :param output_directory: str, path to the output directory where extracted barcode files will be saved
        :param anchor_sequence: str, the anchor sequence used to identify and extract barcodes (default: 'GTACTGCGGCCGCTACCTA')
        :param quality_threshold: int, the minimum average quality score for a barcode to be included (default: 30)
        """
        self.input_directory = input_directory
        self.output_directory = output_directory
        self.anchor_sequence = anchor_sequence
        self.quality_threshold = quality_threshold

    def get_average_quality(self, qualities):
        """
        Calculate the average quality score for a list of quality scores.
        
        :param qualities: list of int, quality scores for each base in a barcode
        :return: float, average quality score
        """
        return sum(qualities) / len(qualities)

    def filter_barcodes(self, barcode_qualities):
        """
        Filter barcodes based on their average quality scores.
        
        :param barcode_qualities: dict, a dictionary mapping barcodes to their list of quality scores
        :return: list of str, filtered barcodes that meet the quality threshold
        """
        return [barcode for barcode, qualities in barcode_qualities.items()
                if self.get_average_quality(qualities) >= self.quality_threshold]

    def extract_barcodes(self, fastq_file, is_gzipped=True): #Modify depending on use case
        """
        Extract barcodes from a FASTQ file using the anchor sequence.
        
        :param fastq_file: str, path to the FASTQ file to be processed
        :param is_gzipped: bool, set to True if the file is compressed with gzip
        :return: list of str, extracted and filtered barcodes
        """
        barcodes = []

        with (gzip.open(fastq_file, 'rt') if is_gzipped else open(fastq_file, 'rt')) as f:
            for record in SeqIO.parse(f, 'fastq'):
                sequence = str(record.seq)
                barcode_start = sequence.find(self.anchor_sequence)

                if barcode_start >= 30:
                    barcode = sequence[barcode_start-30:barcode_start]
                    barcode_qualities = record.letter_annotations['phred_quality'][barcode_start-30:barcode_start]
                    barcodes.append((barcode, barcode_qualities))

        filtered_barcodes = self.filter_barcodes({barcode: qualities for barcode, qualities in barcodes})
        return filtered_barcodes

    def process_fastq_files_helper(self, file_name):
        """
        Helper function to process a single FASTQ file.
        """
        input_file = os.path.join(self.input_directory, file_name)
        output_file = os.path.join(self.output_directory, file_name.split(".")[0] + "_barcodes.txt")

        # Read compressed or uncompressed files as appropriate
        if file_name.endswith(".gz"):
            barcodes = self.extract_barcodes(input_file, True)
        else:
            barcodes = self.extract_barcodes(input_file, False)
        with open(output_file, "w") as f:
            for barcode in barcodes:
                f.write(barcode + "\n")
                
        return barcodes

    def process_fastq_files(self):
        """
        Process all FASTQ files in the input directory, extracting barcodes and saving them to the output directory.
        """
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)

        all_barcodes = []  # Initializing empty list for all_barcodes

        fastq_files = [file_name for file_name in os.listdir(self.input_directory) if file_name.endswith(".fastq") or file_name.endswith(".fastq.gz")]
        
        num_processes= 10 #This number can be changed based on available computational resources. 
        
        with Pool(num_processes) as pool:
            all_barcodes_list = pool.map(self.process_fastq_files_helper, fastq_files)

       # Flatten the list of lists into a single list
        all_barcodes = [barcode for sublist in all_barcodes_list for barcode in sublist]
        
        ''' To profile the code: uncomment the following'''
# Uncomment the following lines to profile the code
# import cProfile, pstats
# profile = cProfile.Profile()
# profile.enable()
# main()
# profile.disable()
# with open("profile_stats.txt", "w") as f:
#     ps = pstats.Stats(profile, stream=f).sort_stats("cumulative")
#     ps.print_stats()

        return all_barcodes
                
       
