import os
import sys
import cProfile
import argparse
import pstats
from tqdm import tqdm
from DataRetrieval import MAPseqDataDownloader
from fastq_data_parsing import BarcodeExtractor
from mapseq_barcode_analysis import MAPseqBarcodeAnalysis

def run_pipeline(email, download_limit, input_directory, output_directory, hamming_distance_threshold=1, user_provided_data=None):
    """
    This function runs the entire pipeline to download, preprocess, and analyze barcode data.
    
    Args:
        email (str, optional): User email address for NCBI Entrez API access. Defaults to None.
        download_limit (int, optional): Maximum number of fastq files to download. Defaults to None.
        input_directory (str): Path to the input directory containing fastq files.
        output_directory (str): Path to the output directory where results will be stored.
        hamming_distance_threshold (int, optional): Maximum Hamming distance to group similar barcodes. Defaults to 1.
        user_provided_data (str, optional): Path to user-provided data (fastq files). Defaults to None.

    Returns:
        None
    """

    # Wrap tqdm around the pipeline steps list
    pipeline_steps = [
        "Downloading fastq files",
        "Extracting and preprocessing barcodes",
        "Analyzing barcodes"
    ]

    for step in tqdm(pipeline_steps, desc="Running pipeline"):
        print("Retrieving files...")
        if step == "Downloading fastq files":
            if not user_provided_data:
                if email is not None and download_limit is not None:
                    fastq_downloader = MAPseqDataDownloader(email, download_limit)
                    fastq_downloader.run(input_directory)
                else:
                    raise ValueError("Email and download_limit are required for downloading fastq files.")
            else:
                input_directory = user_provided_data

        elif step == "Extracting and preprocessing barcodes":
            print("Extracting barcodes...")
            barcode_extractor = BarcodeExtractor(input_directory, output_directory)
            extracted_barcodes = barcode_extractor.process_fastq_files()

        elif step == "Analyzing barcodes":
            print("Validating barcodes...")
            mapseq_analyzer = MAPseqBarcodeAnalysis(extracted_barcodes, hamming_distance_threshold)
            true_barcodes = mapseq_analyzer.get_true_underlying_barcodes(extracted_barcodes)

            with open(os.path.join(output_directory, "true_barcodes.txt"), "w") as f:
                for barcode in true_barcodes:
                    f.write(barcode + "\n")

if __name__ == "__main__":
    email = None
    download_limit = None
    input_directory = sys.argv[1]
    output_directory = "output"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    # The default hamming_distance_threshold is set to 1. Change this value if needed.
    hamming_distance_threshold = 1
    
    user_provided_data = None
    if len(sys.argv) > 1:
        user_provided_data = sys.argv[1]

    profile = cProfile.Profile()
    profile.enable()
    run_pipeline(email, download_limit, input_directory, output_directory, hamming_distance_threshold, user_provided_data)
    profile.disable()

    with open("profile_output.txt", "w") as f:
        ps = pstats.Stats(profile, stream=f).sort_stats("cumulative")
        ps.print_stats()
