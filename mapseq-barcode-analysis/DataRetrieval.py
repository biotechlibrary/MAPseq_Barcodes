'''This script retrieves real data from the National Center for Biotechnology Information'''
  
#I wrote this to take email + number of files to download as input parameters, assuming multiple people will be using this code for their own stuff.
#Instead of having to go in and change 'email' or 'num_of_files_to_download every time. 

#If it's easier for the user to simply edit the variable, then the code can be easily adjusted. 

import os
import subprocess
from Bio import Entrez
import requests

class MAPseqDataDownloader:
    def __init__(self, email, num_files_to_download):
        self.email = email
        self.num_files_to_download = num_files_to_download
        self.accession_ids = []
    
    # Search for MAPseq data in the SRA database
    def search_mapseq_data(self):
        Entrez.email = self.email
        query = "MAPseq"
        search_handle = Entrez.esearch(db="sra", term=query, retmax=self.num_files_to_download)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        self.accession_ids = search_results["IdList"]

    # Download SRA files from the NCBI server
    def download_sra_files(self, output_dir):
        for sra_id in self.accession_ids:
            url = f"https://sra-download.ncbi.nlm.nih.gov/traces/sra/sra-instant/reads/ByRun/sra/SRR/{sra_id[:6]}/{sra_id}/{sra_id}.sra"
            response = requests.get(url, stream=True)

            if response.status_code == 200:
                output_file = os.path.join(output_dir, f"{sra_id}.sra")
                with open(output_file, "wb") as file:
                    for chunk in response.iter_content(chunk_size=8192):
                        file.write(chunk)

    # Convert SRA files to fastq format using fastq-dump from the SRA Toolkit
    def convert_sra_to_fastq(self, sra_files, output_dir):
        for sra_file in sra_files:
            sra_path = os.path.join(output_dir, sra_file)
            fastq_path = os.path.splitext(sra_path)[0] + ".fastq"
            subprocess.run(["fastq-dump", sra_path, "-O", output_dir, "--gzip", "-Z", ">", fastq_path])

    # Run the entire data downloading and conversion process
    def run(self, output_dir):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        self.search_mapseq_data()
        self.download_sra_files(output_dir)

        sra_files = [f"{sra_id}.sra" for sra_id in self.accession_ids]
        self.convert_sra_to_fastq(sra_files, output_dir)

if __name__ == "__main__":
    # Replace with your actual email address
    email = "your_email@example.com"
    # Set the number of SRA files to download
    num_files_to_download = 5
    # Set the output directory for downloaded files and fastq conversions
    output_dir = "mapseq_data"

    # Create a downloader instance and run the process
    downloader = MAPseqDataDownloader(email, num_files_to_download)
    downloader.run(output_dir)

