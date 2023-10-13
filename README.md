# **MAPseq Barcode Analysis**

This repository offers a comprehensive pipeline that facilitates the processing and analysis of barcode data obtained from MAPseq experiments. This pipeline is designed to enable users to download, generate, or use existing fastq sequences for extracting and pre-processing barcode sequences. The pipeline analyzes these apparent barcodes and returns a list of true underlying barcodes.

## **Overview**

>MAPseq is a technique for high-throughput mapping of connections in the brain. Neurons are virally induced to produce random RNA sequences ("barcodes"), which fill both their somas (cell bodies) as well as their axons (outgoing signaling structures). By analyzing the connections between these barcodes, we can learn about the structure and function of the brain.

**These scripts perform the following tasks**:
1. Parse input data (fastq format) to extract barcode sequences. 
2. Calculate the Hamming distance between barcode sequences. 
3. Cluster barcode sequences based on Hamming distance.
4. Filter out the true underlying barcodes based on the parameters provided.
5. Output the list of true underlying barcodes. 

## **Requirements**

To run this project, you need to have Python 3.x installed. Additionally, this project uses a Conda virtual environment called MAPseqRNA. To create this environment, follow the instructions below.

### **Setup**:
1. Install [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) if you haven't already.
2. Create a new Conda virtual environment called MAPseqRNA: 
```
conda create -n MAPseqRNA python=3.x
```
Replace '3.x' with your desired Python version (e.g., '3.8')

3. Activate the environment: 
```
conda activate MAPseqRNA
```
4. Download the repository as a zip file from GitHub. You can do this by clicking the green "Code" button on the main repository page and selecting "Download ZIP".

5. Extract the downloaded zip file inside your main environment directory.
```
unzip MAPseqRNA-main.zip
```
6. Set the required dependencies.
```
pip install -r requirements.txt
```
7. Navigate to the extracted project directory.
```
cd MAPseqRNA-main/mapseq-barcode-analysis
```
8. Edit PATH in bash: 

In a Linux machine, you can add the SRA Toolkit's bin directory to your PATH by modifying the .bashrc or .bash_profile file in your home directory. Add the following line to the file, replacing /path/to/sratoolkit/bin with the actual path to the bin directory in your SRA Toolkit installation:

```
export PATH=$PATH:/path/to/sratoolkit/bin
```
Save the file and restart your terminal, or run source ~/.bashrc or source ~/.bash_profile to reload the configuration. After updating the PATH, you should be able to run fastq-dump and other SRA Toolkit utilities from any directory in your terminal.

- **Other Dependencies**

This project requires the SRA Toolkit to be installed on your system. Please follow the instructions for your operating system on the official website to install it: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
```
sudo apt install sra-toolkit
```

Make sure that the `fastq-dump` command is in your system's PATH before running the project.

_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
## **USAGE**
1. Download the repository and navigate to the project directory

2. If you have your own fastq files, place them in the input directory specified in the main.py script. If you want to download fastq files as part of the pipeline, set the appropriate parameters in the main.py script, such as your email address and download limit. Note that both email and download limit are required if you want to download fastq files.

3. Run the pipeline by executing the main.py script:
```
python3 main.py [path/to/user_provided_data]
```
If you have your own fastq files, provide the path to the directory containing them as a command-line argument. If you want to download fastq files, do not provide the command-line argument.

4. Optional: If you want to change the Hamming distance threshold, open the main.py script and modify the "hamming_distance_threshold" variable in the script. By default, it is set to 1.

**The pipeline will use your provided fastq files or download them if specified, extract and preprocess barcode sequences, and analyze the barcodes to generate a list of true underlying barcodes. The output file containing the true barcodes will be saved in the specified output directory.**
________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
## **Preprocessing and Quality Assurance**

- **fastq_data_parsing.py**

The 'fastq_data_parsing.py' script processes FASTQ files, extracts barcodes based on the anchor sequence, and saves the extracted barcodes to an output directory. To use this script:
1. Run the script using a Python interpreter:
```
python3 fastq_data_parsing.py
```
###**Barcode Statistics and Validation**

This test suite contains two scripts: analyze_barcodes.py and test_barcode_extraction_and_preprocessing.py. They perform various tasks to analyze and test the extracted barcodes.

- **analyze_barcodes.py**

This script reads the extracted barcodes, calculates basic statistics, and checks their length and format. 

To use this script:
1. Update the 'input_directory' variable in the 'main()' function to the path of the directory containing the output files generated by the 'fastq_data_parsing.py' script.
2. Run the script using a Python interpreter:
```
python3 analyze_barcodes.py
```

- **test_barcode_extraction_and_preprocessing.py**

This script contains a test suite that tests the barcode extraction and preprocessing using the 'BarcodeExtractor' class. 

The test suite includes three tests:
1. text_input_output_lengths: This test checks if the number of extracted barcodes is greater than the number of sequenced cells and within a reasonable margin of error.
2. test_barcode_and_anchor_sequence: This test ensures that each extracted barcode is 30 nucleotides long. 
3. test_hamming_1_neighbors: This test checks if the vast majority (at least 95% in this example) of barcodes are Hamming-1 neighbors. 

To use this test suite:
1. Update the 'input_directory' and 'output_directory' variables in the 'setUp()' function to the appropriate paths for the input and output directories.
2. Run the test suite using a Python interpreter: 
```
python3 test_barcpde_extraction_and_preprocessing.py
```

The test suite will automatically run all three tests and displays the results.

### **Profiling**
The profiling feature is useful for users who want to evaluate the performance of the pipeline and identify potential bottlenecks or areas for optimization. By profiling the code, users can gain insight into which parts of the pipeline take the most time and require further improvement.

To profile the code, uncomment the following lines at the end of the script:
```
# Uncomment the following lines to profile the code
# import cProfile, pstats
# profile = cProfile.Profile()
# profile.enable()
# run_pipeline(email, download_limit, input_directory, output_directory, user_provided_data)
# profile.disable()
# with open("profile_output.txt", "w") as f:
#     ps = pstats.Stats(profile, stream=f).sort_stats("cumulative")
#     ps.print_stats()
```
The profiling results will be saved in a file called "profile_stats.txt". Review this file to identify bottlenecks and optimize the code. 

# **Pipeline Overview**

The pipeline is composed of several scripts and classes:

- main.py: The main script that ties all the classes together, runs the pipeline, and coordinates the downloading, preprocessing, and analyzing of barcode data.

- fastq_downloader.py: Contains the FastqDownloader class for downloading fastq files from the NCBI SRA database.

- fastq_data_parsing.py: Contains the BarcodeExtractor class for extracting and preprocessing barcode sequences from fastq files.

- mapseq_barcode_analysis.py: Contains the MAPseqBarcodeAnalysis class for analyzing barcode sequences and generating a list of true underlying barcodes based on Hamming distance.

- test_barcode_extraction_and_preprocessing.py: Contains unit tests for validating barcode extraction and preprocessing steps.

- requirements.txt: Lists the required Python packages for the pipeline.

- README.md: This file, containing an overview and instructions for the pipeline.
