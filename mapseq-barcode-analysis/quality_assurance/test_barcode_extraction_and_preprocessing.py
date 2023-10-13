import sys
sys.path.append('/home/pau/micromamba/envs/cajalneuro/cajalneuro-main/mapseq-barcode-analysis') #Points to the file path containing the function BarcodeExtractor. The path will need to be modified per user. 
#The file fast_data_parsing.py is not in the same directory as this quality assurance script.

''' 
This test suite has three tests:

test_input_output_lengths: This test checks if the number of extracted barcodes is greater than the number of sequenced cells and within a reasonable margin of error.

test_barcode_and_anchor_sequence: This test ensures that each extracted barcode is 30 nucleotides long.

test_hamming_1_neighbors: This test checks if the vast majority (at least 95% in this example)
'''
import io
import os
import unittest
from fastq_data_parsing import BarcodeExtractor
from unittest.runner import TextTestRunner, TextTestResult

class CustomTextTestResult(TextTestResult):
    """
    CustomTextTestResult extends the TextTestResult class to store the test results
    in a list so that they can be accessed and displayed later.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.results = []

    def addSuccess(self, test):
        super().addSuccess(test)
        self.results.append(('pass', test))

    def addFailure(self, test, err):
        super().addFailure(test, err)
        self.results.append(('fail', test, err))

def hamming_distance(s1, s2):
    """
    Calculates the Hamming distance between two strings of equal length.
    """
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

class CustomTextTestRunner(TextTestRunner):
    """
    CustomTextTestRunner extends the TextTestRunner class to use the CustomTextTestResult
    class for storing test results.
    """
    resultclass = CustomTextTestResult

class TestBarcodeExtractor(unittest.TestCase):

    def setUp(self):
        """
        This function is called before each test. It initializes the BarcodeExtractor instance
        and processes the FASTQ files to extract barcodes, which are then stored in the
        'barcodes' list.
        """
        input_directory = "path/to/input_directory" #MAIN.PY OUTPUT DIRECTORY
        output_directory = "path/to/output_directory" #CURRENT WORKING DIRECTORY
        anchor_sequence = "GTACTGCGGCCGCTACCTA"

        self.barcode_extractor = BarcodeExtractor(input_directory, output_directory, anchor_sequence)
        self.barcode_extractor.process_fastq_files()

        self.barcodes = []
        for file_name in os.listdir(output_directory):
            if file_name.endswith("_barcodes.txt"):
                with open(os.path.join(output_directory, file_name), "r") as f:
                    self.barcodes.extend([line.strip() for line in f])

    def test_barcode_and_anchor_sequence(self):
        """
        This test ensures that each extracted barcodes are 30 nucleotides long.
        """   
        for barcode in self.barcodes:
            self.assertEqual(len(barcode), 30)

    def test_hamming_1_neighbors(self):
        """
        This test checks if the vast majority (at least 95% in this example)
        of the barcodes have a Hamming distance of 1 to at least one other barcode.
        """
        hamming_1_neighbors_count = 0
        for i, barcode1 in enumerate(self.barcodes):
            for barcode2 in self.barcodes[i + 1:]:
                if hamming_distance(barcode1, barcode2) == 1:
                    hamming_1_neighbors_count += 1

        self.assertGreater(hamming_1_neighbors_count / len(self.barcodes), 0.95)

    def test_input_output_lengths(self):
        '''
        This test checks if the number of extracted barcodes is greater than the number of sequenced cells and within a reasonable margin of error. (10% in this example)
        '''
        num_cells_sequenced = 100000
        self.assertGreater(len(self.barcodes), num_cells_sequenced)
        self.assertLess(len(self.barcodes), num_cells_sequenced * 1.1)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBarcodeExtractor)
    suite.addTest(TestBarcodeExtractor('test_input_output_lengths'))
    runner = CustomTextTestRunner(verbosity=0, stream=io.StringIO())
    result = runner.run(suite)
    
    print("==================== Test Results ====================")

    for i, (test_result, test_case, *extra) in enumerate(result.results, start=1):
        test_name = test_case._testMethodName
        print(f"{i}. {test_name}")
        if test_result == "pass":
            print("   Status: PASS\n")
        else:
            err = extra[0] if extra else None
            fail_reason = test_case.shortDescription()
            print(f"   Status: FAIL\n   Reason: {fail_reason}\n")
    
    print("=======================================================")


    if not result.wasSuccessful():
        print("Please review the data, processing steps, or thresholds, and try again.")
