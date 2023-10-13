'''
This program has the following objectives:

1. Develop a function to compute the Hamming distance between barcodes.
2. Create a function that groups similar barcodes based on their Hamming distance. 
3. For each group, find the most likely real barcode based on the frequency of each barcode.
4. Consolidate the most likely real barcodes from each group to create a final list of true underlying barcodes. 
'''

import os
import gzip
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict, Counter
from Bio import SeqIO
import numpy as np
from collections import deque

class MAPseqBarcodeAnalysis:

    def __init__(self, input_directory, output_directory, anchor_sequence='GTACTGCGGCCGCTACCTA'):
        self.input_directory = input_directory
        self.output_directory = output_directory
        self.anchor_sequence = anchor_sequence

    def hamming_distance(self, barcode1, barcode2):
        '''
        Computes the Hamming distance between two barcodes.

        :param barcode1: str, first barcode sequence
        :param barcode2: str, second barcode sequence
        :return: int, Hamming distance between the two barcodes
        '''
        return sum(base1 != base2 for base1, base2 in zip(barcode1, barcode2))

    def group_similar_barcodes(self, barcodes, max_hamming_distance=1):
        """
        Groups similar barcodes based on their Hamming distance.

        :param barcodes: list of str, barcode sequences
        :param max_hamming_distance: int, maximum Hamming distance to consider barcodes as similar (default: 1)
        :return: list of lists, grouped barcodes
        """
        def find_neighbors(barcode, barcodes):
            """
            Finds neighbors of a given barcode within the maximum Hamming distance.

            :param barcode: str, barcode sequence
            :param barcodes: list of str, barcode sequences
            :return: list of str, neighboring barcodes
            """
            neighbors = []
            for other_barcode in barcodes:
                if self.hamming_distance(barcode, other_barcode) <= max_hamming_distance:
                    neighbors.append(other_barcode)
            return neighbors
        
        def dfs(visited, graph, node):
            """
            Depth-first search to traverse the barcode similarity graph.

            :param visited: set, visited barcodes
            :param graph: dict, adjacency list of the barcode similarity graph
            :param node: str, current barcode node
            """
            if node not in visited:
                visited.add(node)
                for neighbor in graph[node]:
                    dfs(visited, graph, neighbor)

        graph = {barcode: find_neighbors(barcode, barcodes) for barcode in barcodes}
        visited = set()
        barcode_groups = []

        for barcode in barcodes:
            if barcode not in visited:
                connected_component = set()
                dfs(connected_component, graph, barcode)
                barcode_groups.append(list(connected_component))

        return barcode_groups
    
    def find_most_likely_barcodes(self, barcode_groups):
        """
        Finds the most likely real barcode in each group based on the frequency of each barcode.

        :param barcode_groups: list of lists, grouped barcodes
        :return: list of str, most likely real barcodes
        """
        most_likely_barcodes = []
        for group in barcode_groups:
            counter = Counter(group)
            most_common_barcode, _ = counter.most_common(1)[0]
            most_likely_barcodes.append(most_common_barcode)
        return most_likely_barcodes

    def get_true_underlying_barcodes(self, barcodes):
        """
        Consolidates the most likely real barcodes from each group to create a final list of true underlying barcodes.

        :param barcodes: list of str, barcode sequences
        :return: list of str, true underlying barcodes
        """
        barcode_groups = self.group_similar_barcodes(barcodes)
        true_barcodes = self.find_most_likely_barcodes(barcode_groups)

        return true_barcodes

    
