#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 13:56:07 2024

@author: angela_gao
"""
# import pandas as pd

# SUD1_minP_DNA = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/SUD1_minP_PlasmidTest.txt', sep='\t')

# SUD1_minP_DNA.

# Define the input and output filenames
input_filename = '/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/Test/SUD1_minP_PlasmidTest.txt'
output_filename = '/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/Test/SUD1_minP_PlasmidTest_fasta.txt'

# Open the input file in read mode and the output file in write mode
with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
    # Read the header line
    header = infile.readline().split('\t')
    
    # Read the remaining lines
    for line in infile:
        # Split the line into data
        data = line.split('\t')
        # Process the line if needed (e.g., strip whitespace)
        processed_line = '>' + data[0] + '\n' + data[1]
        
        # Write the processed line to the output file
        outfile.write(processed_line + '\n')