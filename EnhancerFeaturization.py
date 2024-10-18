#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 13:56:07 2024

@author: angela_gao
"""
##################################################################################################
################################# fasta Conversion Test ##########################################
##################################################################################################
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
        # Convert the line from dataframe to fasta format
        processed_line = '>' + data[0] + '\n' + data[1]        
        # Write the processed line to the output file
        outfile.write(processed_line + '\n')
        
##################################################################################################
##################################### fasta Conversion ###########################################
##################################################################################################
# Define the input filenames
SUD1_minP_DNA_InputFilename = '/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831773_ScaleUpDesign1_minP_Plasmid.counts.txt'
SUD1_SV40P_DNA_InputFilename = '/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831774_ScaleUpDesign1_SV40P_Plasmid.counts.txt'
SUD2_minP_DNA_InputFilename = '/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831775_ScaleUpDesign2_minP_Plasmid.counts.txt'
SUD2_SV40P_DNA_InputFilename = '/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831776_ScaleUpDesign2_SV40P_Plasmid.counts.txt'

# Define the output filenames
SUD1_minP_DNA_OutputFilename = '/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/DNA_fasta/GSM1831773_ScaleUpDesign1_minP_Plasmid_fasta.txt'
SUD1_SV40P_DNA_OutputFilename = '/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/DNA_fasta/GSM1831774_ScaleUpDesign1_SV40P_Plasmid_fasta.txt'
SUD2_minP_DNA_OutputFilename = '/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/DNA_fasta/GSM1831775_ScaleUpDesign2_minP_Plasmid_fasta.txt'
SUD2_SV40P_DNA_OutputFilename = '/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/DNA_fasta/GSM1831776_ScaleUpDesign2_SV40P_Plasmid_fasta.txt'

# Open the input file in read mode and the output file in write mode
############################ SUD1 minP DNA ####################################
with open(SUD1_minP_DNA_InputFilename, 'r') as infile, open(SUD1_minP_DNA_OutputFilename, 'w') as outfile:
    # Read the header line
    header = infile.readline().split('\t')
    
    # Read the remaining lines
    for line in infile:
        # Split the line into data
        data = line.split('\t')
        # Convert the line from dataframe to fasta format
        processed_line = '>' + data[0] + '\n' + data[1]        
        # Write the processed line to the output file
        outfile.write(processed_line + '\n')

############################ SUD1 SV40P DNA ####################################
with open(SUD1_SV40P_DNA_InputFilename, 'r') as infile, open(SUD1_SV40P_DNA_OutputFilename, 'w') as outfile:
    # Read the header line
    header = infile.readline().split('\t')
    
    # Read the remaining lines
    for line in infile:
        # Split the line into data
        data = line.split('\t')
        # Convert the line from dataframe to fasta format
        processed_line = '>' + data[0] + '\n' + data[1]        
        # Write the processed line to the output file
        outfile.write(processed_line + '\n')
        
############################ SUD2 minP DNA ####################################
with open(SUD2_minP_DNA_InputFilename, 'r') as infile, open(SUD2_minP_DNA_OutputFilename, 'w') as outfile:
    # Read the header line
    header = infile.readline().split('\t')
    
    # Read the remaining lines
    for line in infile:
        # Split the line into data
        data = line.split('\t')
        # Convert the line from dataframe to fasta format
        processed_line = '>' + data[0] + '\n' + data[1]        
        # Write the processed line to the output file
        outfile.write(processed_line + '\n')

############################ SUD2 SV40P DNA ####################################
with open(SUD2_SV40P_DNA_InputFilename, 'r') as infile, open(SUD2_SV40P_DNA_OutputFilename, 'w') as outfile:
    # Read the header line
    header = infile.readline().split('\t')
    
    # Read the remaining lines
    for line in infile:
        # Split the line into data
        data = line.split('\t')
        # Convert the line from dataframe to fasta format
        processed_line = '>' + data[0] + '\n' + data[1]        
        # Write the processed line to the output file
        outfile.write(processed_line + '\n')



