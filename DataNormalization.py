# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 2024

@author: angela_gao
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

##################################################################################################
####################################### DATA INPUT ###############################################
##################################################################################################

HepG2_SUD1_minP_mRNA_Rep1 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831757_HepG2_ScaleUpDesign1_minP_mRNA_Rep1.counts.txt', sep='\t')
HepG2_SUD1_minP_mRNA_Rep2 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831758_HepG2_ScaleUpDesign1_minP_mRNA_Rep2.counts.txt', sep='\t')
SUD1_minP_DNA = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831773_ScaleUpDesign1_minP_Plasmid.counts.txt', sep='\t')

# Rename the 'Counts' Column
HepG2_SUD1_minP_mRNA_Rep1.rename(columns={'Counts': 'mRNACounts'}, inplace = True)
HepG2_SUD1_minP_mRNA_Rep2.rename(columns={'Counts': 'mRNACounts'}, inplace = True)
SUD1_minP_DNA.rename(columns={'Counts': 'DNACounts'}, inplace = True)

# Only keep ID tags and counts in the DNA file
SUD1_minP_DNA_short = SUD1_minP_DNA[['Tags', 'DNACounts']] # assuming DNA tags are unique

# Merge DNA and RNA counts using sequence tags CellInfo_short_Time1 = CellInfo_short1.merge(ImageInfo_short1, how="left", on= ["ImageNumber"])
HepG2_SUD1_minP_mRNADNA_Rep1 = HepG2_SUD1_minP_mRNA_Rep1.merge(SUD1_minP_DNA_short, how='left', on=['Tags'])
HepG2_SUD1_minP_mRNADNA_Rep2 = HepG2_SUD1_minP_mRNA_Rep2.merge(SUD1_minP_DNA_short, how='left', on=['Tags'])

##################################################################################################
##################################### DATA PROCESSING ############################################
##################################################################################################

# Compute the log2 fold change between two replicates, log2(RNA+1/DNA+1) = log2(RNA+1) - log2(DNA+1)
HepG2_SUD1_minP_mRNADNA_Rep1['log2FC'] = np.log2(HepG2_SUD1_minP_mRNADNA_Rep1['mRNACounts']+1) - np.log2(HepG2_SUD1_minP_mRNADNA_Rep1['DNACounts']+1) # Calculates log2 FC for replicate 1
HepG2_SUD1_minP_mRNADNA_Rep2['log2FC'] = np.log2(HepG2_SUD1_minP_mRNADNA_Rep2['mRNACounts']+1) - np.log2(HepG2_SUD1_minP_mRNADNA_Rep2['DNACounts']+1) # Calculates log2 FC for replicate 2

##################################################################################################
######################################### PLOTING ################################################
##################################################################################################
plt.rcParams['font.size'] = '40'
plt.rcParams.update({'font.family':'Helvetica'})
plt.figure(figsize=(15, 15))

plt.scatter(HepG2_SUD1_minP_mRNADNA_Rep1['log2FC'], HepG2_SUD1_minP_mRNADNA_Rep2['log2FC'], marker='H', s=600, alpha=0.1)
plt.xlabel('Log2 Fold Change, HepG2 minP Replicate 1')
plt.ylabel('Log2 Fold Change, HepG2 minP Replicate 2')
plt.xlim(-10.5, 10.5)
plt.ylim(-10.5, 10.5)