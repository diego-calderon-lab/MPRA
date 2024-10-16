# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 2024

@author: angela_gao
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

##################################################################################################
####################################### DATA INPUT ###############################################
##################################################################################################

############################## Scaled Up Design 1 #############################
################################ Read in DNA data #############################
SUD1_minP_DNA = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831773_ScaleUpDesign1_minP_Plasmid.counts.txt', sep='\t')
SUD1_SV40P_DNA = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831774_ScaleUpDesign1_SV40P_Plasmid.counts.txt', sep='\t')
SUD1_minP_DNA.rename(columns={'Counts': 'DNACounts'}, inplace = True)
SUD1_minP_DNA_short = SUD1_minP_DNA[['Tags', 'DNACounts']] # assuming DNA tags are unique
SUD1_SV40P_DNA.rename(columns={'Counts': 'DNACounts'}, inplace = True)
SUD1_SV40P_DNA_short = SUD1_SV40P_DNA[['Tags', 'DNACounts']] # assuming DNA tags are unique

###################### Read in HepG2 SUD1 minP data ###########################
HepG2_SUD1_minP_mRNA_Rep1 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831757_HepG2_ScaleUpDesign1_minP_mRNA_Rep1.counts.txt', sep='\t')
HepG2_SUD1_minP_mRNA_Rep2 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831758_HepG2_ScaleUpDesign1_minP_mRNA_Rep2.counts.txt', sep='\t')

# Rename the 'Counts' Column
HepG2_SUD1_minP_mRNA_Rep1.rename(columns={'Counts': 'mRNACounts'}, inplace = True)
HepG2_SUD1_minP_mRNA_Rep2.rename(columns={'Counts': 'mRNACounts'}, inplace = True)

# Merge DNA and RNA counts using sequence tags CellInfo_short_Time1 = CellInfo_short1.merge(ImageInfo_short1, how="left", on= ["ImageNumber"])
HepG2_SUD1_minP_mRNADNA_Rep1 = HepG2_SUD1_minP_mRNA_Rep1.merge(SUD1_minP_DNA_short, how='left', on=['Tags'])
HepG2_SUD1_minP_mRNADNA_Rep2 = HepG2_SUD1_minP_mRNA_Rep2.merge(SUD1_minP_DNA_short, how='left', on=['Tags'])

###################### Read in HepG2 SUD1 SV40P data ##########################
HepG2_SUD1_SV40P_mRNA_Rep1 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831759_HepG2_ScaleUpDesign1_SV40P_mRNA_Rep1.counts.txt', sep='\t')
HepG2_SUD1_SV40P_mRNA_Rep2 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831760_HepG2_ScaleUpDesign1_SV40P_mRNA_Rep2.counts.txt', sep='\t')

# Rename the 'Counts' Column
HepG2_SUD1_SV40P_mRNA_Rep1.rename(columns={'Counts': 'mRNACounts'}, inplace = True)
HepG2_SUD1_SV40P_mRNA_Rep2.rename(columns={'Counts': 'mRNACounts'}, inplace = True)

# Merge DNA and RNA counts using sequence tags CellInfo_short_Time1 = CellInfo_short1.merge(ImageInfo_short1, how="left", on= ["ImageNumber"])
HepG2_SUD1_SV40P_mRNADNA_Rep1 = HepG2_SUD1_SV40P_mRNA_Rep1.merge(SUD1_SV40P_DNA_short, how='left', on=['Tags'])
HepG2_SUD1_SV40P_mRNADNA_Rep2 = HepG2_SUD1_SV40P_mRNA_Rep2.merge(SUD1_SV40P_DNA_short, how='left', on=['Tags'])

###################### Read in K562 SUD1 SV40P data ###########################
K562_SUD1_SV40P_mRNA_Rep1 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831767_K562_ScaleUpDesign1_SV40P_mRNA_Rep1.counts.txt', sep='\t')
K562_SUD1_SV40P_mRNA_Rep2 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831768_K562_ScaleUpDesign1_SV40P_mRNA_Rep2.counts.txt', sep='\t')

# Rename the 'Counts' Column
K562_SUD1_SV40P_mRNA_Rep1.rename(columns={'Counts': 'mRNACounts'}, inplace = True)
K562_SUD1_SV40P_mRNA_Rep2.rename(columns={'Counts': 'mRNACounts'}, inplace = True)

# Merge DNA and RNA counts using sequence tags CellInfo_short_Time1 = CellInfo_short1.merge(ImageInfo_short1, how="left", on= ["ImageNumber"])
K562_SUD1_SV40P_mRNADNA_Rep1 = K562_SUD1_SV40P_mRNA_Rep1.merge(SUD1_SV40P_DNA_short, how='left', on=['Tags'])
K562_SUD1_SV40P_mRNADNA_Rep2 = K562_SUD1_SV40P_mRNA_Rep2.merge(SUD1_SV40P_DNA_short, how='left', on=['Tags'])

###################### Read in K562 SUD1 minP data ###########################
K562_SUD1_minP_mRNA_Rep1 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831765_K562_ScaleUpDesign1_minP_mRNA_Rep1.counts.txt', sep='\t')
K562_SUD1_minP_mRNA_Rep2 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831766_K562_ScaleUpDesign1_minP_mRNA_Rep2.counts.txt', sep='\t')

# Rename the 'Counts' Column
K562_SUD1_minP_mRNA_Rep1.rename(columns={'Counts': 'mRNACounts'}, inplace = True)
K562_SUD1_minP_mRNA_Rep2.rename(columns={'Counts': 'mRNACounts'}, inplace = True)

# Merge DNA and RNA counts using sequence tags CellInfo_short_Time1 = CellInfo_short1.merge(ImageInfo_short1, how="left", on= ["ImageNumber"])
K562_SUD1_minP_mRNADNA_Rep1 = K562_SUD1_minP_mRNA_Rep1.merge(SUD1_minP_DNA_short, how='left', on=['Tags'])
K562_SUD1_minP_mRNADNA_Rep2 = K562_SUD1_minP_mRNA_Rep2.merge(SUD1_minP_DNA_short, how='left', on=['Tags'])

############################## Scaled Up Design 2 #############################
############################### Read in DNA data ##############################
SUD2_minP_DNA = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831775_ScaleUpDesign2_minP_Plasmid.counts.txt', sep='\t')
SUD2_SV40P_DNA = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831776_ScaleUpDesign2_SV40P_Plasmid.counts.txt', sep='\t')
SUD2_minP_DNA.rename(columns={'Counts': 'DNACounts'}, inplace = True)
SUD2_minP_DNA_short = SUD2_minP_DNA[['Tags', 'DNACounts']] # assuming DNA tags are unique
SUD2_SV40P_DNA.rename(columns={'Counts': 'DNACounts'}, inplace = True)
SUD2_SV40P_DNA_short = SUD2_SV40P_DNA[['Tags', 'DNACounts']] # assuming DNA tags are unique

###################### Read in HepG2 SUD2 minP data ###########################
HepG2_SUD2_minP_mRNA_Rep1 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831761_HepG2_ScaleUpDesign2_minP_mRNA_Rep1.counts.txt', sep='\t')
HepG2_SUD2_minP_mRNA_Rep2 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831762_HepG2_ScaleUpDesign2_minP_mRNA_Rep2.counts.txt', sep='\t')

# Rename the 'Counts' Column
HepG2_SUD2_minP_mRNA_Rep1.rename(columns={'Counts': 'mRNACounts'}, inplace = True)
HepG2_SUD2_minP_mRNA_Rep2.rename(columns={'Counts': 'mRNACounts'}, inplace = True)

# Merge DNA and RNA counts using sequence tags CellInfo_short_Time1 = CellInfo_short1.merge(ImageInfo_short1, how="left", on= ["ImageNumber"])
HepG2_SUD2_minP_mRNADNA_Rep1 = HepG2_SUD2_minP_mRNA_Rep1.merge(SUD2_minP_DNA_short, how='left', on=['Tags'])
HepG2_SUD2_minP_mRNADNA_Rep2 = HepG2_SUD2_minP_mRNA_Rep2.merge(SUD2_minP_DNA_short, how='left', on=['Tags'])

###################### Read in HepG2 SUD2 SV40P data ###########################
HepG2_SUD2_SV40P_mRNA_Rep1 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831763_HepG2_ScaleUpDesign2_SV40P_mRNA_Rep1.counts.txt', sep='\t')
HepG2_SUD2_SV40P_mRNA_Rep2 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831764_HepG2_ScaleUpDesign2_SV40P_mRNA_Rep2.counts.txt', sep='\t')

# Rename the 'Counts' Column
HepG2_SUD2_SV40P_mRNA_Rep1.rename(columns={'Counts': 'mRNACounts'}, inplace = True)
HepG2_SUD2_SV40P_mRNA_Rep2.rename(columns={'Counts': 'mRNACounts'}, inplace = True)

# Merge DNA and RNA counts using sequence tags CellInfo_short_Time1 = CellInfo_short1.merge(ImageInfo_short1, how="left", on= ["ImageNumber"])
HepG2_SUD2_SV40P_mRNADNA_Rep1 = HepG2_SUD2_SV40P_mRNA_Rep1.merge(SUD2_SV40P_DNA_short, how='left', on=['Tags'])
HepG2_SUD2_SV40P_mRNADNA_Rep2 = HepG2_SUD2_SV40P_mRNA_Rep2.merge(SUD2_SV40P_DNA_short, how='left', on=['Tags'])

###################### Read in K562 SUD2 SV40P data ###########################
K562_SUD2_SV40P_mRNA_Rep1 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831771_K562_ScaleUpDesign2_SV40P_mRNA_Rep1.counts.txt', sep='\t')
K562_SUD2_SV40P_mRNA_Rep2 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831772_K562_ScaleUpDesign2_SV40P_mRNA_Rep2.counts.txt', sep='\t')

# Rename the 'Counts' Column
K562_SUD2_SV40P_mRNA_Rep1.rename(columns={'Counts': 'mRNACounts'}, inplace = True)
K562_SUD2_SV40P_mRNA_Rep2.rename(columns={'Counts': 'mRNACounts'}, inplace = True)

# Merge DNA and RNA counts using sequence tags CellInfo_short_Time1 = CellInfo_short1.merge(ImageInfo_short1, how="left", on= ["ImageNumber"])
K562_SUD2_SV40P_mRNADNA_Rep1 = K562_SUD2_SV40P_mRNA_Rep1.merge(SUD2_SV40P_DNA_short, how='left', on=['Tags'])
K562_SUD2_SV40P_mRNADNA_Rep2 = K562_SUD2_SV40P_mRNA_Rep2.merge(SUD2_SV40P_DNA_short, how='left', on=['Tags'])

###################### Read in K562 SUD2 minP data ###########################
K562_SUD2_minP_mRNA_Rep1 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831769_K562_ScaleUpDesign2_minP_mRNA_Rep1.counts.txt', sep='\t')
K562_SUD2_minP_mRNA_Rep2 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831770_K562_ScaleUpDesign2_minP_mRNA_Rep2.counts.txt', sep='\t')

# Rename the 'Counts' Column
K562_SUD2_minP_mRNA_Rep1.rename(columns={'Counts': 'mRNACounts'}, inplace = True)
K562_SUD2_minP_mRNA_Rep2.rename(columns={'Counts': 'mRNACounts'}, inplace = True)

# Merge DNA and RNA counts using sequence tags CellInfo_short_Time1 = CellInfo_short1.merge(ImageInfo_short1, how="left", on= ["ImageNumber"])
K562_SUD2_minP_mRNADNA_Rep1 = K562_SUD2_minP_mRNA_Rep1.merge(SUD2_minP_DNA_short, how='left', on=['Tags'])
K562_SUD2_minP_mRNADNA_Rep2 = K562_SUD2_minP_mRNA_Rep2.merge(SUD2_minP_DNA_short, how='left', on=['Tags'])

##################################################################################################
##################################### DATA PROCESSING ############################################
##################################################################################################

############################## Scaled Up Design 1 #############################
# Compute the log2 fold change between two replicates, log2(RNA+1/DNA+1) = log2(RNA+1) - log2(DNA+1)
HepG2_SUD1_minP_mRNADNA_Rep1['log2FC'] = np.log2(HepG2_SUD1_minP_mRNADNA_Rep1['mRNACounts']+1) - np.log2(HepG2_SUD1_minP_mRNADNA_Rep1['DNACounts']+1) # Calculates log2 FC for replicate 1
HepG2_SUD1_minP_mRNADNA_Rep2['log2FC'] = np.log2(HepG2_SUD1_minP_mRNADNA_Rep2['mRNACounts']+1) - np.log2(HepG2_SUD1_minP_mRNADNA_Rep2['DNACounts']+1) # Calculates log2 FC for replicate 2

# Compute the log2 fold change between two replicates, log2(RNA+1/DNA+1) = log2(RNA+1) - log2(DNA+1)
HepG2_SUD1_SV40P_mRNADNA_Rep1['log2FC'] = np.log2(HepG2_SUD1_SV40P_mRNADNA_Rep1['mRNACounts']+1) - np.log2(HepG2_SUD1_SV40P_mRNADNA_Rep1['DNACounts']+1) # Calculates log2 FC for replicate 1
HepG2_SUD1_SV40P_mRNADNA_Rep2['log2FC'] = np.log2(HepG2_SUD1_SV40P_mRNADNA_Rep2['mRNACounts']+1) - np.log2(HepG2_SUD1_SV40P_mRNADNA_Rep2['DNACounts']+1) # Calculates log2 FC for replicate 2

# Compute the log2 fold change between two replicates, log2(RNA+1/DNA+1) = log2(RNA+1) - log2(DNA+1)
K562_SUD1_SV40P_mRNADNA_Rep1['log2FC'] = np.log2(K562_SUD1_SV40P_mRNADNA_Rep1['mRNACounts']+1) - np.log2(K562_SUD1_SV40P_mRNADNA_Rep1['DNACounts']+1) # Calculates log2 FC for replicate 1
K562_SUD1_SV40P_mRNADNA_Rep2['log2FC'] = np.log2(K562_SUD1_SV40P_mRNADNA_Rep2['mRNACounts']+1) - np.log2(K562_SUD1_SV40P_mRNADNA_Rep2['DNACounts']+1) # Calculates log2 FC for replicate 2

# Compute the log2 fold change between two replicates, log2(RNA+1/DNA+1) = log2(RNA+1) - log2(DNA+1)
K562_SUD1_minP_mRNADNA_Rep1['log2FC'] = np.log2(K562_SUD1_minP_mRNADNA_Rep1['mRNACounts']+1) - np.log2(K562_SUD1_minP_mRNADNA_Rep1['DNACounts']+1) # Calculates log2 FC for replicate 1
K562_SUD1_minP_mRNADNA_Rep2['log2FC'] = np.log2(K562_SUD1_minP_mRNADNA_Rep2['mRNACounts']+1) - np.log2(K562_SUD1_minP_mRNADNA_Rep2['DNACounts']+1) # Calculates log2 FC for replicate 2

############################## Scaled Up Design 2 #############################
# Compute the log2 fold change between two replicates, log2(RNA+1/DNA+1) = log2(RNA+1) - log2(DNA+1)
HepG2_SUD2_minP_mRNADNA_Rep1['log2FC'] = np.log2(HepG2_SUD2_minP_mRNADNA_Rep1['mRNACounts']+1) - np.log2(HepG2_SUD2_minP_mRNADNA_Rep1['DNACounts']+1) # Calculates log2 FC for replicate 1
HepG2_SUD2_minP_mRNADNA_Rep2['log2FC'] = np.log2(HepG2_SUD2_minP_mRNADNA_Rep2['mRNACounts']+1) - np.log2(HepG2_SUD2_minP_mRNADNA_Rep2['DNACounts']+1) # Calculates log2 FC for replicate 2

# Compute the log2 fold change between two replicates, log2(RNA+1/DNA+1) = log2(RNA+1) - log2(DNA+1)
HepG2_SUD2_SV40P_mRNADNA_Rep1['log2FC'] = np.log2(HepG2_SUD2_SV40P_mRNADNA_Rep1['mRNACounts']+1) - np.log2(HepG2_SUD2_SV40P_mRNADNA_Rep1['DNACounts']+1) # Calculates log2 FC for replicate 1
HepG2_SUD2_SV40P_mRNADNA_Rep2['log2FC'] = np.log2(HepG2_SUD2_SV40P_mRNADNA_Rep2['mRNACounts']+1) - np.log2(HepG2_SUD2_SV40P_mRNADNA_Rep2['DNACounts']+1) # Calculates log2 FC for replicate 2

# Compute the log2 fold change between two replicates, log2(RNA+1/DNA+1) = log2(RNA+1) - log2(DNA+1)
K562_SUD2_SV40P_mRNADNA_Rep1['log2FC'] = np.log2(K562_SUD2_SV40P_mRNADNA_Rep1['mRNACounts']+1) - np.log2(K562_SUD2_SV40P_mRNADNA_Rep1['DNACounts']+1) # Calculates log2 FC for replicate 1
K562_SUD2_SV40P_mRNADNA_Rep2['log2FC'] = np.log2(K562_SUD2_SV40P_mRNADNA_Rep2['mRNACounts']+1) - np.log2(K562_SUD2_SV40P_mRNADNA_Rep2['DNACounts']+1) # Calculates log2 FC for replicate 2

# Compute the log2 fold change between two replicates, log2(RNA+1/DNA+1) = log2(RNA+1) - log2(DNA+1)
K562_SUD2_minP_mRNADNA_Rep1['log2FC'] = np.log2(K562_SUD2_minP_mRNADNA_Rep1['mRNACounts']+1) - np.log2(K562_SUD2_minP_mRNADNA_Rep1['DNACounts']+1) # Calculates log2 FC for replicate 1
K562_SUD2_minP_mRNADNA_Rep2['log2FC'] = np.log2(K562_SUD2_minP_mRNADNA_Rep2['mRNACounts']+1) - np.log2(K562_SUD2_minP_mRNADNA_Rep2['DNACounts']+1) # Calculates log2 FC for replicate 2

##################################################################################################
######################################### PLOTING ################################################
##################################################################################################
plt.rcParams['font.size'] = '40' # Set default font size
plt.rcParams.update({'font.family':'Helvetica'}) # Set default font
plt.rcParams['figure.figsize'] = (18, 15) # Set default figure size

############################## Scaled Up Design 1 #############################
########################## Fig1 HepG2 SUD1 minP HexBin ########################
# Create a hexagonal bin plot
plt.figure(1)
plt.hexbin(HepG2_SUD1_minP_mRNADNA_Rep1['log2FC'], HepG2_SUD1_minP_mRNADNA_Rep2['log2FC'], gridsize=30, cmap='Blues', mincnt=1, bins='log')
# Add a color bar
plt.colorbar(label='Counts in bin')
plt.xlabel('Log2 Fold Change, HepG2 minP Replicate 1')
plt.ylabel('Log2 Fold Change, HepG2 minP Replicate 2')
plt.xlim(-10.5, 10.5)
plt.ylim(-10.5, 10.5)

# Calculate Spearman correlation coefficient
coef_HepG2_SUD1_minP, p_value_HepG2_SUD1_minP = spearmanr(HepG2_SUD1_minP_mRNADNA_Rep1['log2FC'], HepG2_SUD1_minP_mRNADNA_Rep2['log2FC'])

# Annotate the plot with the Spearman coefficient
plt.annotate(f'Spearman coefficient: {coef_HepG2_SUD1_minP:.3f}\nP-value: {p_value_HepG2_SUD1_minP:.1f}', 
              xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top', 
              bbox=dict(boxstyle='round', facecolor='white', alpha=0))
# Save the plot
plt.savefig('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA/MPRA_plots/HepG2_SUD1_minP_HexBin.pdf', format='pdf', bbox_inches='tight', transparent=True)

######################### Fig2 HepG2 SUD1 SV40P HexBin ########################
# Create a hexagonal bin plot
plt.figure(2)
plt.hexbin(HepG2_SUD1_SV40P_mRNADNA_Rep1['log2FC'], HepG2_SUD1_SV40P_mRNADNA_Rep2['log2FC'], gridsize=30, cmap='Blues', mincnt=1, bins='log')
# Add a color bar
plt.colorbar(label='Counts in bin')
plt.xlabel('Log2 Fold Change, HepG2 SV40P Replicate 1')
plt.ylabel('Log2 Fold Change, HepG2 SV40P Replicate 2')
plt.xlim(-10.5, 10.5)
plt.ylim(-10.5, 10.5)

# Calculate Spearman correlation coefficient
coef_HepG2_SUD1_SV40P, p_value_HepG2_SUD1_SV40P = spearmanr(HepG2_SUD1_SV40P_mRNADNA_Rep1['log2FC'], HepG2_SUD1_SV40P_mRNADNA_Rep2['log2FC'])

# Annotate the plot with the Spearman coefficient
plt.annotate(f'Spearman coefficient: {coef_HepG2_SUD1_SV40P:.3f}\nP-value: {p_value_HepG2_SUD1_SV40P:.1f}', 
              xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top', 
              bbox=dict(boxstyle='round', facecolor='white', alpha=0))
# Save the plot
plt.savefig('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA/MPRA_plots/HepG2_SUD1_SV40P_HexBin.pdf', format='pdf', bbox_inches='tight', transparent=True)

########################### Fig3 K562 SUD1 SV40P HexBin #######################
# Create a hexagonal bin plot
plt.figure(3)
plt.hexbin(K562_SUD1_SV40P_mRNADNA_Rep1['log2FC'], K562_SUD1_SV40P_mRNADNA_Rep2['log2FC'], gridsize=30, cmap='Blues', mincnt=1, bins='log')
# Add a color bar
plt.colorbar(label='Counts in bin')
plt.xlabel('Log2 Fold Change, K562 SV40P Replicate 1')
plt.ylabel('Log2 Fold Change, K562 SV40P Replicate 2')
plt.xlim(-10.5, 10.5)
plt.ylim(-10.5, 10.5)

# Calculate Spearman correlation coefficient
coef_K562_SUD1_SV40P, p_value_K562_SUD1_SV40P = spearmanr(K562_SUD1_SV40P_mRNADNA_Rep1['log2FC'], K562_SUD1_SV40P_mRNADNA_Rep2['log2FC'])

# Annotate the plot with the Spearman coefficient
plt.annotate(f'Spearman coefficient: {coef_K562_SUD1_SV40P:.3f}\nP-value: {p_value_K562_SUD1_SV40P:.1f}', 
              xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top', 
              bbox=dict(boxstyle='round', facecolor='white', alpha=0))
# Save the plot
plt.savefig('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA/MPRA_plots/K562_SUD1_SV40P_HexBin.pdf', format='pdf', bbox_inches='tight', transparent=True)

########################### Fig4 K562 SUD1 minP HexBin ########################
# Create a hexagonal bin plot
plt.figure(4)
plt.hexbin(K562_SUD1_minP_mRNADNA_Rep1['log2FC'], K562_SUD1_minP_mRNADNA_Rep2['log2FC'], gridsize=30, cmap='Blues', mincnt=1, bins='log')
# Add a color bar
plt.colorbar(label='Counts in bin')
plt.xlabel('Log2 Fold Change, K562 minP Replicate 1')
plt.ylabel('Log2 Fold Change, K562 minP Replicate 2')
plt.xlim(-10.5, 10.5)
plt.ylim(-10.5, 10.5)

# Calculate Spearman correlation coefficient
coef_K562_SUD1_minP, p_value_K562_SUD1_minP = spearmanr(K562_SUD1_minP_mRNADNA_Rep1['log2FC'], K562_SUD1_minP_mRNADNA_Rep2['log2FC'])

# Annotate the plot with the Spearman coefficient
plt.annotate(f'Spearman coefficient: {coef_K562_SUD1_minP:.3f}\nP-value: {p_value_K562_SUD1_minP:.1f}', 
              xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top', 
              bbox=dict(boxstyle='round', facecolor='white', alpha=0))
# Save the plot
plt.savefig('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA/MPRA_plots/K562_SUD1_minP_HexBin.pdf', format='pdf', bbox_inches='tight', transparent=True)

############################## Scaled Up Design 2 #############################
########################## Fig5 HepG2 SUD1 minP HexBin ########################
# Create a hexagonal bin plot
plt.figure(5)
plt.hexbin(HepG2_SUD2_minP_mRNADNA_Rep1['log2FC'], HepG2_SUD2_minP_mRNADNA_Rep2['log2FC'], gridsize=30, cmap='Blues', mincnt=1, bins='log')
# Add a color bar
plt.colorbar(label='Counts in bin')
plt.xlabel('Log2 Fold Change, HepG2 minP Replicate 1')
plt.ylabel('Log2 Fold Change, HepG2 minP Replicate 2')
plt.xlim(-10.5, 10.5)
plt.ylim(-10.5, 10.5)

# Calculate Spearman correlation coefficient
coef_HepG2_SUD2_minP, p_value_HepG2_SUD2_minP = spearmanr(HepG2_SUD2_minP_mRNADNA_Rep1['log2FC'], HepG2_SUD2_minP_mRNADNA_Rep2['log2FC'])

# Annotate the plot with the Spearman coefficient
plt.annotate(f'Spearman coefficient: {coef_HepG2_SUD2_minP:.3f}\nP-value: {p_value_HepG2_SUD2_minP:.1f}', 
              xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top', 
              bbox=dict(boxstyle='round', facecolor='white', alpha=0))
# Save the plot
plt.savefig('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA/MPRA_plots/HepG2_SUD2_minP_HexBin.pdf', format='pdf', bbox_inches='tight', transparent=True)

######################### Fig6 HepG2 SUD2 SV40P HexBin ########################
# Create a hexagonal bin plot
plt.figure(6)
plt.hexbin(HepG2_SUD2_SV40P_mRNADNA_Rep1['log2FC'], HepG2_SUD2_SV40P_mRNADNA_Rep2['log2FC'], gridsize=30, cmap='Blues', mincnt=1, bins='log')
# Add a color bar
plt.colorbar(label='Counts in bin')
plt.xlabel('Log2 Fold Change, HepG2 SV40P Replicate 1')
plt.ylabel('Log2 Fold Change, HepG2 SV40P Replicate 2')
plt.xlim(-10.5, 10.5)
plt.ylim(-10.5, 10.5)

# Calculate Spearman correlation coefficient
coef_HepG2_SUD2_SV40P, p_value_HepG2_SUD2_SV40P = spearmanr(HepG2_SUD2_SV40P_mRNADNA_Rep1['log2FC'], HepG2_SUD2_SV40P_mRNADNA_Rep2['log2FC'])

# Annotate the plot with the Spearman coefficient
plt.annotate(f'Spearman coefficient: {coef_HepG2_SUD2_SV40P:.3f}\nP-value: {p_value_HepG2_SUD2_SV40P:.1f}', 
              xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top', 
              bbox=dict(boxstyle='round', facecolor='white', alpha=0))
# Save the plot
plt.savefig('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA/MPRA_plots/HepG2_SUD2_SV40P_HexBin.pdf', format='pdf', bbox_inches='tight', transparent=True)

########################### Fig7 K562 SUD2 SV40P HexBin #######################
# Create a hexagonal bin plot
plt.figure(7)
plt.hexbin(K562_SUD2_SV40P_mRNADNA_Rep1['log2FC'], K562_SUD2_SV40P_mRNADNA_Rep2['log2FC'], gridsize=30, cmap='Blues', mincnt=1, bins='log')
# Add a color bar
plt.colorbar(label='Counts in bin')
plt.xlabel('Log2 Fold Change, K562 SV40P Replicate 1')
plt.ylabel('Log2 Fold Change, K562 SV40P Replicate 2')
plt.xlim(-10.5, 10.5)
plt.ylim(-10.5, 10.5)

# Calculate Spearman correlation coefficient
coef_K562_SUD2_SV40P, p_value_K562_SUD2_SV40P = spearmanr(K562_SUD2_SV40P_mRNADNA_Rep1['log2FC'], K562_SUD2_SV40P_mRNADNA_Rep2['log2FC'])

# Annotate the plot with the Spearman coefficient
plt.annotate(f'Spearman coefficient: {coef_K562_SUD2_SV40P:.3f}\nP-value: {p_value_K562_SUD2_SV40P:.1f}', 
              xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top', 
              bbox=dict(boxstyle='round', facecolor='white', alpha=0))
# Save the plot
plt.savefig('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA/MPRA_plots/K562_SUD2_SV40P_HexBin.pdf', format='pdf', bbox_inches='tight', transparent=True)

########################### Fig8 K562 SUD2 minP HexBin ########################
# Create a hexagonal bin plot
plt.figure(8)
plt.hexbin(K562_SUD2_minP_mRNADNA_Rep1['log2FC'], K562_SUD2_minP_mRNADNA_Rep2['log2FC'], gridsize=30, cmap='Blues', mincnt=1, bins='log')
# Add a color bar
plt.colorbar(label='Counts in bin')
plt.xlabel('Log2 Fold Change, K562 minP Replicate 1')
plt.ylabel('Log2 Fold Change, K562 minP Replicate 2')
plt.xlim(-10.5, 10.5)
plt.ylim(-10.5, 10.5)

# Calculate Spearman correlation coefficient
coef_K562_SUD2_minP, p_value_K562_SUD2_minP = spearmanr(K562_SUD2_minP_mRNADNA_Rep1['log2FC'], K562_SUD2_minP_mRNADNA_Rep2['log2FC'])

# Annotate the plot with the Spearman coefficient
plt.annotate(f'Spearman coefficient: {coef_K562_SUD2_minP:.3f}\nP-value: {p_value_K562_SUD2_minP:.1f}', 
              xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top', 
              bbox=dict(boxstyle='round', facecolor='white', alpha=0))
# Save the plot
plt.savefig('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA/MPRA_plots/K562_SUD2_minP_HexBin.pdf', format='pdf', bbox_inches='tight', transparent=True)
