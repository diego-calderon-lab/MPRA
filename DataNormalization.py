# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 2024

@author: angela_gao
"""

import pandas as pd

##################################################################################################
####################################### DATA INPUT ###############################################
##################################################################################################

HepG2_SUD1_minP_mRNA_Rep1 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831757_HepG2_ScaleUpDesign1_minP_mRNA_Rep1.counts.txt', sep='\t')
HepG2_SUD1_minP_mRNA_Rep2 = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831758_HepG2_ScaleUpDesign1_minP_mRNA_Rep2.counts.txt', sep='\t')
SUD1_minP_DNA = pd.read_csv('/Users/angela_gao/Dropbox/Angela Gao - UCSF Tetrad/Calderon Lab Rotation/MPRA_Project/MPRA_dataset/GSE71279_RAW/GSM1831773_ScaleUpDesign1_minP_Plasmid.counts.txt', sep='\t')
