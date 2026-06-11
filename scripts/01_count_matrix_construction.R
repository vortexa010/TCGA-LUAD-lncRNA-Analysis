#=========================================================
# Script: 01_count_matrix_construction.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Downloads and processes TCGA-LUAD RNA-seq count data,
# constructs the count matrix, prepares sample metadata,
# and saves clean objects for downstream analysis.
#
# Input:
#   - GDC HTSeq count files
#   - Clinical/sample metadata
#
# Output:
#   - Count matrix
#   - Sample metadata
#   - Gene annotation
#   - TCGA_LUAD_DESeq2_READY.RData
#=========================================================
