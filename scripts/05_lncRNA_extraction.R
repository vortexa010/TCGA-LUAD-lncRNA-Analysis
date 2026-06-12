#=========================================================
# Script: 05_lncRNA_extraction.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Extracts differentially expressed lncRNAs from annotated
# DESeq2 results and selects strong lncRNA candidates for
# downstream WGCNA analysis.
#
# Input:
# TCGA_LUAD_DESeq2_RESULTS.RData
# TCGA_LUAD_FILTERED_DEGs.RData
#
# Output:
# TCGA_LUAD_lncRNA_RESULTS.RData
# Strong_lncRNAs_TCGA_LUAD.csv
#=========================================================

# Load required objects
load("TCGA_LUAD_DESeq2_RESULTS.RData")
load("TCGA_LUAD_FILTERED_DEGs.RData")

#---------------------------------------------------------
# Extract all differentially expressed lncRNAs
#---------------------------------------------------------

lncRNA_deg <- subset(
  deg_annotated,
  gene_type == "lncRNA"
)

#---------------------------------------------------------
# Rank lncRNAs by adjusted p-value and fold change
#---------------------------------------------------------

lncRNA_ranked <- lncRNA_deg[
  order(
    lncRNA_deg$padj,
    -abs(lncRNA_deg$log2FoldChange)
  ),
]

#---------------------------------------------------------
# Select strong lncRNA candidates
# Criteria:
# adjusted p-value < 0.001
# absolute log2 fold change > 2
#---------------------------------------------------------

lncRNA_strong <- subset(
  lncRNA_deg,
  padj < 0.001 &
    abs(log2FoldChange) > 2
)

lncRNA_strong <- lncRNA_strong[
  order(
    lncRNA_strong$padj,
    -abs(lncRNA_strong$log2FoldChange)
  ),
]

#---------------------------------------------------------
# Save output tables
#---------------------------------------------------------

write.csv(
  lncRNA_deg,
  "All_lncRNA_DEGs_TCGA_LUAD.csv",
  row.names = FALSE
)

write.csv(
  lncRNA_ranked,
  "Ranked_lncRNA_DEGs_TCGA_LUAD.csv",
  row.names = FALSE
)

write.csv(
  lncRNA_strong,
  "Strong_lncRNAs_TCGA_LUAD.csv",
  row.names = FALSE
)

#---------------------------------------------------------
# Save R objects
#---------------------------------------------------------

save(
  lncRNA_deg,
  lncRNA_ranked,
  lncRNA_strong,
  file = "TCGA_LUAD_lncRNA_RESULTS.RData"
)

