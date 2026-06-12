#=========================================================
# Script: 04_DEG_filtering.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Filters annotated DESeq2 results to identify significant
# differentially expressed genes, upregulated genes,
# downregulated genes, and major RNA biotypes.
#
# Input:
# TCGA_LUAD_DESeq2_RESULTS.RData
#
# Output:
# TCGA_LUAD_FILTERED_DEGs.RData
# Upregulated_genes.csv
# Downregulated_genes.csv
# Protein_coding_DEGs.csv
# lncRNA_DEGs.csv
#=========================================================

# Load DESeq2 annotated results
load("TCGA_LUAD_DESeq2_RESULTS.RData")

#---------------------------------------------------------
# Filter biologically meaningful DEGs
# Criteria:
# adjusted p-value < 0.05
# absolute log2 fold change > 1
#---------------------------------------------------------

resSig <- subset(
  deg_annotated,
  padj < 0.05 &
    abs(log2FoldChange) > 1
)

#---------------------------------------------------------
# Separate upregulated and downregulated genes
#---------------------------------------------------------

upGenes <- subset(
  resSig,
  log2FoldChange > 1
)

downGenes <- subset(
  resSig,
  log2FoldChange < -1
)

#---------------------------------------------------------
# Separate major gene biotypes
#---------------------------------------------------------

protein_coding_sig <- subset(
  resSig,
  gene_type == "protein_coding"
)

lncRNA_sig <- subset(
  resSig,
  gene_type == "lncRNA"
)

other_rna <- subset(
  resSig,
  gene_type != "protein_coding" &
    gene_type != "lncRNA"
)

#---------------------------------------------------------
# Extract top upregulated/downregulated genes
#---------------------------------------------------------

top_up <- head(
  protein_coding_sig[
    order(-protein_coding_sig$log2FoldChange),
  ],
  20
)

top_down <- head(
  protein_coding_sig[
    order(protein_coding_sig$log2FoldChange),
  ],
  20
)

top_uplnc <- head(
  lncRNA_sig[
    order(-lncRNA_sig$log2FoldChange),
  ],
  20
)

top_downlnc <- head(
  lncRNA_sig[
    order(lncRNA_sig$log2FoldChange),
  ],
  20
)

#---------------------------------------------------------
# Save tables
#---------------------------------------------------------

write.csv(
  upGenes,
  "Upregulated_genes_TCGA_LUAD.csv",
  row.names = FALSE
)

write.csv(
  downGenes,
  "Downregulated_genes_TCGA_LUAD.csv",
  row.names = FALSE
)

write.csv(
  protein_coding_sig,
  "Protein_coding_DEGs_TCGA_LUAD.csv",
  row.names = FALSE
)

write.csv(
  lncRNA_sig,
  "lncRNA_DEGs_TCGA_LUAD.csv",
  row.names = FALSE
)

write.csv(
  top_uplnc,
  "Top_20_upregulated_lncRNAs_TCGA_LUAD.csv",
  row.names = FALSE
)

write.csv(
  top_downlnc,
  "Top_20_downregulated_lncRNAs_TCGA_LUAD.csv",
  row.names = FALSE
)

#---------------------------------------------------------
# Save R objects
#---------------------------------------------------------

save(
  resSig,
  upGenes,
  downGenes,
  protein_coding_sig,
  lncRNA_sig,
  other_rna,
  top_up,
  top_down,
  top_uplnc,
  top_downlnc,
  file = "TCGA_LUAD_FILTERED_DEGs.RData"
)

