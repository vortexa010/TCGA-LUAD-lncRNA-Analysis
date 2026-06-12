#=========================================================
# Script: 02_DESeq2_analysis.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Performs differential expression analysis using DESeq2
# on the cleaned TCGA-LUAD count matrix.
#
# Input:
# TCGA_LUAD_DESeq2_READY_CLEAN.RData
#
# Output:
# TCGA_LUAD_DESeq2_RESULTS.RData
# TCGA_LUAD_Annotated_DEGs.csv
#=========================================================

library(DESeq2)

# Load clean count matrix and sample metadata
load("TCGA_LUAD_DESeq2_READY_CLEAN.RData")

# Check condition groups
table(coldata$condition)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = tcga_count1,
  colData = coldata,
  design = ~ condition
)

# Run DESeq2 differential expression analysis
# DESeq2 uses the Wald test by default
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Summary of DESeq2 results
summary(res)

# Convert significant results to data frame
sig <- res[!is.na(res$padj) & res$padj < 0.05, ]

deg <- as.data.frame(sig)
deg$gene_id <- rownames(deg)

# Prepare gene annotation table
gene_info2 <- gene_info
gene_info2$gene_id <- sub("\\..*", "", gene_info2$gene_id)

# Remove duplicated gene annotation rows
gene_info2_unique <- gene_info2[
  !duplicated(gene_info2$gene_id),
]

# Merge DESeq2 results with gene annotation
deg_annotated <- merge(
  deg,
  gene_info2_unique,
  by = "gene_id"
)

# Check gene biotypes
table(deg_annotated$gene_type)

# Save annotated DEG table
write.csv(
  deg_annotated,
  "TCGA_LUAD_Annotated_DEGs.csv",
  row.names = FALSE
)

# Save DESeq2 output objects
save(
  dds,
  res,
  deg,
  deg_annotated,
  gene_info2_unique,
  coldata,
  file = "TCGA_LUAD_DESeq2_RESULTS.RData"
)
