#=========================================================
# Script: 03_visualization.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Performs visualization of DESeq2 results including:
# 1. PCA plot
# 2. Volcano plot
# 3. Heatmap of top DEGs
#
# Input:
# TCGA_LUAD_DESeq2_RESULTS.RData
#
# Output:
# PCA plot
# Volcano plot
# Heatmap
#=========================================================

library(DESeq2)
library(ggplot2)
library(pheatmap)

# Load DESeq2 results
load("TCGA_LUAD_DESeq2_RESULTS.RData")

#---------------------------------------------------------
# Variance Stabilizing Transformation
#---------------------------------------------------------

vsd <- vst(dds, blind = FALSE)

#---------------------------------------------------------
# PCA Plot
#---------------------------------------------------------

pcaData <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pcaData, "percentVar")
)

ggplot(
  pcaData,
  aes(
    PC1,
    PC2,
    color = condition
  )
) +
  geom_point(size = 3) +
  xlab(
    paste0(
      "PC1: ",
      percentVar[1],
      "% variance"
    )
  ) +
  ylab(
    paste0(
      "PC2: ",
      percentVar[2],
      "% variance"
    )
  ) +
  theme_minimal()

#---------------------------------------------------------
# Volcano Plot
#---------------------------------------------------------

resDF <- as.data.frame(res)

resDF$threshold <- "Not Significant"

resDF$threshold[
  resDF$padj < 0.05 &
    resDF$log2FoldChange > 1
] <- "Upregulated"

resDF$threshold[
  resDF$padj < 0.05 &
    resDF$log2FoldChange < -1
] <- "Downregulated"

ggplot(
  resDF,
  aes(
    log2FoldChange,
    -log10(padj)
  )
) +
  geom_point(
    aes(color = threshold),
    alpha = 0.6,
    size = 1.5
  ) +
  scale_color_manual(
    values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not Significant" = "grey"
    )
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  )

#---------------------------------------------------------
# Heatmap of Top 50 DEGs
#---------------------------------------------------------

resOrdered <- res[
  order(res$padj),
]

topGenes <- head(
  rownames(resOrdered),
  50
)

vsd_mat <- assay(vsd)[
  topGenes,
]

annotation_col <- data.frame(
  Condition = coldata$condition
)

rownames(annotation_col) <- rownames(coldata)

pheatmap(
  vsd_mat,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = FALSE,
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)
