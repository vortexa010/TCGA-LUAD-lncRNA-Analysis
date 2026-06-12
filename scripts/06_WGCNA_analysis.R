#=========================================================
# Script: 06_WGCNA_analysis.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Performs Weighted Gene Co-expression Network Analysis
# using strong differentially expressed lncRNAs from
# TCGA-LUAD. The input expression matrix is derived from
# variance-stabilized DESeq2 expression values.
#
# Input:
# TCGA_LUAD_DESeq2_RESULTS.RData
# TCGA_LUAD_lncRNA_RESULTS.RData
#
# Output:
# TCGA_LUAD_WGCNA_RESULTS.RData
# WGCNA_sample_clustering.png
#=========================================================

library(DESeq2)
library(WGCNA)

set.seed(123)

# Allow multi-threading
allowWGCNAThreads()

#---------------------------------------------------------
# Load required data
#---------------------------------------------------------

load("TCGA_LUAD_DESeq2_RESULTS.RData")
load("TCGA_LUAD_lncRNA_RESULTS.RData")

#---------------------------------------------------------
# Prepare VST expression matrix
#---------------------------------------------------------

vsd <- vst(dds, blind = FALSE)

expr_mat <- assay(vsd)

# Select strong lncRNAs
lncRNA_ids <- lncRNA_strong$gene_id

lnc_expr <- expr_mat[
  rownames(expr_mat) %in% lncRNA_ids,
]

# WGCNA requires samples as rows and genes as columns
data_Expr <- t(lnc_expr)

dim(data_Expr)

#---------------------------------------------------------
# Quality control
#---------------------------------------------------------

gsg <- goodSamplesGenes(
  data_Expr,
  verbose = 3
)

gsg$allOK

#---------------------------------------------------------
# Sample clustering
#---------------------------------------------------------

sample_Tree <- hclust(
  dist(data_Expr),
  method = "average"
)

png(
  "WGCNA_sample_clustering.png",
  width = 1800,
  height = 1200,
  res = 200
)

plot(
  sample_Tree,
  main = "Sample clustering of LUAD samples",
  labels = FALSE,
  hang = 0.04
)

dev.off()

#---------------------------------------------------------
# Soft-thresholding power selection
#---------------------------------------------------------

powers <- c(1:20)

sft <- pickSoftThreshold(
  data_Expr,
  powerVector = powers,
  verbose = 5
)

# Based on soft-thresholding result, power 8 was selected
softPower <- 8

#---------------------------------------------------------
# Network construction and module detection
#---------------------------------------------------------

net <- blockwiseModules(
  data_Expr,
  power = softPower,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

# Convert numeric module labels to colors
moduleColors <- labels2colors(net$colors)

table(moduleColors)

#---------------------------------------------------------
# Calculate module eigengenes
#---------------------------------------------------------

MEs <- moduleEigengenes(
  data_Expr,
  moduleColors
)$eigengenes

#---------------------------------------------------------
# Correlate modules with tumor/normal trait
#---------------------------------------------------------

traitData <- data.frame(
  Tumor = ifelse(coldata$condition == "Tumor", 1, 0)
)

rownames(traitData) <- rownames(coldata)

moduleTraitCor <- cor(
  MEs,
  traitData,
  use = "p"
)

moduleTraitPvalue <- corPvalueStudent(
  moduleTraitCor,
  nSamples = nrow(data_Expr)
)

moduleTraitCor
moduleTraitPvalue

#---------------------------------------------------------
# Save WGCNA results
#---------------------------------------------------------

save(
  data_Expr,
  sft,
  softPower,
  net,
  moduleColors,
  MEs,
  moduleTraitCor,
  moduleTraitPvalue,
  file = "TCGA_LUAD_WGCNA_RESULTS.RData"
)
