#=========================================================
# Script: 10_enrichment_analysis.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Performs GO Biological Process and KEGG pathway
# enrichment analysis using protein-coding genes
# co-expressed with AC135012.3.
#
# Input:
# TCGA_LUAD_NETWORK_RESULTS.RData
#
# Output:
# GO_results.csv
# KEGG_results.csv
#=========================================================

library(clusterProfiler)
library(org.Hs.eg.db)

#---------------------------------------------------------
# Load network genes
#---------------------------------------------------------

load("TCGA_LUAD_NETWORK_RESULTS.RData")

genes <- unique(top20$gene_name)

#---------------------------------------------------------
# Convert gene symbols to Entrez IDs
#---------------------------------------------------------

gene_convert <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

#---------------------------------------------------------
# GO enrichment
#---------------------------------------------------------

ego_bp <- enrichGO(
  gene = gene_convert$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH"
)

#---------------------------------------------------------
# KEGG enrichment
#---------------------------------------------------------

ekegg <- enrichKEGG(
  gene = gene_convert$ENTREZID,
  organism = "hsa"
)

#---------------------------------------------------------
# Save results
#---------------------------------------------------------

write.csv(
  as.data.frame(ego_bp),
  "GO_BP_results.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ekegg),
  "KEGG_results.csv",
  row.names = FALSE
)

save(
  ego_bp,
  ekegg,
  gene_convert,
  file = "TCGA_LUAD_ENRICHMENT_RESULTS.RData"
)

# Optional plots
dotplot(ego_bp)

barplot(ego_bp)

dotplot(ekegg)
