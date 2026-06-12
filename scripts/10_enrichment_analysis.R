#=========================================================
# Script: 10_enrichment_analysis.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Performs GO Biological Process and KEGG enrichment
# analysis using protein-coding genes positively correlated
# with AC135012.3.
#
# Input:
# TCGA_LUAD_AC1350123_RESULTS.RData
#
# Output:
# GO_BP_results.csv
# KEGG_results.csv
# GO_KEGG_gene_conversion.csv
# TCGA_LUAD_ENRICHMENT_RESULTS.RData
#=========================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Load AC135012.3 characterization results
load("TCGA_LUAD_AC1350123_RESULTS.RData")

# Use the same broader protein-coding correlated gene set
# generated from top100 positive correlations in Script 08
gene_list <- unique(protein_genes$gene_name)

# Remove missing gene names
gene_list <- gene_list[
  !is.na(gene_list)
]

# Convert gene symbols to Entrez IDs
gene_convert <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# GO Biological Process enrichment
ego_bp <- enrichGO(
  gene = gene_convert$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# KEGG pathway enrichment
ekegg <- enrichKEGG(
  gene = gene_convert$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# Convert KEGG IDs to readable gene symbols
ekegg_readable <- setReadable(
  ekegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

# Save enrichment result tables
write.csv(
  as.data.frame(ego_bp),
  "GO_BP_results.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ekegg_readable),
  "KEGG_results.csv",
  row.names = FALSE
)

write.csv(
  gene_convert,
  "GO_KEGG_gene_conversion.csv",
  row.names = FALSE
)

# Save GO dotplot
p_go_dot <- dotplot(
  ego_bp,
  showCategory = 15
)

ggsave(
  "AC135012_GO_BP_dotplot_15terms.png",
  plot = p_go_dot,
  width = 14,
  height = 10,
  dpi = 300
)

# Save GO barplot
p_go_bar <- barplot(
  ego_bp,
  showCategory = 15
)

ggsave(
  "AC135012_GO_BP_barplot_15terms.png",
  plot = p_go_bar,
  width = 14,
  height = 10,
  dpi = 300
)

# Save KEGG dotplot
p_kegg_dot <- dotplot(
  ekegg_readable,
  showCategory = 15
)

ggsave(
  "AC135012_KEGG_dotplot.png",
  plot = p_kegg_dot,
  width = 12,
  height = 8,
  dpi = 300
)

# Save R objects
save(
  protein_genes,
  gene_list,
  gene_convert,
  ego_bp,
  ekegg,
  ekegg_readable,
  file = "TCGA_LUAD_ENRICHMENT_RESULTS.RData"
)
