#=========================================================
# Script: 08_AC1350123_characterization.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Characterizes the candidate hub lncRNA AC135012.3 by
# examining its expression pattern and identifying
# co-expressed genes in TCGA-LUAD.
#
# Input:
# TCGA_LUAD_DESeq2_RESULTS.RData
# TCGA_LUAD_CANDIDATE_lncRNAs.RData
#
# Output:
# AC1350123_expression.csv
# AC1350123_correlated_genes.csv
# TCGA_LUAD_AC1350123_RESULTS.RData
#=========================================================

library(DESeq2)

#---------------------------------------------------------
# Load required objects
#---------------------------------------------------------

load("TCGA_LUAD_DESeq2_RESULTS.RData")
load("TCGA_LUAD_CANDIDATE_lncRNAs.RData")

#---------------------------------------------------------
# Generate variance-stabilized expression matrix
#---------------------------------------------------------

vsd <- vst(dds, blind = FALSE)

expr_mat <- assay(vsd)

#---------------------------------------------------------
# Extract AC135012.3 expression
#---------------------------------------------------------

target_lncRNA <- "ENSG00000268505"

ac_expr <- expr_mat[
  target_lncRNA,
]

#---------------------------------------------------------
# Create expression table
#---------------------------------------------------------

ac_expression_table <- data.frame(
  Sample = colnames(expr_mat),
  Expression = as.numeric(ac_expr),
  Condition = coldata$condition
)

#---------------------------------------------------------
# Correlation analysis
#---------------------------------------------------------

cor_values <- apply(
  expr_mat,
  1,
  function(x)
    cor(
      x,
      ac_expr,
      method = "pearson"
    )
)

# Remove self-correlation
cor_values <- cor_values[
  names(cor_values) != target_lncRNA
]

#---------------------------------------------------------
# Rank correlated genes
#---------------------------------------------------------

top_pos <- sort(
  cor_values,
  decreasing = TRUE
)

top_neg <- sort(
  cor_values,
  decreasing = FALSE
)

#---------------------------------------------------------
# Top 100 positively correlated genes
#---------------------------------------------------------

top100 <- data.frame(
  gene_id = names(top_pos)[1:100],
  correlation = top_pos[1:100]
)

top100 <- merge(
  top100,
  gene_info2_unique,
  by = "gene_id",
  all.x = TRUE
)

#---------------------------------------------------------
# Top 100 negatively correlated genes
#---------------------------------------------------------

top100_neg <- data.frame(
  gene_id = names(top_neg)[1:100],
  correlation = top_neg[1:100]
)

top100_neg <- merge(
  top100_neg,
  gene_info2_unique,
  by = "gene_id",
  all.x = TRUE
)

#---------------------------------------------------------
# Extract protein-coding correlated genes
#---------------------------------------------------------

protein_genes <- subset(
  top100,
  gene_type == "protein_coding"
)

network_genes <- head(
  protein_genes,
  20
)

#---------------------------------------------------------
# Save outputs
#---------------------------------------------------------

write.csv(
  ac_expression_table,
  "AC1350123_expression.csv",
  row.names = FALSE
)

write.csv(
  top100,
  "AC1350123_top100_positive_correlated_genes.csv",
  row.names = FALSE
)

write.csv(
  top100_neg,
  "AC1350123_top100_negative_correlated_genes.csv",
  row.names = FALSE
)

write.csv(
  network_genes,
  "AC1350123_network_genes.csv",
  row.names = FALSE
)

#---------------------------------------------------------
# Save R objects
#---------------------------------------------------------

save(
  ac_expr,
  ac_expression_table,
  cor_values,
  top100,
  top100_neg,
  protein_genes,
  network_genes,
  file = "TCGA_LUAD_AC1350123_RESULTS.RData"
)

