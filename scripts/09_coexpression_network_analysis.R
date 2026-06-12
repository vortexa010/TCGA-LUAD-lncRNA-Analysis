#=========================================================
# Script: 09_coexpression_network_analysis.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Constructs the AC135012.3 co-expression network
# using the top correlated protein-coding genes.
#
# Input:
# TCGA_LUAD_AC1350123_RESULTS.RData
#
# Output:
# AC135012_network.csv
# Network_gene_table.csv
#=========================================================

#---------------------------------------------------------
# Load required objects
#---------------------------------------------------------

load("TCGA_LUAD_AC1350123_RESULTS.RData")

#---------------------------------------------------------
# Select top correlated genes
#---------------------------------------------------------

top20 <- head(
  protein_genes,
  20
)

#---------------------------------------------------------
# Create network table
#---------------------------------------------------------

network_table <- data.frame(
  Source = "AC135012.3",
  Target = top20$gene_name,
  Weight = top20$correlation
)

#---------------------------------------------------------
# Sort by correlation
#---------------------------------------------------------

network_table <- network_table[
  order(
    -network_table$Weight
  ),
]

#---------------------------------------------------------
# Save network files
#---------------------------------------------------------

write.csv(
  network_table,
  "AC135012_network.csv",
  row.names = FALSE
)

write.csv(
  top20,
  "Network_gene_table.csv",
  row.names = FALSE
)

#---------------------------------------------------------
# Save R objects
#---------------------------------------------------------

save(
  top20,
  network_table,
  file = "TCGA_LUAD_NETWORK_RESULTS.RData"
)
