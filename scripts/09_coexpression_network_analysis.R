#=========================================================
# Script: 09_coexpression_network_analysis.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Constructs the AC135012.3 co-expression network using
# the network_genes object generated during AC135012.3
# characterization.
#
# Input:
# TCGA_LUAD_AC1350123_RESULTS.RData
#
# Output:
# AC135012_network.csv
# Network_gene_table.csv
# TCGA_LUAD_NETWORK_RESULTS.RData
#=========================================================

# Load AC135012.3 characterization results
load("TCGA_LUAD_AC1350123_RESULTS.RData")

# Use the same network genes generated in Script 08
network_table <- data.frame(
  Source = "AC135012.3",
  Target = network_genes$gene_name,
  Weight = network_genes$correlation
)

# Sort network by correlation strength
network_table <- network_table[
  order(-network_table$Weight),
]

# Save Cytoscape-ready network file
write.csv(
  network_table,
  "AC135012_network.csv",
  row.names = FALSE
)

# Save the gene table used for the network
write.csv(
  network_genes,
  "Network_gene_table.csv",
  row.names = FALSE
)

# Save R objects
save(
  network_genes,
  network_table,
  file = "TCGA_LUAD_NETWORK_RESULTS.RData"
)
