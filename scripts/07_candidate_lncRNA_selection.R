#=========================================================
# Script: 07_candidate_lncRNA_selection.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Extracts hub lncRNAs from the Blue WGCNA module and
# selects final candidate hub lncRNAs for downstream
# characterization.
#
# Input:
# TCGA_LUAD_WGCNA_RESULTS.RData
# TCGA_LUAD_DESeq2_RESULTS.RData
#
# Output:
# TCGA_LUAD_CANDIDATE_lncRNAs.RData
# Blue_module_all_hub_lncRNAs.csv
# Final_Top10_candidate_hub_lncRNAs_TCGA_LUAD.csv
#=========================================================

#---------------------------------------------------------
# Load required objects
#---------------------------------------------------------

load("TCGA_LUAD_WGCNA_RESULTS.RData")
load("TCGA_LUAD_DESeq2_RESULTS.RData")

set.seed(123)

#---------------------------------------------------------
# Extract Blue module lncRNAs
#---------------------------------------------------------

blue_lncRNAs <- colnames(data_Expr)[
  moduleColors == "blue"
]

blue_module_df <- data.frame(
  gene_id = blue_lncRNAs
)

# Annotate Blue module lncRNAs
blue_module_df <- merge(
  blue_module_df,
  gene_info2_unique,
  by = "gene_id"
)

#---------------------------------------------------------
# Calculate module membership/kME for Blue module
#---------------------------------------------------------

blue_expr <- data_Expr[, blue_lncRNAs]

blue_ME <- MEs$MEblue

blue_kME <- cor(
  blue_expr,
  blue_ME,
  use = "p"
)

blue_kME_df <- data.frame(
  gene_id = colnames(blue_expr),
  kME = as.numeric(blue_kME)
)

blue_hub_lncRNAs <- merge(
  blue_kME_df,
  blue_module_df,
  by = "gene_id"
)

# Add DESeq2 statistics
blue_hub_lncRNAs <- merge(
  blue_hub_lncRNAs,
  deg_annotated[, c(
    "gene_id",
    "log2FoldChange",
    "padj"
  )],
  by = "gene_id"
)

# Rank hub lncRNAs by absolute kME
blue_hub_lncRNAs <- blue_hub_lncRNAs[
  order(-abs(blue_hub_lncRNAs$kME)),
]

#---------------------------------------------------------
# Separate positive and negative hub lncRNAs
#---------------------------------------------------------

blue_positive_hubs <- subset(
  blue_hub_lncRNAs,
  kME > 0
)

blue_negative_hubs <- subset(
  blue_hub_lncRNAs,
  kME < 0
)

#---------------------------------------------------------
# Select final candidate lncRNAs
#These lncRNAs were selected for downstream analysis
# based on differential expression, WGCNA hub status,
# module membership, and biological relevance after
# manual review of the candidate list.
#---------------------------------------------------------

final_lncRNAs <- c(
  "FENDRR",
  "LINC00968",
  "LINC02016",
  "HID1-AS1",
  "LANCL1-AS1",
  "AL606469.1",
  "AC135012.3",
  "AL590226.1",
  "LINC01996",
  "AC093110.1"
)

candidate_lncRNAs <- subset(
  blue_positive_hubs,
  gene_name %in% final_lncRNAs
)

candidate_lncRNAs <- candidate_lncRNAs[
  order(-candidate_lncRNAs$kME),
]

#---------------------------------------------------------
# Save output tables
#---------------------------------------------------------

write.csv(
  blue_hub_lncRNAs,
  "Blue_module_all_hub_lncRNAs.csv",
  row.names = FALSE
)

write.csv(
  blue_positive_hubs,
  "Blue_module_positive_hub_lncRNAs.csv",
  row.names = FALSE
)

write.csv(
  blue_negative_hubs,
  "Blue_module_negative_hub_lncRNAs.csv",
  row.names = FALSE
)

write.csv(
  candidate_lncRNAs,
  "Final_Top10_candidate_hub_lncRNAs_TCGA_LUAD.csv",
  row.names = FALSE
)

#---------------------------------------------------------
# Save R objects
#---------------------------------------------------------

save(
  blue_module_df,
  blue_hub_lncRNAs,
  blue_positive_hubs,
  blue_negative_hubs,
  candidate_lncRNAs,
  file = "TCGA_LUAD_CANDIDATE_lncRNAs.RData"
)
