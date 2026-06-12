#=========================================================
# Script: 13_validation_GEPIA2.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Documents external validation of AC135012.3 using GEPIA2.
# GEPIA2 analyses were performed manually through the web server.
#
# Validation performed:
# 1. Tumor vs normal expression comparison
# 2. Overall survival analysis
# 3. Correlation analysis with endothelial/vascular genes
#
# Output:
# GEPIA2 validation figures saved in figures/validation/
#=========================================================

#---------------------------------------------------------
# GEPIA2 validation summary
#---------------------------------------------------------

validation_summary <- data.frame(
  Analysis = c(
    "Expression validation",
    "Survival validation",
    "Correlation with TEK",
    "Correlation with SLC6A4",
    "Correlation with TMEM100",
    "Correlation with EDNRB",
    "Correlation with CD300LG"
  ),
  Result = c(
    "AC135012.3 showed lower expression in LUAD tumor samples compared with normal lung samples",
    "High AC135012.3 expression was associated with better overall survival",
    "Strong positive correlation",
    "Strong positive correlation",
    "Strong positive correlation",
    "Strong positive correlation",
    "Strong positive correlation"
  ),
  Value = c(
    "Tumor < Normal",
    "Log-rank p = 0.0075; HR = 0.64",
    "R = 0.78",
    "R = 0.76",
    "R = 0.71",
    "R = 0.67",
    "R = 0.64"
  )
)

# Save validation summary table
write.csv(
  validation_summary,
  "GEPIA2_validation_summary.csv",
  row.names = FALSE
)

# Save R object
save(
  validation_summary,
  file = "TCGA_LUAD_GEPIA2_VALIDATION.RData"
)
