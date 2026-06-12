#=========================================================
# Script: 12_cox_regression.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Performs univariate Cox proportional hazards regression
# for AC135012.3 expression using TCGA-LUAD survival data.
#
# Input:
# TCGA_LUAD_SURVIVAL_RESULTS.RData
#
# Output:
# TCGA_LUAD_COX_RESULTS.RData
# Cox_regression_summary.txt
#=========================================================

library(survival)

#---------------------------------------------------------
# Load survival data
#---------------------------------------------------------

load("TCGA_LUAD_SURVIVAL_RESULTS.RData")

#---------------------------------------------------------
# Univariate Cox regression
#---------------------------------------------------------

cox2 <- coxph(
  Surv(OS_time, status) ~ expression,
  data = surv_data2
)

cox_summary <- summary(cox2)

cox_summary

#---------------------------------------------------------
# Save Cox summary as text
#---------------------------------------------------------

sink("Cox_regression_summary.txt")
print(cox_summary)
sink()

#---------------------------------------------------------
# Save R objects
#---------------------------------------------------------

save(
  cox2,
  cox_summary,
  file = "TCGA_LUAD_COX_RESULTS.RData"
)
