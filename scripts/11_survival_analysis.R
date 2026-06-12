#=========================================================
# Script: 11_survival_analysis.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Performs Kaplan-Meier survival analysis for AC135012.3
# using TCGA-LUAD clinical data and expression data.
#
# Input:
# TCGA_LUAD_DESeq2_RESULTS.RData
# clinical.cart.2026-05-26.json
#
# Output:
# TCGA_LUAD_SURVIVAL_RESULTS.RData
#=========================================================

library(DESeq2)
library(jsonlite)
library(survival)
library(survminer)

# Load expression data
load("TCGA_LUAD_DESeq2_RESULTS.RData")

# Set TCGA folder path
TCGA <- "C:\\Users\\hp\\OneDrive\\Desktop\\gdc-client_2.3_Windows_x64"

# Load clinical file
clinical <- fromJSON(
  file.path(TCGA, "clinical.cart.2026-05-26.json")
)

# Create VST expression matrix
vsd <- vst(dds, blind = FALSE)
expr_mat <- assay(vsd)

# Extract patient IDs from TCGA sample barcodes
patient_ids <- substr(
  colnames(expr_mat),
  1,
  12
)

# Extract AC135012.3 expression
ac_surv <- data.frame(
  patient_id = patient_ids,
  expression = as.numeric(
    expr_mat["ENSG00000268505", ]
  )
)

# Extract follow-up days from diagnoses
followup_days <- sapply(
  clinical$diagnoses,
  function(x) {
    if (length(x$days_to_last_follow_up) == 0) {
      return(NA)
    }
    max(x$days_to_last_follow_up, na.rm = TRUE)
  }
)

# Build survival clinical table
surv_clin2 <- data.frame(
  patient_id = clinical$submitter_id,
  vital_status = clinical$demographic$vital_status,
  days_to_death = clinical$demographic$days_to_death,
  days_to_last_follow_up = followup_days
)

# Create overall survival time
surv_clin2$OS_time <- ifelse(
  surv_clin2$vital_status == "Dead",
  surv_clin2$days_to_death,
  surv_clin2$days_to_last_follow_up
)

# Create survival status
surv_clin2$status <- ifelse(
  surv_clin2$vital_status == "Dead",
  1,
  0
)

# Merge expression with clinical survival data
surv_data2 <- merge(
  ac_surv,
  surv_clin2,
  by = "patient_id"
)

# Remove samples with missing survival time
surv_data2 <- surv_data2[
  !is.na(surv_data2$OS_time),
]

# Divide patients into high and low expression groups
surv_data2$group <- ifelse(
  surv_data2$expression >
    median(surv_data2$expression, na.rm = TRUE),
  "High",
  "Low"
)

table(surv_data2$group)

# Kaplan-Meier survival analysis
fit2 <- survfit(
  Surv(OS_time, status) ~ group,
  data = surv_data2
)

km_plot <- ggsurvplot(
  fit2,
  data = surv_data2,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  legend.title = "AC135012.3",
  legend.labs = c("High", "Low")
)

km_plot

# Save results
save(
  surv_clin2,
  surv_data2,
  fit2,
  km_plot,
  file = "TCGA_LUAD_SURVIVAL_RESULTS.RData"
)
