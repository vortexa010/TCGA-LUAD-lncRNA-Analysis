#=========================================================
# Script: 01_count_matrix_construction.R
# Author: Vaishnavi Mangam
# Project: TCGA-LUAD-lncRNA-Analysis
#
# Description:
# Constructs a clean TCGA-LUAD count matrix from GDC HTSeq
# count files and prepares sample metadata for DESeq2.
#
# Output:
# TCGA_LUAD_DESeq2_READY_CLEAN.RData
#=========================================================

# Set GDC download folder path
TCGA <- "C:\\Users\\hp\\OneDrive\\Desktop\\gdc-client_2.3_Windows_x64"

# Find all HTSeq count TSV files
files <- list.files(
  path = TCGA,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

length(files)
head(files)

# Read first count file to extract gene annotation
data1 <- read.delim(files[1], skip = 1)

# Extract gene information
gene_info <- data1[, c("gene_id", "gene_name", "gene_type")]

# Remove first 4 technical rows
gene_info <- gene_info[-c(1:4), ]

# Build count matrix
count_list <- list()

for (i in seq_along(files)) {
  count_file <- read.delim(files[i], skip = 1)
  count_values <- count_file$unstranded
  count_values <- count_values[-c(1:4)]
  count_list[[i]] <- count_values
  print(i)
}

count_matrix_full <- do.call(cbind, count_list)

count_matrix_full <- data.frame(
  gene_id = gene_info$gene_id,
  count_matrix_full
)

# Load sample sheet
sample_sheet <- read.delim(
  file.path(TCGA, "gdc_sample_sheet.2026-05-26.tsv")
)

# Match downloaded file names with sample sheet
file_names <- basename(files)

match_idx <- match(
  file_names,
  sample_sheet$File.Name
)

sample_sheet12 <- sample_sheet[match_idx, ]

# Add sample IDs as column names
colnames(count_matrix_full)[-1] <- sample_sheet12$Sample.ID

# Remove Ensembl version numbers
count_matrix_full$gene_id <- sub("\\..*", "", count_matrix_full$gene_id)
gene_info$gene_id <- sub("\\..*", "", gene_info$gene_id)

# Set gene IDs as rownames
rownames(count_matrix_full) <- count_matrix_full$gene_id

# Remove gene_id column temporarily
count_only <- count_matrix_full[, -1]

# Ensure numeric count matrix
count_only <- as.data.frame(
  lapply(count_only, as.numeric)
)

# Sum duplicated gene IDs safely
count_matrix_sum <- rowsum(
  count_only,
  group = rownames(count_matrix_full)
)

# Bring gene_id back as a column
count_matrix_sum <- data.frame(
  gene_id = rownames(count_matrix_sum),
  count_matrix_sum
)

# Prepare final count matrix for DESeq2
tcga_count1 <- count_matrix_sum
rownames(tcga_count1) <- tcga_count1$gene_id
tcga_count1 <- tcga_count1[, -1]

# Remove duplicated sample IDs from sample sheet
sample_sheet_unique <- sample_sheet12[
  !is.na(sample_sheet12$Sample.ID) &
    sample_sheet12$Sample.ID != "" &
    !duplicated(sample_sheet12$Sample.ID),
]

# Create sample metadata
coldata <- data.frame(
  row.names = sample_sheet_unique$Sample.ID,
  condition = sample_sheet_unique$Tissue.Type
)

# Check sample alignment
all(rownames(coldata) == colnames(tcga_count1))

# Check tumor/normal count
table(coldata$condition)

# Save clean objects for DESeq2
save(
  tcga_count1,
  coldata,
  gene_info,
  sample_sheet_unique,
  file = "TCGA_LUAD_DESeq2_READY_CLEAN.RData"
)
